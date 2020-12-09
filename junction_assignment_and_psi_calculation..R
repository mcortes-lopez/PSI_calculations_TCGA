library(data.table)
new_genes<-fread("data/genes_new.csv",data.table = F)


library(tidyr)
library(dplyr)


# Get searchable coordinates
gene_coord<-new_genes %>% 
  dplyr::select(starts_with("chr"), strand) %>% 
  dplyr::select(-chr_size) %>% 
  dplyr::rename("seqnames" = "chr_name","start" = "chr_start", "end" = "chr_end") %>% 
  GenomicRanges::GRanges() 


# The coordinates are in the hg19 version on the genome. Normally we use the version hg38
# We will proceed to 1) convert the coordinates to hg38 


library(rtracklayer)
path = "data/hg19ToHg38.over.chain"
ch = import.chain(path)
ch

new_exons<-liftOver(gene_coord, ch) %>% 
  unlist()


# Now we need to get the upstream and downstream junctions.
# For this we need to overlap the exon with the introns



library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg38_exons<-exons(txdb, columns= c( "TXID", "TXNAME","EXONID"))
hg38_introns<-intronsByTranscript(txdb, use.names=T)

ovp<-findOverlapPairs(new_exons, hg38_exons, type="any")

# The start of the exons is off by 1 nt 

proper_ovp<-ovp[(width(first(ovp)) - width(second(ovp)))==1]

library(IRanges)


# Getting upstream and downstream exons 
paired_introns<-lapply(seq_along(proper_ovp), function(exon_idx){
  
  introns_un<-unlist(hg38_introns)
  
  upstream_int<-precede(second(proper_ovp[exon_idx]), introns_un, select="all")
  downstream_int<-follow(second(proper_ovp[exon_idx]), introns_un, select="all")

  tx_names_up<-names(introns_un[subjectHits(upstream_int)]) 
  tx_names_down<-names(introns_un[subjectHits(downstream_int)])
  
  match_tx<- intersect(tx_names_up, tx_names_down) 
  
  tx_names_exon<-second(proper_ovp)[exon_idx]$TXNAME


  tx_matched<-intersect(unlist(tx_names_exon), match_tx)

  
  matched_up_introns<- introns_un[subjectHits(upstream_int)[tx_names_up %in% tx_matched]]

  matched_down_introns<- introns_un[subjectHits(downstream_int)[tx_names_down %in% tx_matched]]
  
  
  
  return(list("upstream"=unique(matched_up_introns), 
              "downstream"=unique(matched_down_introns)))
  
})



# Generating the junctions table 
# If "+" -> upstream_ss = start(upstream)  and  downstream_ss = end(downstream)
# If "-" -> upstream_ss = end(downstream)  and  downstream_ss = start(upstream)

library(snapcount)


# Calculation of PSI based on the SNAPTRON manual and using the recently updated data from TCGA (V2)
psi_query_calc<-function(ljx, rjx, skipjx, ex_str){
  print(list(ljx, rjx, skipjx, ex_str))

  lq <- QueryBuilder(compilation="tcgav2", regions=ljx)
      
  lq <- set_row_filters(lq, strand == {{ ex_str }} )
  lq <- set_coordinate_modifier(lq, Coordinates$Exact)
  #right inclusion query
  rq <- QueryBuilder(compilation="tcgav2", regions=rjx)
  rq <- set_row_filters(rq, strand == {{ ex_str}})
  rq <- set_coordinate_modifier(rq, Coordinates$Exact)
  #exclusion query
  ex <- QueryBuilder(compilation="tcgav2", regions=skipjx)
  ex <- set_row_filters(ex, strand == {{ex_str}})
  ex <- set_coordinate_modifier(ex, Coordinates$Exact)
  
  psi <- tryCatch({percent_spliced_in(list(lq), list(rq), list(ex), min_count = 10)}, error= function(e){})
  
  return(psi)
  
}

library(GenomicRanges)

# Function to deal with undefined junctions (more than 1 possibilitie downstream/upstream)
calculate_psi<-function(jxlist){
  exon_str <- as.character(unique(strand(jxlist$upstream)))
  
  upstream_jx <- jxlist[[ifelse(exon_str =="+", 2, 1)]]
  downstream_jx <- jxlist[[ifelse(exon_str =="+", 1 ,2)]]

  upstream_jx<-gsub(":\\+|:\\-", "", as.character(upstream_jx))
  downstream_jx<-gsub(":\\+|:\\-", "", as.character(downstream_jx))
  
  if(length(upstream_jx)>1 ){
    
    psi<-lapply(upstream_jx, function(x){
      sk_jx<-paste0(seqnames(GRanges(downstream_jx)), ":", 
                    start(GRanges(x)),"-",
                          end(GRanges(downstream_jx)))
        
      psi_query_calc(x, downstream_jx, sk_jx, exon_str)
    })
    
    psi<-rbindlist(psi, fill=T, idcol = "TXNAME") # Collapse results of possible isoforms
    
  }
  else if(length(downstream_jx)>1 ){
    
     psi<-lapply(downstream_jx, function(x){
      sk_jx<-paste0(seqnames(GRanges(upstream_jx)), ":", 
                    start(GRanges(upstream_jx)),"-",
                          end(GRanges(x)))
      
     psi_query_calc(upstream_jx,x,sk_jx, exon_str)
    }) 
     psi<-rbindlist(psi, fill=T, idcol = "TXNAME")
  }
  else{ # If there is only one option in both sides, this line is run
    sk_jx<-paste0(seqnames(GRanges(upstream_jx)), ":", 
                  start(GRanges(upstream_jx)),"-",
                  end(GRanges(downstream_jx)))
    
    psi<-psi_query_calc(upstream_jx,downstream_jx,sk_jx, exon_str)
    
  }
  
  return(psi)
    
}

# Apply the function to all the exons
psi_values<-lapply(paired_introns, calculate_psi)     


# Adding the names as GeneSymbol_coordinate
names(psi_values)<-paste0(new_genes$gene_symbol, "_", as.character(gene_coord))

# Generating a single table
psi_values<-rbindlist(psi_values,use.names = T, idcol = "SYMBOL_COORDINATE", fill = T)

# Filtering and transformation to a matrix format
psi_values_w<-psi_values %>% 
  filter(psi>0) %>% 
  mutate(SYMBOL_COORDINATE_TX=paste0(ifelse(is.na(TXNAME), "",paste0( TXNAME, "-")), SYMBOL_COORDINATE)) %>% 
  dplyr::select(SYMBOL_COORDINATE_TX, sample_id, psi) %>% 
  pivot_wider(names_from = sample_id, values_from=psi, values_fill=0)

# Saving the results 
fwrite(psi_values_w, "output/20201209_PSI_table_exons.csv")


