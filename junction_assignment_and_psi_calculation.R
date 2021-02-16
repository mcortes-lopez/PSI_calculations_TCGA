library(data.table)
new_genes <- fread("data/genes_new.csv", data.table = F)


library(tidyr)
library(dplyr)


# Get searchable coordinates
gene_coord <- new_genes %>%
  dplyr::select(starts_with("chr"), strand) %>%
  dplyr::select(-chr_size) %>%
  dplyr::rename("seqnames" = "chr_name", "start" = "chr_start", "end" = "chr_end") %>%
  GenomicRanges::GRanges()


# The coordinates are in the hg19 version on the genome. Normally we use the version hg38
# We will proceed to 1) convert the coordinates to hg38


library(rtracklayer)
path <- "data/hg19ToHg38.over.chain"
ch <- import.chain(path)
ch

new_exons <- liftOver(gene_coord, ch) %>%
  unlist()


# Now we need to get the upstream and downstream junctions.
# For this we need to overlap the exon with the introns



library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg38_exons <- exons(txdb, columns = c("TXID", "TXNAME", "EXONID"))
hg38_introns <- intronsByTranscript(txdb, use.names = T)

ovp <- findOverlapPairs(new_exons, hg38_exons, type = "any")

# The start of the exons is off by 1 nt

proper_ovp <- ovp[(width(first(ovp)) - width(second(ovp))) == 1]

library(IRanges)


# Getting upstream and downstream exons
paired_introns <- lapply(seq_along(proper_ovp), function(exon_idx) {
  introns_un <- unlist(hg38_introns)

  upstream_int <- precede(second(proper_ovp[exon_idx]), introns_un, select = "all")
  downstream_int <- follow(second(proper_ovp[exon_idx]), introns_un, select = "all")

  tx_names_up <- names(introns_un[subjectHits(upstream_int)])
  tx_names_down <- names(introns_un[subjectHits(downstream_int)])

  match_tx <- intersect(tx_names_up, tx_names_down)

  tx_names_exon <- second(proper_ovp)[exon_idx]$TXNAME


  tx_matched <- intersect(unlist(tx_names_exon), match_tx)


  matched_up_introns <- introns_un[subjectHits(upstream_int)[tx_names_up %in% tx_matched]]

  matched_down_introns <- introns_un[subjectHits(downstream_int)[tx_names_down %in% tx_matched]]



  return(list(
    "upstream" = unique(matched_up_introns),
    "downstream" = unique(matched_down_introns)
  ))
})



# Generating the junctions table
# If "+" -> upstream_ss = start(upstream)  and  downstream_ss = end(downstream)
# If "-" -> upstream_ss = end(downstream)  and  downstream_ss = start(upstream)

library(snapcount)


# Calculation of PSI based on the SNAPTRON manual and using the recently updated data from TCGA (V2)  or GTEx (v2)
psi_query_calc <- function(ljx, rjx, skipjx, ex_str, dtb) {
  print(list(ljx, rjx, skipjx, ex_str))

  lq <- QueryBuilder(compilation = dtb, regions = ljx)

  lq <- set_row_filters(lq, strand == {{ ex_str }})
  lq <- set_coordinate_modifier(lq, Coordinates$Exact)
  # right inclusion query
  rq <- QueryBuilder(compilation = dtb, regions = rjx)
  rq <- set_row_filters(rq, strand == {{ ex_str }})
  rq <- set_coordinate_modifier(rq, Coordinates$Exact)
  # exclusion query
  ex <- QueryBuilder(compilation = dtb, regions = skipjx)
  ex <- set_row_filters(ex, strand == {{ ex_str }})
  ex <- set_coordinate_modifier(ex, Coordinates$Exact)

  psi <- tryCatch(
    {
      percent_spliced_in(list(lq), list(rq), list(ex), min_count = 10)
    },
    error = function(e) {}
  )

  return(psi)
}

library(GenomicRanges)

# Function to deal with undefined junctions (more than 1 possibility downstream/upstream)
calculate_psi <- function(jxlist, dtbs) {
  exon_str <- as.character(unique(strand(jxlist$upstream)))

  upstream_jx <- jxlist[[ifelse(exon_str == "+", 2, 1)]]
  downstream_jx <- jxlist[[ifelse(exon_str == "+", 1, 2)]]

  upstream_jx <- gsub(":\\+|:\\-", "", as.character(upstream_jx))
  downstream_jx <- gsub(":\\+|:\\-", "", as.character(downstream_jx))

  if (length(upstream_jx) > 1) {
    psi <- lapply(upstream_jx, function(x) {
      sk_jx <- paste0(
        seqnames(GRanges(downstream_jx)), ":",
        start(GRanges(x)), "-",
        end(GRanges(downstream_jx))
      )

      psi_query_calc(x, downstream_jx, sk_jx, exon_str, dtbs)
    })

    psi <- rbindlist(psi, fill = T, idcol = "TXNAME") # Collapse results of possible isoforms
  }
  else if (length(downstream_jx) > 1) {
    psi <- lapply(downstream_jx, function(x) {
      sk_jx <- paste0(
        seqnames(GRanges(upstream_jx)), ":",
        start(GRanges(upstream_jx)), "-",
        end(GRanges(x))
      )

      psi_query_calc(upstream_jx, x, sk_jx, exon_str, dtbs)
    })
    psi <- rbindlist(psi, fill = T, idcol = "TXNAME")
  }
  else { # If there is only one option in both sides, this line is run
    sk_jx <- paste0(
      seqnames(GRanges(upstream_jx)), ":",
      start(GRanges(upstream_jx)), "-",
      end(GRanges(downstream_jx))
    )

    psi <- psi_query_calc(upstream_jx, downstream_jx, sk_jx, exon_str, dtbs)
  }

  return(psi)
}




# TCGA DATA ------------------------------------------------------------------
# Apply the function to all the exons
psi_values <- lapply(paired_introns, function(intronpair) {
  calculate_psi(intronpair, "tcgav2")
})


# Adding the names as GeneSymbol_coordinate
names(psi_values) <- paste0(new_genes$gene_symbol, "_", as.character(gene_coord))

# Generating a single table
psi_values <- rbindlist(psi_values, use.names = T, idcol = "SYMBOL_COORDINATE", fill = T)

# Filtering

psi_values_w <- psi_values %>%
  filter(psi > 0) %>%
  mutate(SYMBOL_COORDINATE_TX = paste0(ifelse(is.na(TXNAME), "", paste0(TXNAME, "-")), SYMBOL_COORDINATE))

# Intron coordinates ---------------------------------------------------------


names(paired_introns) <- paste0(new_genes$gene_symbol, "_", as.character(gene_coord))

introns_table<-lapply(paired_introns, function(x){
  event_df <- lapply(x, as.data.frame)
  event_df <- rbindlist(event_df, use.names = T, idcol = "tx")
  event_df$tx_name <- unlist(sapply(x, names))
  return(event_df)
  }) %>% 
  rbindlist(., idcol = "exon_id")


tx_one<-psi_values_w %>% 
  filter(is.na(psi_values_w$TXNAME)) %>% 
  pull(SYMBOL_COORDINATE) %>% 
  unique()


hg38_convert<- data.frame("exon_hg38_coordinate" = as.character(second(proper_ovp)), 
           "exon_id" = names(paired_introns))


introns_table %>% 
  filter(tx_name %in% unique(psi_values_w$TXNAME)  | exon_id %in% tx_one) %>%
  left_join(., hg38_convert, by = "exon_id") %>% 
  fwrite("output/intron_information.tab", sep = "\t", quote = F, row.names = F)



####


# Adding TCGA IDs in the shape of: TCGA-XX-XX-XX

metadata <- fread("http://snaptron.cs.jhu.edu/data/tcgav2/samples.tsv", header = T)

ids <- metadata[, c("rail_id", "cgc_sample_id")]
ids <- ids %>%
  mutate(sample_id_red = gsub("[A-Z]$", "", cgc_sample_id)) %>%
  filter(sample_id_red != "")

# Merge and transform to matrix format
psi_values_w <- psi_values_w %>%
  inner_join(., ids, by = c("sample_id" = "rail_id")) %>%
  group_by(sample_id_red, SYMBOL_COORDINATE_TX) %>%
  arrange(desc(psi)) %>% # Some ids have more than one sample, so we select the one with a biggest psi
  slice_head() %>%
  ungroup()


# Generate a new file for the individual read counts 

junctions_reads<-psi_values_w %>% 
  dplyr::select(SYMBOL_COORDINATE_TX,
                sample_id_red,
                inclusion_group1_coverage,
                inclusion_group2_coverage,
                exclusion_group_coverage) %>% 
  pivot_longer(cols = c(inclusion_group1_coverage,
                        inclusion_group2_coverage,
                        exclusion_group_coverage), 
               values_to = "junction_reads") %>% 
  pivot_wider(names_from = sample_id_red, 
              values_from = junction_reads, 
              values_fill = 0)


# PSI counts
psi_values_w <- psi_values_w %>%
  dplyr::select(SYMBOL_COORDINATE_TX, sample_id_red, psi) %>%
  pivot_wider(
    names_from = sample_id_red,
    values_from = psi,
    values_fill = 0
  )


# GTExs DATA ------------------------------------------------------------------
# Apply the function to all the exons
psi_values <- lapply(paired_introns, function(intronpair) {
  calculate_psi(intronpair, "gtexv2")
})


# Adding the names as GeneSymbol_coordinate
names(psi_values) <- paste0(new_genes$gene_symbol, "_", as.character(gene_coord))

# Generating a single table
psi_values <- rbindlist(psi_values, use.names = T, idcol = "SYMBOL_COORDINATE", fill = T)

# Filtering

psi_values_w_gtex <- psi_values %>%
  filter(psi > 0) %>%
  mutate(SYMBOL_COORDINATE_TX = paste0(ifelse(is.na(TXNAME), "", paste0(TXNAME, "-")), SYMBOL_COORDINATE))


# Adding GTEx IDs in the shape of: TCGA-XX-XX-XX

metadata <- fread("http://snaptron.cs.jhu.edu/data/gtexv2/samples.tsv", header = T)

ids <- metadata[, c("rail_id", "SAMPID")]


# Merge and transform to matrix format
psi_values_w_gtex <- psi_values_w_gtex %>%
  inner_join(., ids, by = c("sample_id" = "rail_id")) %>%
  group_by(SAMPID, SYMBOL_COORDINATE_TX) %>%
  arrange(desc(psi)) %>% # Some ids have more than one sample, so we select the one with a biggest psi
  slice_head() %>%
  ungroup()



junctions_reads_gtex<-psi_values_w_gtex %>% 
  dplyr::select(SYMBOL_COORDINATE_TX,
                SAMPID,
                inclusion_group1_coverage,
                inclusion_group2_coverage,
                exclusion_group_coverage) %>% 
  pivot_longer(cols = c(inclusion_group1_coverage,
                        inclusion_group2_coverage,
                        exclusion_group_coverage), 
               values_to = "junction_reads") %>% 
  pivot_wider(names_from = SAMPID, 
              values_from = junction_reads, 
              values_fill = 0)



psi_values_w_gtex <- psi_values_w_gtex %>%
  dplyr::select(SYMBOL_COORDINATE_TX, SAMPID, psi) %>%
  pivot_wider(
    names_from = SAMPID,
    values_from = psi,
    values_fill = 0
  )


# Merge TCGA and GTEx --------------------------------------------
psi_mtx <- left_join(psi_values_w, psi_values_w_gtex, by = "SYMBOL_COORDINATE_TX")

# Saving the results
#fwrite(psi_mtx, "output/20210104_PSI_table_exons.csv.gz")


# Junction reads 
jxreads_mtx <- left_join(junctions_reads, junctions_reads_gtex, by = c("SYMBOL_COORDINATE_TX", "name"))

jxlist<-split.data.frame(jxreads_mtx[, -1], 
                 jxreads_mtx$SYMBOL_COORDINATE_TX)


lapply(names(jxlist), function(jx_file_name){
  fwrite(jxlist[[jx_file_name]], paste0("output/jx_reads/", jx_file_name, ".csv"))
  })




