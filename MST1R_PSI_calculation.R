# Get MST1R junctions
library(snapcount)
library(dplyr)


# Function from the previous script
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


# RON data
# exon_id	tx	seqnames	start	end	width	strand	tx_name	exon_hg38_coordinate
# MST1R_chr3:49895961-49896107:-	upstream	chr3	49895881	49895960	-	ENST00000296474.8	chr3:49895961-49896107:-
# MST1R_chr3:49895961-49896107:-	downstream	chr3	49896108	49896194	-	ENST00000296474.8	chr3:49895961-49896107:-

ron_tcga <- psi_query_calc("chr3:49895881-49895960", "chr3:49896108-49896194", "chr3:49895881-49896194", "-", "tcgav2")
ron_gtex <- psi_query_calc("chr3:49895881-49895960", "chr3:49896108-49896194", "chr3:49895881-49896194", "-", "gtexv2")



library(data.table)
library(tidyr)
# TCGA -------------------------------------------------------------------------------------------------------------------
metadata <- fread("http://snaptron.cs.jhu.edu/data/tcgav2/samples.tsv", header = T)

ids <- metadata[, c("rail_id", "cgc_sample_id")]
ids <- ids %>%
  mutate(sample_id_red = gsub("[A-Z]$", "", cgc_sample_id)) %>%
  filter(sample_id_red != "")

# Merge and transform to matrix format
psi_tab <- ron_tcga %>%
  filter(psi > 0)


psi_tab <- psi_tab %>%
  inner_join(., ids, by = c("sample_id" = "rail_id")) %>%
  group_by(sample_id_red) %>%
  arrange(desc(psi)) %>% # Some ids have more than one sample, so we select the one with a biggest psi
  slice_head() %>%
  ungroup()


# Generate a new file for the individual read counts

junctions_reads <- psi_tab %>%
  dplyr::select(
    sample_id_red,
    inclusion_group1_coverage,
    inclusion_group2_coverage,
    exclusion_group_coverage
  ) %>%
  pivot_longer(
    cols = c(
      inclusion_group1_coverage,
      inclusion_group2_coverage,
      exclusion_group_coverage
    ),
    values_to = "junction_reads"
  ) %>%
  pivot_wider(
    names_from = sample_id_red,
    values_from = junction_reads,
    values_fill = 0
  )


# PSI counts
psi_tab <- psi_tab %>%
  dplyr::select(sample_id_red, psi) %>%
  pivot_wider(
    names_from = sample_id_red,
    values_from = psi,
    values_fill = 0
  )


# GTEx ------------------------------------------------------------------------------------------------------------
# Adding GTEx IDs in the shape of: TCGA-XX-XX-XX

metadata <- fread("http://snaptron.cs.jhu.edu/data/gtexv2/samples.tsv", header = T)

ids <- metadata[, c("rail_id", "SAMPID")]


# Merge and transform to matrix format
psi_gtex <- ron_gtex %>%
  filter(psi > 0) %>%
  inner_join(., ids, by = c("sample_id" = "rail_id")) %>%
  group_by(SAMPID) %>%
  arrange(desc(psi)) %>% # Some ids have more than one sample, so we select the one with a biggest psi
  slice_head() %>%
  ungroup()



junctions_reads_gtex <- psi_gtex %>%
  dplyr::select(
    SAMPID,
    inclusion_group1_coverage,
    inclusion_group2_coverage,
    exclusion_group_coverage
  ) %>%
  pivot_longer(
    cols = c(
      inclusion_group1_coverage,
      inclusion_group2_coverage,
      exclusion_group_coverage
    ),
    values_to = "junction_reads"
  ) %>%
  pivot_wider(
    names_from = SAMPID,
    values_from = junction_reads,
    values_fill = 0
  )



psi_values_w_gtex <- psi_gtex %>%
  dplyr::select(SAMPID, psi) %>%
  pivot_wider(
    names_from = SAMPID,
    values_from = psi,
    values_fill = 0
  )



# Merge matrices -----------------------------------------------------------------------------------------
ron_id <- data.frame("SYMBOL_COORDINATE_TX" = "MST1R_chr3:49895961-49896107:-")
psi_mtx <- cbind.data.frame(ron_id, psi_tab)
psi_mtx <- cbind.data.frame(psi_mtx, psi_values_w_gtex)

jx_mtx <- merge(junctions_reads, junctions_reads_gtex, by = "name")

# fwrite(psi_mtx, "output/20210219_PSI_table_MST1R.csv.gz")
# fwrite(jx_mtx, "output/jx_reads/MST1R_chr3:49895961-49896107:-.csv")
