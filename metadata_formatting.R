# Metadata processing

library(dplyr)
library(data.table)

# TCGA ----------------------------------------------------------------------------

tcga_meta <- fread("http://snaptron.cs.jhu.edu/data/tcgav2/samples.tsv")

tcga_meta <- tcga_meta %>%
  mutate(sample = gsub("[A-Z]$", "", gdc_cases.samples.submitter_id)) %>%
  select(
    sample,
    gdc_cases.project.name,
    gdc_cases.samples.sample_type,
    gdc_cases.project.primary_site
  ) %>%
  rename(
    description = gdc_cases.project.name,
    tissue = gdc_cases.project.primary_site,
    type = gdc_cases.samples.sample_type
  ) %>%
  mutate(
    set = "TCGA",
    simple_type = "Tumor"
  )



# GTEx ----------------------------------------------------------------------------
gtex_meta <- fread("http://snaptron.cs.jhu.edu/data/gtexv2/samples.tsv")

gtex_meta <- gtex_meta %>%
  select(SAMPID, SMTS, SMTSD) %>%
  rename(
    sample = SAMPID,
    tissue = SMTS,
    description = SMTSD
  ) %>%
  mutate(
    set = "GTEX",
    simple_type = "Normal",
    type = "Normal"
  )


gtex_meta[!grepl("GTEX", gtex_meta$sample)]$simple_type <- "Tumor"
gtex_meta[!grepl("GTEX", gtex_meta$sample)]$type <- "Tumor"

gtex_meta <- gtex_meta %>%
  filter(sample != "") # Drop empty samples


# Merge metadata and save ------------------------------------------------------------

meta_snaptron <- rbind(tcga_meta, gtex_meta)
meta_snaptron %>%
  select(
    sample,
    type,
    tissue,
    description,
    set,
    simple_type
  )


#fwrite(meta_snaptron, "output/20210205_metadata_tcga_gtex.tab")
