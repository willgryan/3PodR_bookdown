### Load packages
suppressPackageStartupMessages({
  #Packages have to be loaded individually for renv static analysis
  library(magrittr)
  library(tidyverse)
  library(HGNChelper)
  library(purrr)
  library(ggpubr)
  library(ggrepel)
  library(pheatmap)
  library(knitr)
  library(BioPathNet)
  library(drugfindR)
  library(jsonlite)
  library(httr)
  library(enrichR)
  library(factoextra)
  library(babelgene)
  library(ggVennDiagram)
  library(PAVER)
  library(ggupset)
  library(DT)
  library(circlize)
  library(seriation)
  library(randomcoloR)
  library(ComplexHeatmap)
  
  options(repos = BiocManager::repositories())
  
  source("R/utils.R")
  
  set.seed(123)
  
})

### Functions for reading in transcriptomic data

read_deg <- function(file) {
  file %>%
    read_csv() %>%
    select(1:3) %>%
    rename(Symbol = 1, log2FoldChange = 2, pvalue = 3) %>%
    filter(if_all(everything(), ~ . != '')) %>%
    drop_na() %>%
    mutate(Symbol = str_trim(Symbol)) %>%
    distinct(Symbol, .keep_all = TRUE)
}

read_counts <- function(file) {
  file %>%
    read_csv() %>%
    rename(Symbol = 1) %>%
    filter(if_all(everything(), ~ . != '')) %>%
    drop_na() %>%
    mutate(Symbol = str_trim(Symbol)) %>%
    distinct(Symbol, .keep_all = TRUE)
}

#Input: table where 1st column is gene symbols to correct
fix_hgnc <- function(X, map) {
  HGNChelper::checkGeneSymbols(X %>% pull(1), map = map) %>%
    select(Symbol = 1, Suggested.Symbol) %>%
    mutate(Suggested.Symbol = if_else(Suggested.Symbol == "" | is.na(Suggested.Symbol), Symbol, Suggested.Symbol),
           Suggested.Symbol = if_else(str_detect(Suggested.Symbol, ".+?(?= ///)"), str_extract(Suggested.Symbol, ".+?(?= ///)"), Suggested.Symbol)) %>%
    inner_join(X %>% rename(Symbol = 1), by = "Symbol") %>%
    select(-Symbol, Symbol = Suggested.Symbol)
}

###Setup global report state###
global_state <- yaml::yaml.load_file("extdata/_variables.yml")

#This is the HGNCHelper precomputed mapping file
global_state$humanmap <-
  read_csv(paste0("extdata/", global_state$humanmap)) %>%
  filter(Symbol != "NA" & Approved.Symbol != "NA") #TODO investigate later why HGNChelper is mapping an "NA" gene symbol?

#These are gene annotations for each species
global_state$hgnc <- read_csv(paste0("extdata/", global_state$hgnc)) %>%
  select(Symbol = symbol, Name = name)

global_state$mgi <- read_csv(paste0("extdata/", global_state$mgi)) %>%
  select(Symbol = `Marker Symbol`, Name = `Marker Name`)

global_state$rgd <- read_csv(paste0("extdata/", global_state$rgd)) %>%
  select(Symbol = SYMBOL, Name = NAME)

global_state$map <- case_when(global_state$species == "human" ~ list(global_state$hgnc),
                              global_state$species == "mouse" ~ list(global_state$mgi),
                              global_state$species == "rat" ~ list(global_state$rgd)) %>% pluck(1)
#This is the scraped LINCS FDA data
global_state$lincs_fda <- read_csv(paste0("extdata/", global_state$lincs_fda))

#Theses are files for PAVER
global_state$embeddings <- readRDS(paste0("extdata/", global_state$embeddings))

global_state$ontology_index <- readRDS(paste0("extdata/", global_state$ontology_index))

#Read counts and design if specified
if (!is_empty(global_state$design) & !is_empty(global_state$counts)) {
  global_state$using_counts <- TRUE
  
  global_state$design <-
    read_csv(paste0("extdata/", global_state$design))

  global_state$counts <-
    read_counts(paste0("extdata/", global_state$counts))
  
  if (global_state$species == "human") {
    #Correct HGNC symbols of the counts if human data
    global_state$counts %<>% fix_hgnc(global_state$humanmap)
  }
} else {
  global_state$using_counts <- FALSE
}

#Load input DEG data
global_state$data %<>%
  map(~ update_list(., name = paste0(.$group1, " vs ", .$group2), #Generate comparison name
                    data = read_deg(paste0("extdata/", .$file)),
                    results = list())) %>% #Read DEG data
  set_names(map_chr(., "name")) %>% #Set names on data list
  map_if(~ global_state$species == "human", modify_at, "data", fix_hgnc, global_state$humanmap) %>% #Correct HGNC symbols #TODO: How to fix Mouse or rat symbols?
  map(~ update_list(., bpn = BioPathNet::prepare_data(.$data$Symbol, .$data$log2FoldChange, .$data$pvalue))) #Create BPN object for each dataset

#Create global results store
global_state$results <- list()