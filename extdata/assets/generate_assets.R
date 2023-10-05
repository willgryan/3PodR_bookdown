library(tidyverse)
library(jsonlite)
library(httr)
library(HGNChelper)
library(XML)
library(fgsea)

date = Sys.Date()

###HGNCHelper map
HGNChelper_currentHumanMap <- HGNChelper::getCurrentHumanMap()

HGNChelper_currentHumanMap %>%
  write_csv(paste0("HGNChelper_currentHumanMap_", date, ".csv"))
  
###HGNC Gene Symbol Annotations
hgnc_complete_set <- read_tsv("ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt")

hgnc_complete_set %>%
  write_csv(paste0("hgnc_complete_set_", date, ".csv"))

###Mouse Gene Symbol Annotations
MRK_List2 <- read_tsv("http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt")
MRK_List2 %>%
  write_csv(paste0("MRK_List2_", date, ".csv"))

###Rat Gene Symbol Annotations
GENES_RAT <- read_tsv("https://download.rgd.mcw.edu/pub/data_release/GENES_RAT.txt", comment = "#")
GENES_RAT %>%
  write_csv(paste0("GENES_RAT_", date, ".csv"))

###GMT files

##Human
url <- "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/GO/"
gmt <- XML::getHTMLLinks(url) %>%
  .[str_detect(., "GOALL_no_GO")]

fgsea::gmtPathways(paste0(url, gmt)) %>% 
  fgsea::writeGmtPathways(gmt)

##Rat
url <- "http://download.baderlab.org/EM_Genesets/current_release/Rat/symbol/GO/"
gmt <- XML::getHTMLLinks(url) %>%
  .[str_detect(., "GOALL_no_GO")]

fgsea::gmtPathways(paste0(url, gmt)) %>% 
  fgsea::writeGmtPathways(gmt)

##Mouse
url <- "http://download.baderlab.org/EM_Genesets/current_release/Mouse/symbol/GO/"
gmt <- XML::getHTMLLinks(url) %>%
  .[str_detect(., "GOALL_no_GO")]

fgsea::gmtPathways(paste0(url, gmt)) %>% 
  fgsea::writeGmtPathways(gmt)

###PAVER files
embeddings = readRDS(url("https://github.com/willgryan/PAVER_embeddings/raw/main/2023-03-06/embeddings_2023-03-06.RDS"))

term2name = readRDS(url("https://github.com/willgryan/PAVER_embeddings/raw/main/2023-03-06/term2name_2023-03-06.RDS"))

###SCRAPE FOR FDA METADATA
scrape=F

get_data <- function(X) {
  url = "https://lincsportal.ccs.miami.edu/sigc-api/small-molecule/fetch"
  
  r = RETRY("GET", url, query = list(
    limit = 100,
    returnSignatures = FALSE,
    page=X
  ), times = 1000)
  
  metadata = r %>%
    content(as = "text") %>%
    fromJSON()%>%
    pluck("data")
  
  metadata
}
if(scrape) {
  data <- tibble()
  for(page in 1:(round(21230/100)+1)) {
    tmp = get_data(page) %>% mutate(page = page)
    data <- data %>% bind_rows(tmp)
    cat(c("PAGE: ", page))
    system("sleep 1")
  }
  
  processed <- data %>%
    select(c("canonical_inchi_key", "canonical_inchi", "perturbagen_id", "sm_name", "canonical_smiles", "max_fda_phase")) %>%
    write_csv("LINCS_FDA_Phases.csv")
}

