quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#Input: Vector of Gene symbols
#Output: Tibble of enrichR results
do_enrichr <- function(X) {
  dbs = c("GO_Biological_Process_2023", 
           "GO_Molecular_Function_2023", 
           "GO_Cellular_Component_2023")
  
  columns = c("Biological_Process",
               "Molecular_Function",
               "Cellular_Component")
  
  quiet(X %>%
    enrichr(databases = dbs) %>% 
    map2(columns, ~ mutate(.x, namespace = .y)) %>%
    bind_rows %>%
    filter(Adjusted.P.value <= .05) %>%
    extract(Term, "GOID", "(GO:\\d+)", remove = FALSE) %>%
    extract(Term, "Term", "(.*?)\\(GO:\\d+\\)"))
}

#Input: Vector of LINCS signature ids from drugFindr
#Output: Tibble of iLINCS signature metadata
get_ilincs_metadata <- function(X) {
  url <- "https://www.ilincs.org/api/SignatureMeta/findMany"
  body <- list(signatures = toJSON(X))
  
  metadata <- POST(url, body = body, encode = "json") %>%
    content(as = "text") %>%
    fromJSON() %>%
    pluck("data") %>%
    select(TargetSignature = signatureid, tissue, integratedMoas, GeneTargets) #where(~ any(!is.na(.x)))
}

make_table <- function(X, caption=NULL) {
  DT::datatable(
    X,
    rownames = FALSE,
    caption = htmltools::tags$caption(
      style = 'text-align: left;',
      caption
    ),
    options = list(
      scrollX = TRUE,
      scrollY = TRUE,
      paging = TRUE,
      fixedHeader = TRUE,
      pageLength = 10
    )
  )
}