test_get_enrichment <- function(){

  string_db <- STRINGdb$new( version="11", species=511145, score_threshold=800 ) 

  proteins <- string_db$get_proteins()
  hits <- proteins$protein_external_id[1:50]
  clusters <- string_db$get_clusters(hits)
  checkTrue(length(clusters) > 1)

}
