test_get_homologs_besthits <- function() {

  string_db <- STRINGdb$new( version="11", species=9606, score_threshold=400 ) 
  stringid = string_db$mp('CDK1') 
  homologs = string_db$get_homologs_besthits( stringid )
  checkTrue(nrow(homologs) > 0)

}
