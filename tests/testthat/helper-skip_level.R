skip_level <- function(test_lvl){
  lvl <- if (nzchar(s <- Sys.getenv("LAV_TEST_LEVEL")) &&
             is.finite(s <- as.numeric(s))) s
  else 1
  
  if (test_lvl > lvl) testthat::skip(paste("test level", test_lvl, ">", lvl))
}

