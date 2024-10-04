#' convert markdown doc to readme file for github
#'
#' rmd2md()
#'
#' For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#'
#' @param rmd_fl name of the markdown file (without the .Rmd)
#' @param md_fl  name of the output; default is 'README'
#'
#' @examples
#' # convert the mymarkdown.Rmd to README.md
#' rmd2md(rmd_fl="mymarkdown",md_fl = "README")
#'
#' @export
rmd2md <- function(rmd_fl ,md_fl = "README"){
  library(rmarkdown)
  render(paste0(rmd_fl,".Rmd"), md_document(variant = "markdown_github"),params=F)
  file.copy(from=paste0(rmd_fl,".md"),to=paste0(md_fl,".md"),overwrite=T)
}

