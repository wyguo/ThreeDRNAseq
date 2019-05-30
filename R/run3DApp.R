#' Run 3D RNA-seq App
#' @param wd directory to run the App. Default is the project working directory.
#' @param ... additional arguments passed to \code{runApp} in \code{\link{shiny}} R package.
#' @export
#' @examples 
#' run3DApp()
run3DApp <- function(wd=getwd(),...){
  app.from <- system.file("app", package = "ThreeDRNAseq")
  app.to <- paste0(wd,'/app')
  if(!all(list.files(app.from) %in% list.files(app.to)))
  message(paste0('Downloading 3D RNA-seq App:\n ',wd))
    file.copy(from=app.from, to=wd, 
              overwrite = T, recursive = T, 
              copy.mode = T)
    message('Done!!!')
  shiny::runApp(app.to)
}

#' Download example data for 3D RNA-seq analysis
#' @param wd directory to save example data. Default is the project working directory.
#' @export
#' @examples 
#' getExample()
getExample <- function(wd=getwd()){
  exp.from <- system.file("example_data", package = "ThreeDRNAseq")
  exp.to <- paste0(wd,'/example_data')
  if(!all(list.files(exp.from) %in% list.files(exp.to)))
    message(paste0('Downloading example data to:\n ',wd))
  file.copy(from=exp.from, to=wd, 
            overwrite = T, recursive = T, 
            copy.mode = T)
  message('Done!!!')
}