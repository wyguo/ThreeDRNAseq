#' Generate html, pdf and word report from Rmd (RMardkown) file
#' @param input.Rmd the Rmd file used to generate the report.
#' @param report.folder the folder to save the reports.
#' @param type the ouput document type. Options are "all", "html_document", "pdf_document" and "word_document".
#' @return the generated html, pdf or word report will be saved into the provided \code{report.folder}.
#' @export
#' 
generate.report <- function(input.Rmd,report.folder,
                            type=c('all','html_document','pdf_document','word_document')){
  if(!file.exists(report.folder))
    dir.create(report.folder,recursive = T)
  type <- match.arg(type,c('all','html_document','pdf_document','word_document'))
  if(type=='all'){
    try(rmarkdown::render(input = input.Rmd,
                          params = list(doc.type = "html_document"
                          ),
                          output_format = "html_document",
                          output_dir = report.folder))
    try(rmarkdown::render(input = input.Rmd, 
                          params = list(doc.type = "pdf_document"
                          ),
                          output_format = "pdf_document",output_dir = report.folder))
    
    try(rmarkdown::render(input = input.Rmd, 
                          params = list(doc.type = "word_document"
                          ),
                          output_format = "word_document",output_dir = report.folder))
  } else {
    try(rmarkdown::render(input = input.Rmd, 
                          params = list(doc.type = type
                          ),
                          output_format = type,output_dir = report.folder))
  }
}
