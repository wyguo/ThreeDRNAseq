plotSave <- function(p,filename,folder2save=getwd(),
                     width,height,units = 'in',res = 150,
                     type=c('all','png','pdf'),...){
  
  type <- match.arg(type,c('all','png','pdf'))
  if(type == 'all' | type == 'png'){
    png(filename = paste0(folder2save,'/',filename,'.png'),width = width,height = height,units = units,res = res,...)
    grid.draw(p())
    dev.off()
  }
  if(type == 'all' | type == 'pdf'){
    pdf(file = paste0(folder2save,'/',filename,'.pdf'),width = width,height = height)
    grid.draw(p())
    dev.off()
  }
}
