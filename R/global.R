
##Functions######################################################
# ========================== Functions ======================== #
##pairwise
selectorUI <- function(id){
  ns = NS(id)
  rev_label <- function(label){
    paste0(rev(unlist(strsplit(label,split = '-'))),collapse = '-')
  }
  tags$div(
    fluidRow(
      div(style="display:inline-block;vertical-align:middle;margin-left:15px",
          HTML(paste0(rev_label(ns('Contrast group'))))
      ),
      br(),
      div(style="display:inline-block;vertical-align:middle;margin-left:15px;margin-right:10px",
          uiOutput(ns('Treatment'))),
      div(style="display:inline-block;vertical-align:middle;height:45px",
          HTML('VS')
      ),
      div(style="display:inline-block;vertical-align:middle;margin-left:10px;margin-right:15px",
          uiOutput(ns('Control'))
      )
    ),
    id = paste0('param', id)
  )
}

selectorServer <- function(input, output, session, thisList, multiple=F){
  ns = session$ns
  output$Treatment <- renderUI({
    selectInput(
      inputId = ns('Treatment'),multiple = multiple,
      label = NULL,
      width = '250px',
      choices = thisList, selected = NULL)
  })
  
  output$Control <- renderUI({
    selectInput(
      inputId = ns('Control'),multiple = multiple,
      label = NULL,
      width = '250px',
      choices = thisList, selected = NULL)
  })
}

##mean
selectorUI_mean <- function(id){
  ns = NS(id)
  rev_label <- function(label){
    paste0(rev(unlist(strsplit(label,split = '-'))),collapse = '-')
  }
  tags$div(
    fluidRow(
      div(style="display:inline-block;vertical-align:middle;margin-left:15px",
          HTML(paste0(rev_label(ns('Contrast group'))))
      ),
      br(),
      div(style="display:inline-block;vertical-align:middle;margin-left:15px;margin-right:10px",
          uiOutput(ns('Treatment_mean'))),
      div(style="display:inline-block;vertical-align:middle;height:45px",
          HTML('VS')
      ),
      div(style="display:inline-block;vertical-align:middle;margin-left:10px;margin-right:15px",
          uiOutput(ns('Control_mean'))
      )
    ),
    id = paste0('param_mean', id)
  )
}

selectorServer_mean <- function(input, output, session, thisList, multiple=T){
  ns = session$ns
  output$Treatment_mean <- renderUI({
    selectInput(
      inputId = ns('Treatment_mean'),multiple = multiple,
      label = NULL,
      width = '250px',
      choices = thisList, selected = NULL)
  })
  
  output$Control_mean <- renderUI({
    selectInput(
      inputId = ns('Control_mean'),multiple = multiple,
      label = NULL,
      width = '250px',
      choices = thisList, selected = NULL)
  })
}

##pairwise difference
selectorUI_pgdiff <- function(id){
  ns = NS(id)
  rev_label <- function(label){
    paste0(rev(unlist(strsplit(label,split = '-'))),collapse = '-')
  }
  tags$div(
    fluidRow(
      div(style="display:inline-block;vertical-align:middle;margin-left:15px",
          HTML(paste0(rev_label(ns('Contrast group'))))
      ),
      br(),
      div(style="display:inline-block;vertical-align:middle;margin-left:15px;margin-right:10px",
          uiOutput(ns('Treatment_pgdiff'))),
      div(style="display:inline-block;vertical-align:middle;height:45px",
          HTML('VS')
      ),
      div(style="display:inline-block;vertical-align:middle;margin-left:10px;margin-right:15px",
          uiOutput(ns('Control_pgdiff'))
      )
    ),
    id = paste0('param_pgdiff', id)
  )
}

selectorServer_pgdiff <- function(input, output, session, thisList, multiple=T){
  ns = session$ns
  output$Treatment_pgdiff <- renderUI({
    selectInput(
      inputId = ns('Treatment_pgdiff'),multiple = multiple,
      label = NULL,
      width = '250px',
      choices = thisList, selected = NULL)
  })
  
  output$Control_pgdiff <- renderUI({
    selectInput(
      inputId = ns('Control_pgdiff'),multiple = multiple,
      label = NULL,
      width = '250px',
      choices = thisList, selected = NULL)
  })
}


###=====>TS trend
selectorUI_tstrend <- function(id){
  ns = NS(id)
  rev_label <- function(label){
    paste0(rev(unlist(strsplit(label,split = '-'))),collapse = '-')
  }
  tags$div(
    fluidRow(
      div(style="display:inline-block;vertical-align:middle;margin-left:15px",
          HTML(paste0(rev_label(ns('Time series groups'))))
      ),
      br(),
      div(style="display:inline-block;vertical-align:middle;margin-left:15px;margin-right:10px",
          uiOutput(ns('tstrend_groups')))
    ),
    id = paste0('param_tstrend', id)
  )
}

selectorServer_tstrend  <- function(input, output, session, thisList,selected=NULL){
  ns = session$ns
  output$tstrend_groups <- renderUI({
    selectInput(
      inputId = ns('tstrend_groups'),multiple = T,
      label = NULL,
      width = '250px',
      choices = thisList, selected = selected)
  })
}



create.folders <- function(wd=getwd()){
  folder.name <- paste0(wd,'/',c('data','figure','result','report'))
  for(i in folder.name){
    if(!file.exists(i))
      dir.create(i,recursive = T)
  }
}

showmessage <- function(text='Done! Plots in pdf and png formats are saved to figure folder in App working directory',
                        duration = 10, type='default',showNoteify=T,showTime=T,...){
  
  if(showTime){
    if(showNoteify)
      text <- paste0(text,' (click "x" to close) (',Sys.time(),')') else text <- paste0(text,' (',Sys.time(),')')
  }
    
  
  cat(text,'\n')
  # if(grepl('\n',text))
  #   cat(text) else cat(text,'\n')
  if(showNoteify)
    showNotification(text,duration = duration,type=type,...)
}

startmessage <- function(text,id = 'message_id'){
  showmessage(paste(text,'...'),
              action = HTML("<i style='font-size:35px;' class='fas fa-coffee'> ... ...</i>"),
              # action = HTML("<i style='font-size:35px;' class='fas fa-dizzy'> ... ...</i>"),
              duration = NULL,id = id )
}

endmessage <- function(text,id = 'message_id'){
  showmessage(paste(text,'-> Done!!!'),
              action = HTML("<i style='font-size:35px;' class='fas fa-coffee'> ... ...</i>"),
              duration = NULL,showNoteify = F)
  cat('\n')
  removeNotification(id = id)
}

createLog <- function(logfile, path=getwd()) {
  if (file.exists(paste0(path, logfile))) {
    file.remove(paste0(path, logfile))
  }
  fid <- file(file.path(path, logfile), open = "wt")
  sink(fid, type = "message", split = F)
  sink(fid, append = T, type = "output", split = T)
  # warning("Use closeAllConnections() in the end of the script")
}
