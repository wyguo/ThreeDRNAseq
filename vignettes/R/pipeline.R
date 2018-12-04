#' Dfferential expression (DE) and alternative splicing (DAS) analysis
#' @description Two pipelines are provided to perform differential expression (DE), differential
#' alternative splicing (DAS) gene and differential transcript usage (DTU) transcript analysis.
#' @param dge a numeric matrix of read counts or a \code{\link[edgeR]{DGEList}} object.
#' Counts must be non-negative and NAs are not permitted.
#' @param method method to use in the \code{\link{edgeR.pipeline}}. Options are "glm" or "glmQL".
#' @param contrast a vector of contrast groups for expression comparisons,
#' e.g. c('B-A','C-A') compares conditions "B" and "C" to condition "A".
#' @param span width of the \code{\link{lowess}} smoothing window
#' as a proportion, which is passed to the \code{\link[limma]{voom}} function to estimate
#' expression mean-variance trend.
#' @param design design matrix with rows of samples and columns of conditions to fit a linear model.
#' @param deltaPS a matrix of \eqn{\Delta}PS values, with rows of transcripts and columns of contrast groups.
#' @param diffAS logical, whether to perform DAS analysis. If \code{TRUE}, the \code{\link[limma]{diffSplice}} (limma pipeline)
#' or \code{\link[edgeR]{diffSpliceDGE}} (edgeR pipeline) function is used to peform DAS gene and DTU transcript analysis.
#' @param adjust.method method the adjust p-values for multiple comparisons. Default is "BH" and full details can be found in \code{\link{p.adjust}.}
#' @details
#' \bold{Abbreviation}
#' \itemize{
#' \item{DE: differential expressed genes/transcripts.}
#' \item{DAS: differential alternative spliced genes.}
#' \item{DTU: differential transcript usage transcripts.}
#' }
#'
#' \bold{Statistics}
#' \itemize{
#'     \item{
#'     Contrast groups: refer to the conditions to compare, e.g. "B-A" and "C-A" compares the conditions
#'     "B" and "C" to condition "A".
#'     }
#'     \item{
#'     Adjusted p-value: in the expression comparitive analysis, each gene/transcript is fitted with a statistical model, reporting a p-value.
#'     However the testing is performed over a substantial number of genes/transcripts and each testing has a probability of making errors.
#'     Hence, all p-values in the analsysis are adjusted by controlling the false discovery rate (FDR). Please see the R function \code{\link{p.adjust}}
#'     for details.
#'     }
#'     \item{
#'     \eqn{L_2FC}: \eqn{\log_2} fold change of expression based on contrast groups.
#'    }
#'    \item{
#'    \eqn{\Delta}PS: Transcript level PS (percent of splice) is defined as the ratio of transcript average abundance of conditions divided by the
#'    average gene abundance. \eqn{\Delta}PS is the PS differences of conditions based on the contrast groups. The abundance to calculate PS values
#'    can be TPMs or read counts.
#'    }
#' }
#' \bold{DE gene/transcript analysis}
#'
#' A gene/transcript is identified as DE in a contrast group if \eqn{L_2FC} of expression \eqn{\geq} a cut-off and with adjusted p-value < a cut-off.
#' 
#' Two pipelines are provided for DE analysis: 
#' 
#' ``limma-voom'' and ``edgeR''. The ``edgeR'' pipeline includes two methods: ``glmQL'' (Genewise Negative Binomial Generalized Linear Models with 
#' Quasi-likelihood Tests) and ``glm'' (Genewise Negative Binomial Generalized Linear Models). In limma pipeline, L2FCs are calculated based count 
#' per million reads (CPMs) while in edgeR, L2FCs are based on read counts. From several real RNA-seq data analyses, high consensus is achieved between 
#' ``limma-voom'' and ``glmQL'' (>90%) and they have more stringent controls of false discovery rate than the ``glm'' method.
#'
#' \bold{DAS gene and DTU transcript analysis}
#'
#' To test differential expression, the gene level expression is compared to transcript level expression in the contrast groups.
#' A gene is DAS in a contrast group if adjusted p-value < cut-off and at least one transcript of the gene with \eqn{\Delta}PS \eqn{\geq} a cut-off.
#' A transcript is DTU if adjusted p-value < a cut-off and \eqn{\Delta}PS \eqn{\geq} a cut-off.
#'
#' Two pipelines for AS analysis:
#' \itemize{
#'     \item{
#'     limma pipeline: To identify DAS genes, the expression of each transcript is compared to the weighted average expression of all the transcripts
#'     for the same gene (weight on transcript expression variance; the average expression can be treated as gene level expression). P-value of each
#'     test is converted to genewise p-value by using F-test or Simes method. To identify DTU transcript, each transcript is compared to the weighted
#'     average of all the other transcripts for the same gene. See the \code{\link{diffSplice}} function in \code{\link{limma}} package for details.
#'     }
#'     \item{
#'     edgeR pipeline: To identify DTU transcripts, log-fold-change of each transcript is compared to log-fold-change of entire gene
#'     (sum of all transcripts). The DTU transcript test p-values are summarised to genewise p-value by using F-test or Simes method
#'     to report DAS genes. See \code{\link{diffSpliceDGE}} function in \code{\link{edgeR}} R package for details.
#'     }
#' }
#' @return a list object with elements:
#' \itemize{
#'    \item{
#'    \code{voom.object:} the results of \code{\link[limma]{voom}} function to estimate expression mean-variance trend.
#'    }
#'    \item{
#'    \code{fit.lmFit:} the results of \code{\link[limma]{lmFit}} to fit a linear regression model on conditions.
#'    }
#'    \item{
#'    \code{fit.glm:} the results of \code{\link[edgeR]{glmFit}} in the edgeR pipeline.
#'    }
#'    \item{
#'    \code{fit.glmQL:} the  results of \code{\link[edgeR]{glmQLFit}} in the edgeR pipeline.
#'    }
#'    \item{
#'    \code{fit.contrast:} the results of \code{\link[limma]{contrasts.fit}}, i.e. incorporate the contrast groups to the linear regression model.
#'    }
#'    \item{
#'    \code{fit.eBayes:} the results of \code{\link[limma]{eBayes}}.
#'    }
#'    \item{
#'    \code{DE.pval.list:} a list object of which each element is the result of applying \code{\link[limma]{topTable}} or \code{\link[edgeR]{topTags}}
#'    to individual contrast groups.
#'    }
#'    \item{
#'    \code{DE.pval:} a data.frame object of adjusted p-values, with rows of genes/transcripts and columns of contrast groups.
#'    }
#'    \item{
#'    \code{DE.lfc:} a data.frame object of \eqn{L_2FC}, with rows of genes/transcripts and columns of contrast groups.
#'    }
#'    \item{
#'    \code{DE.stat:} a data.frame object with first column of target (genes/transcripts), second column of contrast (contrast groups), third column of
#'    adj.pval and fourth column of log2FC.
#'    }
#'    \item{
#'    \code{DE.stat.overalltest:} a data.frame object of testing statistics over all contrast groups rather than testing in individual contrast groups.
#'    }
#' }
#' If alternative splicing analysis is performed (\code{diffAS=TRUE}), more results are returned:
#' \itemize{
#'    \item{
#'    \code{fit.splice:} the results of \code{\link[limma]{diffSplice}} or \code{\link[edgeR]{diffSpliceDGE}} for individual contrast groups.
#'    }
#'    \item{
#'    \code{DTU.pval.list:} a list object of which each element is the DTU transcript result of \code{\link[limma]{topSplice}} or
#'    \code{\link[edgeR]{topSpliceDGE}} for individual contrast groups.
#'    }
#'    \item{
#'    \code{DTU.pval:} a data.frame object of adjusted p-values, with rows of transcripts and columns of contrast groups.
#'    }
#'    \item{
#'    \code{DTU.deltaPS:} a data.frame object of \eqn{\Delta}PS, with rows of transcripts and columns of contrast groups.
#'    }
#'    \item{
#'    \code{DTU.stat:} a data.frame object with first column of target (transcripts), second column of contrast groups, third column of adj.pval and
#'     fourth column of deltaPS.
#'    }
#'    \item{
#'    \code{maxdeltaPS:} a data.frame object with rows of genes and columns of maximum \eqn{\Delta}PS of the transcripts in the same gene in different
#'    contrast groups.
#'    }
#'    \item{
#'    \code{DAS.pval.F.list:} a list object of which each element is the DAS gene result of "F-test" using \code{\link[limma]{topSplice}} or
#'    \code{\link[edgeR]{topSpliceDGE}} in individual contrast groups.
#'    }
#'    \item{
#'    \code{DAS.pval.F:} a data.frame object of adjusted p-value with rows of genes and columns of contrast groups.
#'    }
#'    \item{
#'    \code{DAS.F.stat:} the DAS gene statistics using "F-test".
#'    }
#'    \item{
#'    \code{DAS.pval.simes.list:} a list object of "Simes-test" using \code{\link[limma]{topSplice}} or
#'    \code{\link[edgeR]{topSpliceDGE}} in individual contrast groups.
#'    }
#'    \item{
#'    \code{DAS.pval.simes:} a data.frame object of DAS gene adjusted p-values using "Simes" method.
#'    }
#'    \item{
#'    \code{DAS.simes.stat:} the DAS gene statistics using "Simes-test".
#'    }
#' }
#' @export
#' @rdname limma.pipeline
#' @seealso R packages: \code{\link{limma}} and \code{\link{edgeR}}.

limma.pipeline <- function(dge,
                           contrast,
                           span=0.5,
                           design,
                           deltaPS=NULL,
                           diffAS=F,
                           adjust.method='BH'){
  start.time <- Sys.time()
  results <- list()
  ##########################################################
  ##--limma voom
  message('Limma-voon to estimate mean-vriance trend ...')
  voom.object<-voom(dge,design,plot=F,span = span)
  results$voom.object <- voom.object
  targets <- rownames(voom.object$E)

  ##########################################################
  ##--Fit a basic linear model
  message('Fit a basic linear model ...')
  fit.lmFit <- lmFit(voom.object, design)
  results$fit.lmFit <- fit.lmFit

  ##########################################################
  ##--Fit a basic linear model
  message('Fit the contrast model ...')
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)
  print(paste0('Contrast groups: ',paste0(contrast,collapse = '; ')))
  fit.contrast<-contrasts.fit(fit.lmFit, contrast.matrix)
  results$fit.contrast<-fit.contrast

  ##########################################################
  ##--Fit a eBayes model
  message('Fit a eBayes model ...')
  fit.eBayes<-eBayes(fit.contrast)
  results$fit.eBayes<-fit.eBayes

  ##########################################################
  ##--Testing statistics for each contrast group
  message('Testing for each contrast group ...')
  DE.pval.list <- lapply(contrast,function(i){
    x <- topTable(fit.eBayes,
                  coef=i,
                  adjust.method =adjust.method,
                  number = Inf)
    x <- x[targets,]
    x
  })
  names(DE.pval.list) <- contrast
  results$DE.pval.list<-DE.pval.list

  ###---DE pval and lfc
  DE.pval <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$adj.P.Val))
  DE.lfc <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$logFC))
  rownames(DE.pval) <- rownames(DE.lfc) <- targets
  colnames(DE.pval) <- colnames(DE.lfc) <- contrast

  DE.stat <- summaryStat (x = DE.pval,y = DE.lfc,
                          target = rownames(DE.pval),
                          contrast = contrast,
                          stat.type = c('adj.pval','log2FC'))

  # DE.stat <- cbind(DE.pval,DE.lfc)
  # colnames(DE.stat) <- c(paste0('pval:',contrast),paste0('lfc:',contrast))
  results$DE.pval<-DE.pval
  results$DE.lfc<-DE.lfc
  results$DE.stat<-DE.stat

  ##########################################################
  ##--Testing statistics for across all contrast groups
  message('Testing across all contrast groups ...')
  DE.stat.overalltest<-topTable(fit.eBayes,number = Inf,
                            coef = contrast,adjust.method =adjust.method )
  # DE.stat.overalltest <- DE.stat.overalltest[targets,]
  col.idx <- gsub('-','.',contrast)
  col.idx <- grep(paste0(col.idx,collapse = '|'),colnames(DE.stat.overalltest))
  DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                    contrast='overall',DE.stat.overalltest,row.names = NULL)
  colnames(DE.stat.overalltest)[col.idx] <- gsub('[.]','-',colnames(DE.stat.overalltest)[col.idx])
  results$DE.stat.overalltest<-DE.stat.overalltest

  if(diffAS){
    message('Fit a splicing model ...')
    if(is.null(deltaPS))
      stop('Please provide deltaPS for DAS analysis...')
    fit.splice<-diffSplice(fit.contrast, geneid = 'GENEID')
    results$fit.splice<-fit.splice

    # ##########################################################
    ##---DTU transcripts
    ##DTU transcript pval list
    genes.idx <- unique(fit.splice$genes$GENEID)
    trans.idx <- unique(fit.splice$genes$TXNAME)

    DTU.pval.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="t", number=Inf, FDR=10000)
      rownames(y)<-y$TXNAME
      z <- y[trans.idx,]
      z
    })
    names(DTU.pval.list) <- contrast

    ##DTU transcript pvals
    DTU.pval <- lapply(contrast,function(i){
      x <- DTU.pval.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    DTU.pval <- do.call(cbind,DTU.pval)
    colnames(DTU.pval) <- contrast

    ##DTU transcript deltaPS
    DTU.deltaPS <- deltaPS[rownames(DTU.pval),,drop=F]

    DTU.stat <- summaryStat (x = DTU.pval,y = DTU.deltaPS,
                            target = rownames(DTU.pval),
                            contrast = contrast,
                            stat.type = c('adj.pval','deltaPS'))

    # DTU.stat <- cbind(DTU.pval,DTU.deltaPS)
    # colnames(DTU.stat) <- c(paste0('pval:',contrast),paste0('deltaPS:',contrast))

    results$DTU.pval.list<-DTU.pval.list
    results$DTU.pval<-DTU.pval
    results$DTU.deltaPS<-DTU.deltaPS
    results$DTU.stat<-DTU.stat

    # ##########################################################
    ##---DAS genes
    ##---max deltaPS
    maxdeltaPS <- by(deltaPS[fit.splice$genes$TXNAME,,drop=F],
                     INDICES = fit.splice$genes$GENEID,
                     function(x){
                       apply(x,2,function(i) i[abs(i)==max(abs(i))][1])
                     },simplify = F)
    maxdeltaPS <- do.call(rbind,maxdeltaPS)
    maxdeltaPS <- maxdeltaPS[genes.idx,,drop=F]
    results$maxdeltaPS<-maxdeltaPS

    ##---F test
    DAS.pval.F.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="F", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      y[genes.idx,]
    })
    names(DAS.pval.F.list) <- contrast

    DAS.pval.F <- lapply(contrast,function(i){
      x <- DAS.pval.F.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })

    DAS.pval.F <- do.call(cbind,DAS.pval.F)
    DAS.pval.F <- DAS.pval.F[genes.idx,,drop=F]
    colnames(DAS.pval.F) <- contrast

    DAS.F.stat <- summaryStat (x = DAS.pval.F,y = maxdeltaPS,
                             target = rownames(DAS.pval.F),
                             contrast = contrast,
                             stat.type = c('adj.pval','maxdeltaPS'))

    # DAS.F.stat <- cbind(DAS.pval.F,maxdeltaPS)
    # colnames(DAS.F.stat) <- c(paste0('pval:',contrast),paste0('MaxdeltaPS:',contrast))
    #
    results$DAS.pval.F.list<-DAS.pval.F.list
    results$DAS.pval.F<-DAS.pval.F
    results$DAS.F.stat<-DAS.F.stat

    ##---simes test
    DAS.pval.simes.list<-lapply(contrast,function(i){
      y<-topSplice(fit.splice, coef=i, test="simes", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      y[genes.idx,]
    })
    names(DAS.pval.simes.list) <- contrast

    DAS.pval.simes <- lapply(contrast,function(i){
      x <- DAS.pval.simes.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })

    DAS.pval.simes <- do.call(cbind,DAS.pval.simes)
    DAS.pval.simes <- DAS.pval.simes[genes.idx,,drop=F]
    colnames(DAS.pval.simes) <- contrast

    DAS.simes.stat <- summaryStat (x = DAS.pval.simes,y = maxdeltaPS,
                               target = rownames(DAS.pval.simes),
                               contrast = contrast,
                               stat.type = c('adj.pval','maxdeltaPS'))
    #
    # DAS.simes.stat <- cbind(DAS.pval.simes,maxdeltaPS)
    # colnames(DAS.simes.stat) <- c(paste0('pval:',contrast),paste0('MaxdeltaPS:',contrast))

    results$DAS.pval.simes.list<-DAS.pval.simes.list
    results$DAS.pval.simes<-DAS.pval.simes
    results$DAS.simes.stat<-DAS.simes.stat
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste0('Time for analysis: ',round(time.taken,3)))
  message('Done!!! ')
  return(results)
}



#############################################################################################
##-----------------------------------------------------------------------------------------
#' @rdname limma.pipeline
#' @export
edgeR.pipeline <- function(dge,
                           method=c('glm','glmQL'),
                           design,
                           contrast,
                           diffAS=F,
                           deltaPS=NULL,
                           adjust.method='BH'){

  start.time <- Sys.time()
  results <- list()
  method <- match.arg(method,c('glm','glmQL'))
  targets <- rownames(dge$counts)
  message('Estimate dispersion ...')
  Disp <- estimateGLMCommonDisp(dge, design)
  Disp <- estimateGLMTagwiseDisp(Disp, design)
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)

  switch(method,
         glm = {
           message('Fit Genewise Negative Binomial Generalized Linear Models...')
           fit <- glmFit(Disp, design = design)
           results$fit.glm <- fit

           ###individual test
           message('Fit the contrast model ...')
           DE.pval.list <- lapply(contrast,function(i){
             x <- glmLRT(fit, contrast=contrast.matrix[,i,drop=F])
             y <- topTags(x,n = Inf,adjust.method = adjust.method)
             y <- y[targets,]
           })
           ###overall test
           x <- glmLRT(fit, contrast=contrast.matrix)
           DE.stat.overalltest <- topTags(x,n = Inf,adjust.method = adjust.method)
           # DE.stat.overalltest <- DE.stat.overalltest[targets,]
           DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                             contrast='overall',DE.stat.overalltest,row.names = NULL)
           results$DE.stat.overalltest<-DE.stat.overalltest
         },
         glmQL = {
           message('Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests...')
           fit <- glmQLFit(Disp,design = design)
           results$fit.glmQL <- fit

           ###individual test
           message('Fit the contrast model ...')
           DE.pval.list <- lapply(contrast,function(i){
             x <- glmQLFTest(fit, contrast=contrast.matrix[,i,drop=F])
             y <- topTags(x,n = Inf,adjust.method = adjust.method)
             y <- y[targets,]
           })

           ###overall test
           x <- glmQLFTest(fit, contrast=contrast.matrix)
           DE.stat.overalltest <- topTags(x,n = Inf,adjust.method = adjust.method)
           # DE.stat.overalltest <- DE.stat.overalltest[targets,]
           DE.stat.overalltest <- data.frame(target=rownames(DE.stat.overalltest),
                                             contrast='overall',DE.stat.overalltest,row.names = NULL)
           results$DE.stat.overalltest<-DE.stat.overalltest
         }
  )
  names(DE.pval.list) <- contrast
  results$DE.pval.list <- DE.pval.list

  DE.pval <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$table$FDR))
  DE.lfc <- do.call(cbind,lapply(DE.pval.list,FUN = function(x) x$table$logFC))
  rownames(DE.pval) <- rownames(DE.lfc) <- targets
  colnames(DE.pval) <- colnames(DE.lfc) <- contrast

  DE.stat <- summaryStat (x = DE.pval,y = DE.lfc,
                          target = rownames(DE.pval),
                          contrast = contrast,
                          stat.type = c('adj.pval','log2FC'))
  results$DE.pval<-DE.pval
  results$DE.lfc<-DE.lfc
  results$DE.stat<-DE.stat
  ##################################################################
  ###---DAS
  if(diffAS){
    fit.splice <- lapply(contrast,function(i){
      message(paste0('\ndiffSpliceDGE of contrast: ', i))
      diffSpliceDGE(fit, contrast=contrast.matrix[,i,drop=F], geneid="GENEID")
    })
    names(fit.splice) <- contrast
    results$fit.splice<-fit.splice

    genes.idx <- unique(fit.splice[[1]]$genes$GENEID)
    trans.idx <- unique(fit.splice[[1]]$genes$TXNAME)

    DTU.pval.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="exon", number=Inf, FDR=10000)
      rownames(y)<-y$TXNAME
      z <- y[trans.idx,]
      z
    })
    names(DTU.pval.list) <- contrast

    ##DTU transcript pvals
    DTU.pval <- lapply(contrast,function(i){
      x <- DTU.pval.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })
    DTU.pval <- do.call(cbind,DTU.pval)
    colnames(DTU.pval) <- contrast

    DTU.deltaPS <- deltaPS[rownames(DTU.pval),,drop=F]
    DTU.stat <- summaryStat (x = DTU.pval,y = DTU.deltaPS,
                             target = rownames(DTU.pval),
                             contrast = contrast,
                             stat.type = c('adj.pval','deltaPS'))

    results$DTU.pval.list<-DTU.pval.list
    results$DTU.pval<-DTU.pval
    results$DTU.deltaPS<-DTU.deltaPS
    results$DTU.stat<-DTU.stat

    # ##########################################################
    ##---DAS genes
    ##---max deltaPS
    maxdeltaPS <- by(deltaPS[fit.splice[[1]]$genes$TXNAME,,drop=F],
                     INDICES = fit.splice[[1]]$genes$GENEID,
                     function(x){
                       apply(x,2,function(i) i[abs(i)==max(abs(i))][1])
                     },simplify = F)
    maxdeltaPS <- do.call(rbind,maxdeltaPS)
    maxdeltaPS <- maxdeltaPS[genes.idx,,drop=F]
    results$maxdeltaPS<-maxdeltaPS

    ##---F test
    DAS.pval.F.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="gene", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      z <- y[genes.idx,]
      z
    })
    names(DAS.pval.F.list) <- contrast

    DAS.pval.F <- lapply(contrast,function(i){
      x <- DAS.pval.F.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })

    DAS.pval.F <- do.call(cbind,DAS.pval.F)
    DAS.pval.F <- DAS.pval.F[genes.idx,,drop=F]
    colnames(DAS.pval.F) <- contrast

    DAS.F.stat <- summaryStat (x = DAS.pval.F,y = maxdeltaPS,
                               target = rownames(DAS.pval.F),
                               contrast = contrast,
                               stat.type = c('adj.pval','maxdeltaPS'))

    results$DAS.pval.F.list<-DAS.pval.F.list
    results$DAS.pval.F<-DAS.pval.F
    results$DAS.F.stat<-DAS.F.stat

    ##---Simes test
    DAS.pval.simes.list<-lapply(contrast,function(i){
      fit.splice.i <- fit.splice[[i]]
      y<-topSpliceDGE(fit.splice.i, test="Simes", number=Inf, FDR=10000)
      rownames(y)<-y$GENEID
      z <- y[genes.idx,]
      z
    })
    names(DAS.pval.simes.list) <- contrast

    DAS.pval.simes <- lapply(contrast,function(i){
      x <- DAS.pval.simes.list[[i]]
      y <- data.frame(pval=x$FDR)
      rownames(y) <- rownames(x)
      y
    })

    DAS.pval.simes <- do.call(cbind,DAS.pval.simes)
    DAS.pval.simes <- DAS.pval.simes[genes.idx,,drop=F]
    colnames(DAS.pval.simes) <- contrast

    DAS.simes.stat <- summaryStat (x = DAS.pval.simes,y = maxdeltaPS,
                                   target = rownames(DAS.pval.simes),
                                   contrast = contrast,
                                   stat.type = c('adj.pval','maxdeltaPS'))

    results$DAS.pval.simes.list<-DAS.pval.simes.list
    results$DAS.pval.simes<-DAS.pval.simes
    results$DAS.simes.stat<-DAS.simes.stat
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste0('Time for analysis: ',round(time.taken,3)))
  message('Done!!! ')
  return(results)
}


