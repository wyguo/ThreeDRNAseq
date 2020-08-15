#' Isoform switch analysis for time-series data
#'
#' This function is used to search and score transcript isoform switch in time-series expression.
#'
#' @details The detailed steps:
#'
#' \if{html}{\figure{figures001.png}{options: width="80\%" alt="Figure: figures_001.png"}}
#' \if{latex}{\figure{figures001.pdf}{options: width=7cm}}
#'
#' \bold{Figure 1:} Isoform switch analysis methods. Expression data with 3 replicates for
#' each condition/time point is simulated for isoforms \eqn{iso_i} and \eqn{iso_j}.
#' The points in the plots represent the samples and the black lines connect the average
#' of samples. (A) is the iso-kTSP algorithm for comparisons of two conditions \eqn{c_1} and \eqn{c_2}.
#' Time-Series Isoform Switch (TSIS) tool is designed for detection and characterization of
#' isoform switches for time series data shown in figure (B). The time-series with 6 time
#' points is divided into 4 intervals by the intersection points of average expression.
#'
#' \bold{Step 1: determine switch points.}
#'
#' Given that a pair of isoforms \eqn{iso_i} and \eqn{iso_j} may have a number of switches
#' in a time-series, we have offered two approaches to search for the switch time points in TSIS:
#'
#' \itemize{
#'     \item{The first approach takes the average values of the replicates for each time point
#'     for each transcript isoform. Then it searches for the cross points of the average
#'     value of two isoforms across the time points (seen in Figure 1(B)).
#'     }
#'     \item{The second approach uses natural spline curves to fit the time-series data
#'     for each transcript isoform and find cross points
#'       of the fitted curves for each pair of isoforms.
#'    }
#' }
#'
#' In most cases, these two methods produce very similar results. However, average values
#' of expression may lose precision without having information of backward and forward time points.
#' The spline method fit time-series of expression with control points (depending on spline
#' degree of freedom provided) and weights of several neighbours to obtain designed precision
#' (Hastie and Tibshirani, 1990). The spline method is useful to find global trends in the time-series data
#' when the data is very noisy. But it may lack details of isoform switch in the local region.
#' It is recommended that users use both average and spline method to search for the switch
#' points and examine manually when inconsistent results were produced by the above two methods.
#'
#'
#' \bold{Step 2: define switch scoring measures}
#'
#' We define each transcript isoform switch by 1) the switch point \eqn{P_i},
#' 2) time points between switch points \eqn{P_(i-1)}  and \eqn{P_i}  as interval  \eqn{I_1}
#' before switch \eqn{P_i}  and 3) time points between switch points \eqn{P_i}  and \eqn{P_(i+1)}
#' as interval \eqn{I_2}  after the switch \eqn{P_i}  (see Figure 1(B)). We defined five measurements
#' to characterize each isoform switch. The first two are the probability/frequency of switch and
#' the sum of average sample differences before and after switch, which are similar to Score 1 and
#' Score 2 in iso-kTSP method (Sebestyen, et al., 2015) (see Figure 1(A))).
#'
#' \itemize{
#'   \item{
#'     Firstly, for a switch point \eqn{P_i}  of two isoforms \eqn{iso_i}  and \eqn{iso_j}  with interval \eqn{I_1}  before the switch and  interval \eqn{I_2}  after the switch ( Figure1 (B)), Score 1 is defined as
#'     \deqn{S_1 (iso_i,iso_j |I_1,I_2)=|p(iso_i>iso_j |I_1)+p(iso_i<iso_j |I_2)-1|,}
#'     Where \eqn{p(iso_i>iso_j |I_1)} and \eqn{p(iso_i<iso_j |I_2)} are the frequencies/probabilities that the samples of one isoform is greater or less than in the other in corresponding intervals.
#'   }
#'   \item{
#'     Secondly, instead of rank differences as in iso-kTSP to avoid possible ties, we directly use the average abundance differences. The sum of mean differences of samples in intervals \eqn{I_1} and \eqn{I_2} are calculated as
#'     \deqn{S_2 (iso_i,iso_j |I_1,I_2)=d(iso_i,iso_j |I_1)+d(iso_i,iso_j |I_2)}
#'     Where \eqn{d(iso_i,iso_j |I_k)} is the average expression difference in interval \eqn{I_k,k=1,2} defined as
#'     \deqn{d(iso_i, iso_j |I_k)=\frac{1}{|I_k|}\sum_{m_{I_k}} |exp(iso_i |s_{m_{I_k}},I_k)-exp(iso_j |s_{m_{I_k}},I_k) |}
#'     \eqn{|I_k |} is the number of samples in interval \eqn{I_k} and \eqn{exp(iso_i |s_(m_(I_k ) ),I_k)} is the expression of \eqn{iso_i} of sample \eqn{s_(m_(I_k ) )} in interval \eqn{I_k}. }
#'   \item{
#'     Thirdly, paired t-tests were carried out to test if there are significant differences for the two switched isoforms within the intervals before and after the switch.  }
#'   \item{
#'     Fourthly, the numbers of time points in intervals \eqn{I_1} and \eqn{I_2} were also provided, which indicate whether this switch is transient or a long lived change.
#'   }
#'   \item{
#'     Finally, Isoforms with high negative correlations across the time points may point to the splicing regulation have antagonistic effects on the switched isoforms. Thus they are of great interest for investigation of alternative splicing regulations. As an additional score, we calculated the Pearson correlation of two isoforms across the whole time series.
#'   }
#' }
#'
#' For TSIS analysis details and examples, please go the the user manual on Github: \url{https://github.com/wyguo/TSIS}.
#' 
#' Note: \code{\link{iso.switch.shiny}} is the same to \code{\link{iso.switch}} but used in the shiny App function \code{\link{TSIS.app}} to show progress
#' bar. \code{\link{iso.switch.pairwise}} is used to identify switches between pairwise conditions. The parameter \code{min.t.points} is set to 1.
#'
#' @param data.exp  Time-series isoform expression data with first row indicating the replicate labels and second row indicating the time points. The remained lines are isoform names in the first column followed by the expression values. All the replicates for each time point should be grouped together and the time points follow the sequential order.
#' @param mapping gene and isoform mapping table with gene names in first column  and transcript isoform names in the second column
#' @param times a numeric vector of time labellings of all the relicated samples, e.g. 1,1,1,2,2,2,3,3,3,... If a vector of characters is
#' provided, this vector will be converted into numeric.
#' @param rank logical (TRUE or FALSE). Should isoform expression be converted to rank of isoform expression in sample basis?
#' @param min.t.points if the number of time points in all intervals < \code{min.t.points},
#'   this pair of isoforms are not switch candidates since they only have transient switches.
#' @param min.difference If the mean of differences of average isoform expression or spline fitted expression < \code{min.difference},
#'   this pair of isoforms are supposed to tied together and they are not considered as switch candidates.
#' @param spline logical, whether to use spline method to fit isoform expression (TRUE) or mean expression of time points (FALSE).
#' @param spline.df the degree of freedom used in spline method. See \code{\link{ns}} in \code{\link{splines}} for details.
#' @param verbose logical, to track the running progressing (TRUE) or not (FALSE).
#'
#' @references
#'   1. Sebestyen E, Zawisza M, Eyras E: Detection of recurrent alternative splicing switches in tumor samples reveals novel signatures of cancer.
#'   Nucleic Acids Res 2015, 43(3):1345-1356.
#'
#'   2. Hastie, T.J. and Tibshirani, R.J. Generalized additive models. Chapter 7 of Statistical Models in S eds. Wadsworth & Brooks/Cole 1990.
#'
#' @return  a data frame of scores. The column names:
#'     \itemize{
#'       \item{ iso1,iso2: }{ the isoform pairs.}
#'       \item{ iso1.mean.ratio, iso2.mean.ratio: }{ the mean ratios of isoforms to their gene.}
#'       \item{ left.interval, left.interval: }{ The intervals before and after switch points.}
#'       \item{ x.value, y.value: }{ The values of x axis (time) and y axis (expression) coordinates of the switch points.}
#'       \item{ left.prob, right.prob: }{ the frequencies/probabilities that the samples of an isoform is greater or less than the other in left and right intervals, respectively.}
#'       \item{ left.dist, right.diff: }{ the average sample differences in intervals before and after switch, respectively.}
#'       \item{ left.pval, right.pval: }{ the paired t-test p-values of the samples in the intervals before and after switch points.}
#'       \item{ left.t.points, right.t.points: }{ the number of time points in intervals before and after the switch points.}
#'       \item{ prob: }{ Score1: the probability/frequency of switch. }
#'       \item{ diff: }{ Score2: the sum of average sample differences before and after switch.}
#'       \item{ cor: }{ Pearson correlation of two isoforms. }
#'     }
#' @export
#' @rdname iso.switch


iso.switch<-function(data.exp,
                     mapping,
                     times,
                     rank=F,
                     min.t.points=2,
                     min.difference=1,
                     spline=F,
                     spline.df=NULL,
                     verbose=T){
  
  if(nrow(data.exp)!=nrow(mapping))
    stop("Gene-isoform mapping table does not match to isoform expression table")
  
  times0 <- unique(times)
  if(!is.numeric(times))
    times <- as.numeric(factor(times,levels = unique(times)))
  
    
  
  start.time <- Sys.time()
  ##Step 1: Checking data information
  ##genes more than 2 transcripts
  genes0<-as.vector(mapping[,1])
  isoforms0<-as.vector(mapping[,2])
  
  ##expression ratio aboundance
  data2intersect.ratio<-apply(data.exp,1,mean)
  data2intersect.ratio<-rowratio(x = data2intersect.ratio,group = genes0)
  
  #choose genes with more than 2 transcripts
  idx<-table(genes0)[unique(genes0)]
  genes<-names(idx)[idx>1]
  isoforms<-isoforms0[which(genes0 %in% genes)]
  data.exp<-data.exp[isoforms,]
  if(rank)
    data2switch<-apply(data.exp,2,rank) else data2switch=data.exp
  
  
  
  message('Summarising input information ... ')
  message(' Input genes: ', length(unique(genes0)))
  message(' Genes with more than 2 isoforms: ', length(genes))
  message(' Average isoforms per gene for switch analysis: ', round(length(isoforms)/length(genes),3))
  
  
  
  ##step 2: Pickign isoform-pairs with intersection points...
  message(paste0('Step 1: Search for intersection points with ',if(spline) 'Spline method...' else 'Mean expression...'))
  ##Average values of time points
  ##data for searching isoform swtich points
  if(spline){
    message(' Spline fitting expression ...')
    if(is.null(spline.df))
      spline.df<-floor(2*(length(unique(times)-1))/3)
    data2intersect<-t(apply(data.exp,1,function(x) ts.spline(x,times  = times,df = spline.df)))
  } else {
    data2intersect<-t(rowmean(t(data.exp),group = times))
  }
  
  
  ##Find intersection points
  message(' Searching ...')
  iso.intersections<-list()
  t.start<-min(times)
  t.end<-max(times)
  
  if(verbose)
    pb <- txtProgressBar(min = 0, max = length(genes), style = 3,width = 75)
  for(i in 1:length(genes)){
    Sys.sleep(0)
    
    sub.gene<-genes[i]
    sub.isoform<-as.vector(isoforms0[which(genes0==sub.gene)])
    
    onegene<-data2intersect[sub.isoform,]
    comb<-combn(sub.isoform,2)
    # iso.intersections<-list()
    for(j in 1:ncol(comb)){
      iso1<-comb[1,j]
      iso2<-comb[2,j]
      ##
      x1=onegene[iso1,]
      x2=onegene[iso2,]
      
      if(max(abs(x1-x2)) < min.difference)
        next
      
      ##intersection points
      iso.inter <- ts.intersection(x1 = x1, x2 = x2,times = times)
      if(is.null(iso.inter))
        next
      iso.inter<-iso.inter[iso.inter$x.points<t.end & iso.inter$x.points>t.start,]
      ##check if have intersection points
      if(nrow(iso.inter)==0)
        next
      
      ##check if the switching lasting more than 2 consecutive time points
      check.consecutive<-c(t.start,iso.inter$x.points,t.end)
      check.consecutive<-data.frame(low=check.consecutive[-length(check.consecutive)],up=check.consecutive[-1])
      check.consecutive$low<-ceiling(check.consecutive$low)
      check.consecutive$up<-floor(check.consecutive$up)
      check.consecutive<-check.consecutive$up-check.consecutive$low+1
      
      if(max(check.consecutive) < min.t.points)
        next
      
      iso.inter<-data.frame(iso1=iso1,iso2=iso2,t(iso.inter))
      colnames(iso.inter)<-c('iso1','iso2',paste0('switch',1:(ncol(iso.inter)-2)))
      
      iso.intersections<-c(iso.intersections,setNames(object = list(iso.inter),nm = paste0(iso1,'_vs_',iso2)))
    }
    if(verbose)
      setTxtProgressBar(pb, i)
  }
  if(verbose)
    close(pb)
  
  if(length(iso.intersections)==0)
    return(NULL)
  
  message(' ',paste0(length(iso.intersections),' pairs of isoforms have intersection points.'))
  
  
  
  
  
  message('Step 2: Calculate scores for isoform switch ...')
  message(' Score 1: Switch frequencies/probabilities')
  message(' Score 2: Sum of average sample differences before and after switch.')
  message(' Score 3: P-values of sample differences before and after switch')
  message(' Score 4: Time points in each intervals')
  message(' Score 5: Pearson correlation of isoforms')
  ##step 2: calculate probablities
  
  iso.names<-names(iso.intersections)
  
  iso.scores<-data.frame()
  
  if(verbose)
    pb <- txtProgressBar(min = 0, max = length(iso.names), style = 3,width = 75)
  for(i in 1:length(iso.names)){
    Sys.sleep(0)
    #  i=1
    # i=which(grepl('AT1G01060.2_vs_AT1G01060.3',x = iso.names))
    iso<-unlist(strsplit(iso.names[i],'_vs_'))
    # message(iso)
    iso1<-iso[1]
    iso2<-iso[2]
    #generate dataset for iso for isoform switch analysis
    n.inters<-data.frame(iso.intersections[[iso.names[i]]][,-c(1:2)])
    colnames(n.inters)<-colnames(iso.intersections[[iso.names[i]]])[-c(1:2)]
    inters=as.numeric(n.inters[1,])
    
    interval=cut(times,unique(c(t.start,as.numeric(inters),t.end)),include.lowest = T,dig.lab = 5)
    # sub.data2swith=data.frame(interval=interval,t(data.exp.rank[c(iso1,iso2),]),iso.idx=c(-1,1)[as.numeric(interval)%%2+1])
    sub.data2swith=data.frame(interval=interval,t(data2switch[c(iso1,iso2),]))
    ##add a column of sign of differences
    
    diff.sign<-sign(sub.data2swith[,2]-sub.data2swith[,3])
    interval.idx<-as.numeric(interval) %%2
    interval.idx[interval.idx==0]<--1
    diff.sign<-diff.sign*interval.idx
    sub.data2swith$diff.sign<-diff.sign
    x0<-diff.sign[diff.sign!=0][1]
    #Define scores
    
    
    s<-by(data = sub.data2swith,INDICES =sub.data2swith[,1],FUN = function(x){
      ##score1:
      prob<-length(x[,4][x[,4]==x0])/length(x[,4])
      prob[is.na(prob)]<-0
      ##score2:
      x.diff<-x[,2]-x[,3]
      diff.mean<-mean(x.diff)
      # ##score3:
      
      pval <- tryCatch({
        t.test(x[,3],x[,2],paired = T)$p.value
      },error=function(e) return(0))
      pval[is.na(pval)]<-1
        
      #score4:
      # inter.length<-nrow(x)/nrep
      y<-as.character(x$interval[1])
      y<-as.numeric(unlist(strsplit(substr(x = y,start = 2,stop = nchar(y)-1),',',fixed = T)))
      inter.length <- length(which(unique(times)>=ceiling(min(y)) & unique(times)<= floor(max(y))))
      x<-data.frame(prob=prob,diff=diff.mean,pval=pval,inter.length=inter.length)
      x
    },simplify = T)
    s<-t(do.call(rbind,s))
    
    if(ncol(s)<3)
      lf.idx<-1:2 else lf.idx<-c(1,rep(2:(ncol(s)-1),each=2),ncol(s))
    
    ##score 1: prob
    score1.2side<-matrix(s[1,lf.idx],ncol=2,byrow = T)
    colnames(score1.2side)<-c('before.prob','after.prob')
    score1<-zoo::rollapply(as.numeric(s[1,]), width = 2, FUN = function(x) abs(sum(x)-1))
    
    ##score 2: difference
    score2<-abs(diff(s[2,]))
    score2.2side<-matrix(s[2,lf.idx],ncol=2,byrow = T)
    colnames(score2.2side)<-c('before.diff','after.diff')
    
    ##score 3: p-value
    score3<-matrix(s[3,lf.idx],ncol=2,byrow = T)
    colnames(score3)<-c('before.pval','after.pval')
    
    ##score3.adj: adjusted p-value
    score3.adj <- p.adjust(p = s[3,lf.idx],method = 'BH')
    score3.adj<-matrix(score3.adj,ncol=2,byrow = T)
    colnames(score3.adj)<-c('before.FDR','after.FDR')
    
    ##score 4: interval length
    score4<-matrix(s[4,lf.idx],ncol=2,byrow = T)
    colnames(score4)<-c('before.t.points','after.t.points')
    
    
    ###average ratio
    iso.mean.ratio<-data.frame(iso1.mean.ratio=data2intersect.ratio[iso1],iso2.mean.ratio=data2intersect.ratio[iso2])
    
    ###before and after intervals
    inter.lr<-matrix(unique(interval)[lf.idx],ncol=2,byrow = T)
    colnames(inter.lr)<-c('before.interval','after.interval')
    
    
    score<-data.frame(iso1=iso1,
                      iso2=iso2,
                      iso.mean.ratio,
                      inter.lr,
                      x.value=as.numeric(n.inters[1,]),
                      y.value=as.numeric(n.inters[2,]),
                      score1.2side,
                      score2.2side,
                      score3,
                      score3.adj,
                      score4,
                      prob=score1,
                      diff=score2,
                      cor=suppressWarnings(cor(as.numeric(data.exp[iso1,]),as.numeric(data.exp[iso2,]))),
                      row.names = NULL)
    
    iso.scores<-rbind(score,iso.scores)
    if(verbose)
      setTxtProgressBar(pb, i)
  }
  if(verbose)
    close(pb)
  ##sort
  iso.scores<-iso.scores[order(iso.scores$prob,decreasing = T),]
  iso.scores <- data.frame(iso.scores,row.names = NULL)
  if(!is.numeric(times0)){
    iso.scores$numeric.coordinates <- NA
    iso.scores$numeric.coordinates[1] <- paste0(paste0(times0,': ',1:length(times0)),collapse = '; ')
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste0('Time for analysis: ',round(time.taken,3),' ',attributes(time.taken)$units))
  
  message('Done!!! ')
  return(iso.scores)
}


#' @export
#' @rdname iso.switch

iso.switch.shiny<-function(data.exp,mapping,times,rank=F,
                           min.t.points=2,min.difference=1,spline=F,spline.df=NULL,verbose=T){
  
  if(nrow(data.exp)!=nrow(mapping))
    stop("Gene-isoform mapping table does not match to isoform expression table")
  
  times0 <- unique(times)
  if(!is.numeric(times))
    times <- as.numeric(factor(times,levels = unique(times)))
  
  start.time <- Sys.time()
  ##Step 1: Checking data information
  ##genes more than 2 transcripts
  genes0<-as.vector(mapping[,1])
  isoforms0<-as.vector(mapping[,2])
  
  ##expression ratio aboundance
  data2intersect.ratio<-apply(data.exp,1,mean)
  data2intersect.ratio<-rowratio(x = data2intersect.ratio,group = genes0)
  
  #choose genes with more than 2 transcripts
  idx<-table(genes0)[unique(genes0)]
  genes<-names(idx)[idx>1]
  isoforms<-isoforms0[which(genes0 %in% genes)]
  data.exp<-data.exp[isoforms,]
  if(rank)
    data2switch<-apply(data.exp,2,rank) else data2switch=data.exp
  
  
  
  message('Summarising input information ... ')
  message(' Input genes: ', length(unique(genes0)))
  message(' Genes with more than 2 isoforms: ', length(genes))
  message(' Average isoforms per gene for switch analysis: ', round(length(isoforms)/length(genes),3))
  
  
  
  ##step 2: Pickign isoform-pairs with intersection points...
  message(paste0('Step 1: Search for intersection points with ',if(spline) 'Spline method...' else 'Mean expression...'))
  ##Average values of time points
  ##data for searching isoform swtich points
  if(spline){
    message(' Spline fitting expression ...')
    if(is.null(spline.df))
      spline.df<-floor(2*(length(unique(times)-1))/3)
    data2intersect<-t(apply(data.exp,1,function(x) ts.spline(x,times  = times,df = spline.df)))
  } else {
    data2intersect<-t(rowmean(t(data.exp),group = times))
  }
  
  
  ##Find intersection points
  message(' Searching ...')
  iso.intersections<-list()
  t.start<-min(times)
  t.end<-max(times)
  
  withProgress(message = 'Searching intersections: ',value=0,{
    if(verbose)
      pb <- txtProgressBar(min = 0, max = length(genes), style = 3,width = 75)
    for(i in 1:length(genes)){
      
      sub.gene<-genes[i]
      sub.isoform<-as.vector(isoforms0[which(genes0==sub.gene)])
      
      onegene<-data2intersect[sub.isoform,]
      comb<-combn(sub.isoform,2)
      # iso.intersections<-list()
      for(j in 1:ncol(comb)){
        iso1<-comb[1,j]
        iso2<-comb[2,j]
        ##
        x1=onegene[iso1,]
        x2=onegene[iso2,]
        
        if(max(abs(x1-x2))<min.difference)
          next
        
        ##intersection points
        iso.inter <- ts.intersection(x1 = x1, x2 = x2,times = times)
        if(is.null(iso.inter))
          next
        iso.inter<-iso.inter[iso.inter$x.points<t.end & iso.inter$x.points>t.start,]
        ##check if have intersection points
        if(nrow(iso.inter)==0)
          next
        
        ##check if the switching lasting more than 2 consecutive time points
        check.consecutive<-c(t.start,iso.inter$x.points,t.end)
        check.consecutive<-data.frame(low=check.consecutive[-length(check.consecutive)],up=check.consecutive[-1])
        check.consecutive$low<-ceiling(check.consecutive$low)
        check.consecutive$up<-floor(check.consecutive$up)
        check.consecutive<-check.consecutive$up-check.consecutive$low+1
        
        if(max(check.consecutive)<min.t.points)
          next
        
        iso.inter<-data.frame(iso1=iso1,iso2=iso2,t(iso.inter))
        colnames(iso.inter)<-c('iso1','iso2',paste0('switch',1:(ncol(iso.inter)-2)))
        
        iso.intersections<-c(iso.intersections,setNames(object = list(iso.inter),nm = paste0(iso1,'_vs_',iso2)))
      }
      
      incProgress(1/length(genes), detail = paste(i, ' of ', length(genes)))
      Sys.sleep(0)
      if(verbose)
        setTxtProgressBar(pb, i)
    }
    if(verbose)
      close(pb)
  })
  
  if(length(iso.intersections)==0)
    return(NULL)
  message(' ',paste0(length(iso.intersections),' pairs of isoforms have intersection points.'))
  
  
  
  message('Step 2: Calculate scores for isoform switch ...')
  message(' Score 1: Switch frequencies/probabilities')
  message(' Score 2: Sum of average sample differences before and after switch.')
  message(' Score 3: P-values of sample differences before and after switch')
  message(' Score 4: Time points in each intervals')
  message(' Score 5: Pearson correlation of isoforms')
  ##step 2: calculate probablities
  
  iso.names<-names(iso.intersections)
  
  iso.scores<-data.frame()
  
  withProgress(message = 'Scoring isoforms: ',value=0,{
    if(verbose)
      pb <- txtProgressBar(min = 0, max = length(iso.names), style = 3,width = 75)
    for(i in 1:length(iso.names)){
      
      #  i=1
      # i=which(grepl('AT1G01060.2_vs_AT1G01060.3',x = iso.names))
      iso<-unlist(strsplit(iso.names[i],'_vs_'))
      iso1<-iso[1]
      iso2<-iso[2]
      
      # message(iso)
      #generate dataset for iso for isoform switch analysis
      n.inters<-data.frame(iso.intersections[[iso.names[i]]][,-c(1:2)])
      colnames(n.inters)<-colnames(iso.intersections[[iso.names[i]]])[-c(1:2)]
      inters=as.numeric(n.inters[1,])
      
      interval=cut(times,unique(c(t.start,as.numeric(inters),t.end)),include.lowest = T,dig.lab = 5)
      # sub.data2swith=data.frame(interval=interval,t(data.exp.rank[c(iso1,iso2),]),iso.idx=c(-1,1)[as.numeric(interval)%%2+1])
      sub.data2swith=data.frame(interval=interval,t(data2switch[c(iso1,iso2),]))
      ##add a column of sign of differences
      
      diff.sign<-sign(sub.data2swith[,2]-sub.data2swith[,3])
      interval.idx<-as.numeric(interval) %%2
      interval.idx[interval.idx==0]<--1
      diff.sign<-diff.sign*interval.idx
      sub.data2swith$diff.sign<-diff.sign
      x0<-diff.sign[diff.sign!=0][1]
      #Define scores
      
      
      s<-by(data = sub.data2swith,INDICES =sub.data2swith[,1],FUN = function(x){
        ##score1:
        prob<-length(x[,4][x[,4]==x0])/length(x[,4])
        prob[is.na(prob)]<-0
        ##score2:
        x.diff<-x[,2]-x[,3]
        diff.mean<-mean(x.diff)
        # ##score3:
        
        pval <- tryCatch({
          t.test(x[,3],x[,2],paired = T)$p.value
        },error=function(e) return(0))
        pval[is.na(pval)]<-1
        
        #score4:
        # inter.length<-nrow(x)/nrep
        
        y<-as.character(x$interval[1])
        y<-as.numeric(unlist(strsplit(substr(x = y,start = 2,stop = nchar(y)-1),',',fixed = T)))
        inter.length <- length(which(unique(times)>=ceiling(min(y)) & unique(times)<= floor(max(y))))
        x<-data.frame(prob=prob,diff=diff.mean,pval=pval,inter.length=inter.length)
        x
      },simplify = T)
      s<-t(do.call(rbind,s))
      
      if(ncol(s)<3)
        lf.idx<-1:2 else lf.idx<-c(1,rep(2:(ncol(s)-1),each=2),ncol(s))
      
      ##score 1: prob
      score1.2side<-matrix(s[1,lf.idx],ncol=2,byrow = T)
      colnames(score1.2side)<-c('before.prob','after.prob')
      score1<-zoo::rollapply(as.numeric(s[1,]), width = 2, FUN = function(x) abs(sum(x)-1))
      
      ##score 2: difference
      score2<-abs(diff(s[2,]))
      score2.2side<-matrix(s[2,lf.idx],ncol=2,byrow = T)
      colnames(score2.2side)<-c('before.diff','after.diff')
      
      ##score 3: p-value
      score3<-matrix(s[3,lf.idx],ncol=2,byrow = T)
      colnames(score3)<-c('before.pval','after.pval')
      
      ##score3.adj: adjusted p-value
      score3.adj <- p.adjust(p = s[3,lf.idx],method = 'BH')
      score3.adj<-matrix(score3.adj,ncol=2,byrow = T)
      colnames(score3.adj)<-c('before.FDR','after.FDR')
      
      ##score 4: interval length
      score4<-matrix(s[4,lf.idx],ncol=2,byrow = T)
      colnames(score4)<-c('before.t.points','after.t.points')
      
      
      ###average ratio
      iso.mean.ratio<-data.frame(iso1.mean.ratio=data2intersect.ratio[iso1],iso2.mean.ratio=data2intersect.ratio[iso2])
      
      ###before and after intervals
      inter.lr<-matrix(unique(interval)[lf.idx],ncol=2,byrow = T)
      colnames(inter.lr)<-c('before.interval','after.interval')
      
      
      score<-data.frame(iso1=iso1,
                        iso2=iso2,
                        iso.mean.ratio,
                        inter.lr,
                        x.value=as.numeric(n.inters[1,]),
                        y.value=as.numeric(n.inters[2,]),
                        score1.2side,
                        score2.2side,
                        score3,
                        score3.adj,
                        score4,
                        prob=score1,
                        diff=score2,
                        cor=suppressWarnings(cor(as.numeric(data.exp[iso1,]),as.numeric(data.exp[iso2,]))),
                        row.names = NULL)
      
      iso.scores<-rbind(score,iso.scores)
      
      incProgress(1/length(iso.names), detail = paste(i, ' of ', length(iso.names)))
      Sys.sleep(0)
      if(verbose)
        setTxtProgressBar(pb, i)
    }
    if(verbose)
      close(pb)
  })
  ##sort
  iso.scores<-iso.scores[order(iso.scores$prob,decreasing = T),]
  iso.scores <- data.frame(iso.scores,row.names = NULL)
  if(!is.numeric(times0)){
    iso.scores$numeric.coordinates <- NA
    iso.scores$numeric.coordinates[1] <- paste0(paste0(times0,': ',1:length(times0)),collapse = '; ')
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste0('Time for analysis: ',round(time.taken,3),' ',attributes(time.taken)$units))
  
  message('Done!!! ')
  return(iso.scores)
}

#' @export
#' @rdname iso.switch
iso.switch.pairwise <- function(data.exp=TSIS.tpm,
                                mapping = TSIS.mapping,
                                times=times,
                                rank=F,
                                min.difference = 1,
                                spline =F,
                                spline.df = 9,
                                verbose = T,shiny=F){
  if(shiny){
    scores<-iso.switch.shiny(data.exp = data.exp,
                             mapping = mapping,
                             times = times,
                             rank = rank,
                             min.t.points = 1,
                             min.difference = min.difference,
                             spline = spline,
                             spline.df = spline.df,
                             verbose = verbose)
  } else {
    scores<-iso.switch(data.exp = data.exp,
                       mapping = mapping,
                       times = times,
                       rank = rank,
                       min.t.points = 1,
                       min.difference = min.difference,
                       spline = spline,
                       spline.df = spline.df,
                       verbose = verbose)
  }
  if(is.null(scores))
    return(NULL)
  if(!is.numeric(times)){
    time.idx <- unique(times)
    time.label <- data.frame(rep(time.idx[1],nrow(scores)),rep(time.idx[2],nrow(scores)))
    colnames(time.label) <- c('Condition1(x=1)','Condition2(x=2)')
    scores <- cbind(scores[,1:4],time.label,scores[,-c(1:4)])
  }
  scores
}