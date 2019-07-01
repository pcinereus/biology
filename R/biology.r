###########################################################################
#                                                                         #
#  biology.r is a set of functions written by Murray Logan for biologists #
#  and ecologist                                                          #
#                                                                         #
#                                                                         #
#                                                                         #
#                                                                         #
#                                                                         #
#                                                                         #
#                                                                         #
###########################################################################


is.balanced<-function(form, data){
  !is.list(replications(form, data))
}

odds.ratio<-function(object){
    ss <-summary(object) 
    mm<-cbind(exp(ss$coef[-1,"Estimate"])
              ,exp(ss$coef[-1,"Estimate"]-(qnorm(0.5 * (1 + .95))*ss$coef[-1,"Std. Error"]))
              ,exp(ss$coef[-1,"Estimate"]+(qnorm(0.5 * (1 + .95))*ss$coef[-1,"Std. Error"]))
              )
    rownames(mm) <- dimnames(ss$coef)[[1]][-1]
    colnames(mm)<-c("Odds ratio", "Lower 95% CI","Upper 95% CI")
    mm
  }

oddsratios<-function(tab, rev="neither", correct=FALSE){
require(epitools)
  or <- data.frame(matrix(1,choose(length(rownames(tab)),2),7))
  #rnames<-vector(mode="character",length=nrow(tab))
  k<-0
  for (i in 1:(nrow(tab)-1)) {
    for (j in (i+1):nrow(tab)) {
      if (correct==TRUE) {
        or[k<-k+1,2:4]<-(oddsratio.wald(tab[c(i,j),]+.5, rev=rev)$measure[2,]) #added 0.5 to deal with 0's
        or[k,5:7]<-(oddsratio.wald(tab[c(i,j),]+.5, rev=rev)$p.value[2,])
      }
      if (correct==FALSE) {
        or[k<-k+1,2:4]<-(oddsratio.wald(tab[c(i,j),], rev=rev)$measure[2,]) #added 0.5 to deal with 0's
        or[k,5:7]<-(oddsratio.wald(tab[c(i,j),], rev=rev)$p.value[2,])
      }
                                        #rnames[k]<-paste(rownames(tab[c(i,j),]), collapse=" vs ")
      or[k,1]<-paste(rownames(oddsratio.wald(tab[c(i,j),], rev=rev)$measure), collapse=" vs ") 
    }
  }
  #rownames(or)<-rnames
  colnames(or)<-c("Comparison",colnames(oddsratio.wald(tab[c(1,2),], rev=rev)$measure),colnames(oddsratio.wald(tab[c(1,2),], rev=rev)$p.value))
  or[,1]<-as.factor(or[,1])
  or
#  as.data.frame(or)
}

#Pete Hurd's version of a g.test 
      g.test <- function(x, y = NULL, correct="none",
                         p = rep(1/length(x), length(x)), simulate.p.value = FALSE, B = 2000)
        {
          DNAME <- deparse(substitute(x))
          if (is.data.frame(x)) x <- as.matrix(x)
          if (is.matrix(x)) {
            if (min(dim(x)) == 1) 
              x <- as.vector(x)
          }
          if (!is.matrix(x) && !is.null(y)) {
            if (length(x) != length(y)) 
              stop("x and y must have the same length")
            DNAME <- paste(DNAME, "and", deparse(substitute(y)))
            OK <- complete.cases(x, y)
            x <- as.factor(x[OK])
            y <- as.factor(y[OK])
            if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
              stop("x and y must have at least 2 levels")
            x <- table(x, y)
          }
          if (any(x < 0) || any(is.na(x))) 
            stop("all entries of x must be nonnegative and finite")
          if ((n <- sum(x)) == 0) 
            stop("at least one entry of x must be positive")
                                        #If x is matrix, do test of independence
          if (is.matrix(x)) {
                                        #Test of Independence
            nrows<-nrow(x)
            ncols<-ncol(x)
            if (correct=="yates"){ # Do Yates' correction?
              if(dim(x)[1]!=2 || dim(x)[2]!=2) # check for 2x2 matrix
                stop("Yates' correction requires a 2 x 2 matrix")
              if((x[1,1]*x[2,2])-(x[1,2]*x[2,1]) > 0)
                {
                  x[1,1] <- x[1,1] - 0.5
                  x[2,2] <- x[2,2] - 0.5
                  x[1,2] <- x[1,2] + 0.5
                  x[2,1] <- x[2,1] + 0.5
                }
              else
                {
                  x[1,1] <- x[1,1] + 0.5
                  x[2,2] <- x[2,2] + 0.5
                  x[1,2] <- x[1,2] - 0.5
                  x[2,1] <- x[2,1] - 0.5
                }
            }
            
            sr <- apply(x,1,sum)
            sc <- apply(x,2,sum)
            E <- outer(sr,sc, "*")/n
                                        # are we doing a monte-carlo?
                                        # no monte carlo GOF?
            if (simulate.p.value){
              METHOD <- paste("Log likelihood ratio (G-test) test of independence\n\t with simulated p-value based on", B, "replicates")
              tmp <- .C("gtestsim", as.integer(nrows), as.integer(ncols),
                        as.integer(sr), as.integer(sc), as.integer(n), as.integer(B),
                        as.double(E), integer(nrows * ncols), double(n+1),
                        integer(ncols), results=double(B), PACKAGE= "ctest")
              g <- 0
              for (i in 1:nrows){
                for (j in 1:ncols){
                  if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
                }
              }
              STATISTIC <- G <- 2 * g
              PARAMETER <- NA
              PVAL <- sum(tmp$results >= STATISTIC)/B
            }
            else {
                                        # no monte-carlo
                                        # calculate G
              g <- 0
              for (i in 1:nrows){
                for (j in 1:ncols){
                  if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
                }
              }
              q <- 1
              if (correct=="williams"){ # Do Williams' correction
                row.tot <- col.tot <- 0    
                for (i in 1:nrows){ row.tot <- row.tot + 1/(sum(x[i,])) }
                for (j in 1:ncols){ col.tot <- col.tot + 1/(sum(x[,j])) }
                q <- 1+ ((n*row.tot-1)*(n*col.tot-1))/(6*n*(ncols-1)*(nrows-1))
              }
              STATISTIC <- G <- 2 * g / q
              PARAMETER <- (nrow(x)-1)*(ncol(x)-1)
              PVAL <- 1-pchisq(STATISTIC,df=PARAMETER)
              if(correct=="none")
                METHOD <- "Log likelihood ratio (G-test) test of independence without correction"
              if(correct=="williams")
                METHOD <- "Log likelihood ratio (G-test) test of independence with Williams' correction"
              if(correct=="yates")
                METHOD <- "Log likelihood ratio (G-test) test of independence with Yates' correction"
            }
          }
          else {
                                        # x is not a matrix, so we do Goodness of Fit
            METHOD <- "Log likelihood ratio (G-test) goodness of fit test"
            if (length(x) == 1) 
              stop("x must at least have 2 elements")
            if (length(x) != length(p)) 
              stop("x and p must have the same number of elements")
            E <- n * p
            
            if (correct=="yates"){ # Do Yates' correction
              if(length(x)!=2)
                stop("Yates' correction requires 2 data values")
              if ( (x[1]-E[1]) > 0.25) {
                x[1] <- x[1]-0.5
                x[2] <- x[2]+0.5
              }
              else if ( (E[1]-x[1]) > 0.25){
                x[1] <- x[1]+0.5
                x[2] <- x[2]-0.5
              }
            }
            names(E) <- names(x)
            g <- 0
            for (i in 1:length(x)){
              if (x[i] != 0) g <- g + x[i] * log(x[i]/E[i])
            }
            q <- 1
            if (correct=="williams"){ # Do Williams' correction
              q <- 1+(length(x)+1)/(6*n)
            }
            STATISTIC <- G <- 2*g/q
            PARAMETER <- length(x) - 1
            PVAL <- pchisq(STATISTIC, PARAMETER, lower = FALSE)
          }
          names(STATISTIC) <- "Log likelihood ratio statistic (G)"
          names(PARAMETER) <- "X-squared df"
          names(PVAL) <- "p.value"
          structure(list(statistic=STATISTIC,parameter=PARAMETER,p.value=PVAL,
                         method=METHOD,data.name=DNAME, observed=x, expected=E),
                    class="htest")
        }

##The following are used for the Johnson-Neyman Wilcox procedure
##1. JN - is an internal
##2. wilcox.JN - is the external
JN <-  function(N1,M1,SSx1,b01,b11,SEres1, N2,M2,SSx2,b02,b12,SEres2,tmts, type="H") {
    ## calculate the MSresiduals
    MSres1<-SEres1/(N1-2)
    MSres2<-SEres2/(N2-2)

    ## determine the degrees of freedom
    FNcalcD<-function(samplesize, covmean, covval, covss){
      ##calculates the value of 'd' in the equations
      (1/samplesize)+(((covmean-covval)^2)/covss)
    }
    d1<-FNcalcD(N1,M1,M1,SSx1)
    d2<-FNcalcD(N2,M2,M2,SSx2)
    MSres1<-SEres1/(N1-2)
    MSres2<-SEres2/(N2-2)
    mu1=MSres1*d1 
    mu2=MSres2*d2
    (nu<-as.integer((((mu1+mu2)^2)/((mu1^2/(N1-2))+(mu2^2/(N2-2))))))

    interpval<-function(lowdf,lowval,highdf,highval, compdf){
      while (1==1){
        newdf<-((highdf-lowdf)/2)+lowdf
        newval<-1/(0.5*((1/lowval)+(1/highval)))
        if (abs(newdf-compdf)<0.1) break 
        if (newdf>compdf) {
          highdf<-newdf
          highval<-newval
        }
        if (newdf<compdf) {
          lowdf<-newdf
          lowval<-newval
        }
      }
      newval
    }

    critvalue<-function(tmatrix, tmts, nu){
      Hmatrix <-
        structure(c(4.38, 4.13, 3.97, 3.85, 3.76, 3.69, 3.64, 3.59, 3.56, 
                    3.53, 3.49, 3.48, 3.46, 3.44, 3.43, 3.41, 3.37, 3.33, 3.28, 3.23, 
                    3.21, 5.38, 5.03, 4.79, 4.64, 4.52, 4.42, 4.35, 4.29, 4.24, 4.19, 
                    4.16, 4.12, 4.09, 4.07, 4.05, 4.03, 3.97, 3.91, 3.85, 3.78, 3.74, 
                    6.01, 5.59, 5.33, 5.13, 4.99, 4.88, 4.79, 4.71, 4.65, 4.59, 4.55, 
                    4.52, 4.48, 4.45, 4.43, 4.39, 4.33, 4.26, 4.18, 4.09, 4.05, 6.47, 
                    6.01, 5.71, 5.49, 5.33, 5.19, 5.09, 5.02, 4.95, 4.89, 4.84, 4.79, 
                    4.76, 4.73, 4.69, 4.67, 4.59, 4.51, 4.42, 4.33, 4.27, 6.83, 6.34, 
                    6.01, 5.77, 5.59, 5.46, 5.35, 5.26, 5.19, 5.12, 5.07, 5.02, 4.98, 
                    4.94, 4.91, 4.88, 4.79, 4.69, 4.61, 4.51, 4.44, 7.12, 6.59, 6.25, 
                    5.99, 5.82, 5.67, 5.55, 5.46, 5.38, 5.31, 5.25, 5.19, 5.16, 5.12, 
                    5.08, 5.05, 4.96, 4.86, 4.76, 4.66, 4.58, 7.37, 6.83, 6.46, 6.19, 
                    5.99, 5.85, 5.73, 5.63, 5.54, 5.47, 5.41, 5.36, 5.31, 5.27, 5.23, 
                    5.19, 5.09, 4.99, 4.89, 4.78, 4.69),
                  .Dim = c(21L, 7L),
                  .Dimnames = list(
                    c("5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
                      "15", "16", "17", "18", "19", "20", "24", "30", "40", "60", 
                      "120"), c("2", "3", "4", "5", "6", "7", "8")))
      ##  matrix of H values - 1st line is data for matrix size
      ## 	1st dimension of matrix is rows - degrees of freedom
      ## 	2nd			is cols - number of treatments
      ## 	3rd				- probability - currently 0.05 only
      SMMmatrix <-
        structure(c(5.57, 3.96, 3.38, 3.09, 2.92, 2.8, 2.72, 2.66, 2.61, 
                    2.57, 2.54, 2.49, 2.46, 2.43, 2.41, 2.38, 2.35, 2.32, 2.29, 2.24, 
                    6.34, 4.43, 3.74, 3.4, 3.19, 3.06, 2.96, 2.89, 2.83, 2.78, 2.75, 
                    2.69, 2.65, 2.62, 2.59, 2.56, 2.52, 2.49, 2.45, 2.39, 6.89, 4.76, 
                    4, 3.62, 3.39, 3.24, 3.13, 3.05, 2.98, 2.93, 2.89, 2.83, 2.78, 
                    2.75, 2.72, 2.68, 2.64, 2.6, 2.56, 2.49, 7.31, 5.02, 4.2, 3.79, 
                    3.54, 3.38, 3.26, 3.17, 3.1, 3.05, 3, 2.94, 2.89, 2.85, 2.82, 
                    2.77, 2.73, 2.69, 2.65, 2.57, 7.65, 5.23, 4.37, 3.93, 3.66, 3.49, 
                    3.36, 3.27, 3.2, 3.14, 3.09, 3.02, 2.97, 2.93, 2.9, 2.85, 2.8, 
                    2.76, 2.72, 2.63),
                  .Dim = c(20L, 5L),
                  .Dimnames = list(c("2", 
                    "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "14", "16", 
                    "18", "20", "24", "30", "40", "60", "10000"), c("2", "3", "4", 
                                                                    "5", "6")))
      ##  matrix of SMM values - 1st line is data for matrix size
      ## 	1st dimension of matrix is rows - degrees of freedom
      ## 	2nd			is cols - number of contrasts
      ## 	3rd				- probability - currently 0.05 only
      
      ## need to put in a catch to get number of categories greater than are in the H or SMM table
      if (tmatrix=="Hmatrix") hmatrix<-Hmatrix
      if (tmatrix=="SMMmatrix")  hmatrix<-SMMmatrix
      lowval<-0
      lowdof<-0
      highval<-0
      highdof<-0
      for (i in 1: nrow(hmatrix)){
                                        #print(as.numeric(rownames(hmatrix)[i]))
                                        #print(nu)
        if (as.numeric(rownames(hmatrix)[i])<nu){
          lowdof<-as.numeric(rownames(hmatrix)[i])
          lowval<-hmatrix[i,as.character(tmts)]
        }
        if (as.numeric(rownames(hmatrix)[i])>nu){
          highdof<-as.numeric(rownames(hmatrix)[i])
          highval<-hmatrix[i,as.character(tmts)]
          value<-interpval(lowdof,lowval,highdof,highval,nu)
          break
        }
        if (as.numeric(rownames(hmatrix)[i])==nu) {
          value<-hmatrix[i,as.character(tmts)]
          break
        }
      }
      class(value)<-"H"
      value
    }
    ## calculate the critical h value
    if (type=="H") {
      h<-critvalue("Hmatrix",tmts,nu)
    }
    else {
      h<-(2^.5)*critvalue("SMMmatrix",tmts,nu)
    }
    
    #A<-((b11-b12)^2)-((h^2)/2)*((MSres1/SSx1)+(MSres2/SSx2))
    A<-((h^2)/2)*((MSres1/SSx1)+(MSres2/SSx2))
    A<-((b11-b12)^2)-A
    B<-(h^2)*((MSres1*M1/SSx1)+(MSres2*M2/SSx2))
	B<-B+(2*((b01-b02)*(b11-b12)))
    #B<-2*((b11-b12)*(b01-b02))
    #B2<-MSres1*((M1)/SSx1)
    #B3<-MSres2*((M2)/SSx2)
    #B<-B+((h^2)*(B2+B3))
    C<- -((h^2)/2)*((MSres1/N1)+(MSres2/N2))
		firsthalf<-MSres1*((M1^2)/SSx1)
		secondhalf<-MSres2*((M2^2)/SSx2)
	    C<-C-(((h^2)/2)*(firsthalf+secondhalf))
		C<-C+((b01-b02)^2)

    #E <-((h^2)/2)*((MSres1/N1)+(MSres2/N2))
    #C <-(b01-b02)^2
    #C2<-MSres1*((M1^2)/SSx1)
    #C3<-MSres2*((M2^2)/SSx2)
    #C<-C-E
    #C<-C-((h^2/2)*(C2+C3))
    XLOWER<-(-B-sqrt((B^2)-(4*A*C)))/(2*A)
    XUPPER<-(-B+sqrt((B^2)-(4*A*C)))/(2*A)
    
    c(nu,h,XLOWER,XUPPER)
  }


## wilcox.JN<-function(object, type="H"){
##   ## determine how many treatment levels there are
##   #whichfn<-#which(attr(constable.lm$terms,"dataClasses")=="factor")
##   #whichf<-#names(which(attr(constable.lm$terms,"dataClasses")=="factor"))
##   object
##   NumLevels <- length(unlist(object$xlevels)) #length(object$xlevels[[1]])
##   combN <- choose(NumLevels,2)
##   combS <- combn(unlist(object$xlevels),2)
##   tab<-matrix(1,nrow=combN,ncol=4)
##   rownames(tab)<-apply(combS,2,paste, collapse=" vs ")
##   colnames(tab)<-c("df", "critical value", "lower","upper")
##   for (i in 1:combN) {
##     LEV1 <- combS[1,i]
##     N1<-length(object$model[2][object$model[3]==LEV1])
##     M1<-mean(object$model[2][object$model[3]==LEV1])
##     SSx1<-var(object$model[2][object$model[3]==LEV1])*(N1-1)
##     #Update model
##     paste("update(object, ~object$model[[2]], subset=",names(object$xlevels),"=='",LEV1,"')", sep="")
##     lmm <- eval(parse(text=paste("update(object, ~object$model[[2]], subset=",names(object$xlevels),"=='",LEV1,"')", sep="")))
##     b01<-summary(lmm)$coef[1,1]
##     b11<-summary(lmm)$coef[2,1]
##     SEres1 <- anova(lmm)[2,2]
  
##     LEV2 <- combS[2,i]
##     N2<-length(object$model[2][object$model[3]==LEV2])
##     M2<-mean(object$model[2][object$model[3]==LEV2])
##     SSx2<-var(object$model[2][object$model[3]==LEV2])*(N2-1)
##     lmm <- eval(parse(text=paste("update(object, ~object$model[[2]], subset=",names(object$xlevels),"=='",LEV2,"')", sep="")))
##     b02<-summary(lmm)$coef[1,1]
##     b12<-summary(lmm)$coef[2,1]
##     SEres2 <- anova(lmm)[2,2]
##     tab[i,]<-(JN(N1,M1,SSx1,b01,b11,SEres1,N2,M2,SSx2,b02,b12,SEres2,NumLevels, type))
##   }
##   tab
## }
##################################
wilcox.JN<-function(object, type="H"){
  NumLevels <- length(unlist(object$xlevels))
  combN <- choose(NumLevels,2)
  combS <- combn(unlist(object$xlevels),2)
  tab<-matrix(1,nrow=combN,ncol=4)
  rownames(tab)<-apply(combS,2,paste, collapse=" vs ")
  colnames(tab)<-c("df", "critical value", "lower","upper")
  #determine the covariate
  covName <- attr(object$terms, "dataClasses")[-attr(object$terms, "response")]
  covName <-names(covName[covName!="factor"])
  for (i in 1:combN) {
    LEV1 <- combS[1,i]
    N1<-length(object$model[2][object$model[3]==LEV1])
    M1<-mean(object$model[2][object$model[3]==LEV1])
    SSx1<-var(object$model[2][object$model[3]==LEV1])*(N1-1)
    lmm <- eval(parse(text=paste("update(object, ~",covName,", subset=",names(object$xlevels),"=='",LEV1,"')", sep="")))
    b01<-summary(lmm)$coef[1,1]
    b11<-summary(lmm)$coef[2,1]
    SEres1 <- anova(lmm)[2,2]
    LEV2 <- combS[2,i]
    N2<-length(object$model[2][object$model[3]==LEV2])
    M2<-mean(object$model[2][object$model[3]==LEV2])
    SSx2<-var(object$model[2][object$model[3]==LEV2])*(N2-1)
    lmm <- eval(parse(text=paste("update(object, ~",covName,", subset=",names(object$xlevels),"=='",LEV2,"')", sep="")))
    b02<-summary(lmm)$coef[1,1]
    b12<-summary(lmm)$coef[2,1]
    SEres2 <- anova(lmm)[2,2]
    tab[i,]<-(JN(N1,M1,SSx1,b01,b11,SEres1,N2,M2,SSx2,b02,b12,SEres2,NumLevels, type))
  }
 tab 
}



## epsi.GG.HF<-function (object) {
##   require(car)
##   gg.hf <- NULL
##   call <- attr(object, "call")
##   formula <- formula(attr(object,"call"))
##   Terms <- attr(object, "terms")
##   indError <- attr(Terms, "specials")$Error
##                                         #put something in here incase there are no error strata
##                                         #
##   allTerms <- Terms
##                                         #Determine the error term   
##   errorterm <- attr(Terms, "variables")[[1 + indError]]
##   eTerm <- deparse(errorterm[[2]], width.cutoff = 500, 
##                    backtick = TRUE)
##                                         #Response variable
##   resp <- as.character(formula[[2]])   
##   block <- names(object)[2]
##                                         #Within subject effects   
##   rept <- row.names(summary(object[[3]])[[1]])
##   rept_lab <- trim.blanks(rept[-grep("Residuals", rept)])
##   rept <- trim.blanks(rept[-grep(":|Residuals", rept)])
##   rept <- sub("C\\(([A-Za-z_\\.]*).*\\)", "\\1", rept)
##   drop<-NULL
##   data <- eval(call[[3]])
##   eval(parse(text=paste("data$",paste(rept,collapse=""),"<-with(data, interaction(", paste(rept,collapse=","),"))",sep="")))
##   if (length(rept)>1) {drop<-rept}
##                                         #idata <- eval(parse(text=paste("expand.grid(",paste(drop,"=levels(data$",drop,")", sep="",collapse=","),")", sep="")))
##   idata <- eval(parse(text=paste("expand.grid(",paste(rept,"=levels(data$",rept,")", sep="",collapse=","),")", sep="")))
##   within <- paste(rept,collapse="*")
##   rept<-paste(rept,collapse="")
##   between <- row.names(summary(object[[2]])[[1]])
##   between_lab <- trim.blanks(between[-grep("Residuals", between)])
##   between <- trim.blanks(between[-grep(":|Residuals", between)])
##                                         #Reshape the data for MANOVA
##   data.rm <- reshape(data, direction = "wide", timevar = rept, 
##                      idvar = block, v.names = resp, drop=drop)
##   if (length(between_lab) == 0) {
##     data.manova <- lm(as.matrix(data.rm[grep(resp, colnames(data.rm))]) ~ 
##                       1)
##   }
##   if(length(between_lab) > 0){
##     between <- eval(parse(text = paste("data.rm$", between, 
##                             sep = "")))
##     data.manova <- lm(as.matrix(data.rm[grep(resp, colnames(data.rm))]) ~ 
##                       between)
##   }
##   data.ok <- eval(parse(text=paste("Anova(data.manova, idata = idata, idesign =~",within,", 
##         type = 'III')",sep="")))
##   gg.hf <- GGHF(data.ok)
##   colnames(gg.hf) <- c("Greenhouse-Geisser", "Hunyh-Feldt")
##   rownames(gg.hf) <- c(between_lab, rept_lab)
##   cat("Epsilon values\n")
##   gg.hf
## }


## epsi.GG.HF<-function (object) 
## {
##     require(car)
##     gg.hf <- NULL
##     call <- attr(object, "call")
##  Terms <- attr(object, "terms")
##  indError <- attr(Terms, "specials")$Error
##  #put something in here incase there are no error strata
##  #
##  allTerms <- Terms
##  #Determine the error term   
##  errorterm <- attr(Terms, "variables")[[1 + indError]]
##  eTerm <- deparse(errorterm[[2]], width.cutoff = 500, 
##             backtick = TRUE)
##  #Response variable
##  resp <- formula[[2]]   
##     block <- names(object)[2]
##     #resp <- as.character(call[[2]][[2]])
##  #Within subject effects   
##     rept <- row.names(summary(object[[3]])[[1]])
##     rept_lab <- trim.blanks(rept[-grep("Residuals", rept)])
##     rept <- trim.blanks(rept[-grep(":|Residuals", rept)])
##     rept <- sub("C\\(([A-Za-z_\\.]*).*\\)", "\\1", rept)
##     data <- eval(call[[3]])
##     between <- row.names(summary(object[[2]])[[1]])
##     between_lab <- trim.blanks(between[-grep("Residuals", between)])
##     between <- trim.blanks(between[-grep(":|Residuals", between)])
## #Reshape the data for MANOVA
##     data.rm <- reshape(data, direction = "wide", timevar = rept, 
##         idvar = block, v.names = resp)
##     if (length(between_lab) == 0) {
##         data.manova <- lm(as.matrix(data.rm[grep(resp, colnames(data.rm))]) ~ 
##             1)
##     }
##     else {
##         between <- eval(parse(text = paste("data.rm$", between, 
##             sep = "")))
##         data.manova <- lm(as.matrix(data.rm[grep(resp, colnames(data.rm))]) ~ 
##             between)
##     }
##     idata <- data.frame(within = as.factor(levels(eval(parse(text = paste(as.character(call[[3]]), 
##         "$", rept, sep = ""))))))
##     data.ok <- Anova(data.manova, idata = idata, idesign = ~within, 
##         type = "III")
##     gg.hf <- GGHF(data.ok)
##     colnames(gg.hf) <- c("Greenhouse-Geisser", "Hunyh-Feldt")
##     rownames(gg.hf) <- c(between_lab, rept_lab)
##     cat("Epsilon values\n")
##     gg.hf
## }

epsi.GG.HF<-function (object) {
  require(car)
  gg.hf <- NULL
  call <- attr(object, "call")
  formula <- formula(attr(object,"call"))
  Terms <- attr(object, "terms")
  indError <- attr(Terms, "specials")$Error
                                        #put something in here incase there are no error strata
                                        #
  allTerms <- Terms
                                        #Determine the error term   
  errorterm <- attr(Terms, "variables")[[1 + indError]]
  eTerm <- deparse(errorterm[[2]], width.cutoff = 500, 
                   backtick = TRUE)
                                        #Response variable
  resp <- as.character(formula[[2]])   
  block <- names(object)[2]
                                        #Within subject effects   
  rept <- row.names(summary(object[[3]])[[1]])
  rept_lab <- trim.blanks(rept[-grep("Residuals", rept)])
  rept <- trim.blanks(rept[-grep(":|Residuals", rept)])
  rept <- sub("C\\(([A-Za-z_\\.]*).*\\)", "\\1", rept)
  drop<-NULL
  data <- eval(call[[3]])
  eval(parse(text=paste("data$",paste(rept,collapse=""),"<-with(data, interaction(", paste(rept,collapse=","),"))",sep="")))
  if (length(rept)>1) {drop<-rept}
                                        #idata <- eval(parse(text=paste("expand.grid(",paste(drop,"=levels(data$",drop,")", sep="",collapse=","),")", sep="")))
  idata <- eval(parse(text=paste("expand.grid(",paste(rept,"=levels(data$",rept,")", sep="",collapse=","),")", sep="")))
  within <- paste(rept,collapse="*")
  rept<-paste(rept,collapse="")
  between <- row.names(summary(object[[2]])[[1]])
  between_lab <- trim.blanks(between[-grep("Residuals", between)])
  between <- trim.blanks(between[-grep("Residuals", between)])
  between <- paste(between,collapse="+")
                                        #Reshape the data for MANOVA
  data.rm <- reshape(data, direction = "wide", timevar = rept, 
                     idvar = block, v.names = resp, drop=drop)
  if (length(between_lab) == 0) {
    #data.manova <- lm(as.matrix(data.rm[grep(resp, colnames(data.rm))]) ~ 
    #                  1)
    epsi <- stats:::sphericity(SSD(lm(as.matrix(data.rm[grep(resp, colnames(data.rm))])~1, data.rm)), X=~1, idata=idata)
  }
  if(length(between_lab) > 0){
    #between <- eval(parse(text = paste("data.rm$", between, 
    #                        sep = "")))
    #data.manova <- lm(as.matrix(data.rm[grep(resp, colnames(data.rm))]) ~ 
    #                  between)
    epsi <- eval(parse(text=paste("stats:::sphericity(SSD(lm(as.matrix(data.rm[grep(resp, colnames(data.rm))])~",between,", data.rm)), X=~1, idata=idata)", sep="")))
#    epsi <- eval(parse(text=paste("stats:::sphericity(SSD(lm(as.matrix(",data.rm[grep(resp, colnames(data.rm))],",)~",between,", data.rm)), X=~1, idata=idata)", sep="")))
  }
#  epsi <- stats:::sphericity(SSD(lm(as.matrix(data.rm[grep(resp, colnames(data.rm))])~between, data.rm)), X=~1, idata=idata)
#  data.ok <- eval(parse(text=paste("Anova(data.manova, idata = idata, idesign =~",within,", 
#        type = 'III')",sep="")))
#  gg.hf <- GGHF(data.ok)
#  colnames(gg.hf) <- c("Greenhouse-Geisser", "Hunyh-Feldt")
#  rownames(gg.hf) <- c(between_lab, rept_lab)
#  cat("Epsilon values\n")
  epsi
}


mainEffects <- function(object, at, ...) {
  #performs simple main effects
  ## object - name of fitted model (aov, aovlist, lm, lme)
  ## at - how the analysis should be split (at=FACTA=="level1")
  ## should be used along with summary
  INT<-calculateInteraction(object, deparse(substitute(at)))
  if(class(object)=="mer") data_name<-attr(object, "call")[[3]] #lmer
  else {
    if(!is.null(object$call)) data_name<-object$call[[3]] #aov or lme
    else data_name<-attr(object, "call")[[3]]            #aovlist
  }
  dn<<-eval(data_name) #this has a <<- so as to use the global scope
  dn<<-data.frame(dn, INT) #it would be more elegant to modify the actual data set than make a new global one that is not linked to the original
  if (any(class(object)=="aovlist")) {
    form <- update(formula(attr(object,"call")), ~INT+.)
    aa <- paste("aov(",form[2],"~",form[3],",dn)")
    command<-aa
    #object1<<-aov(update(formula(attr(object,"call")), ~INT+.), dn)
    #command<-paste("object1")
  }
  if (any(class(object)=="aov")) {
    object<-update(object, data=dn)
    command <- paste("aov(update(object,~INT+.),dn)")
  }
  if (any(class(object)=="lm")) {
    object<-update(object, data=dn)
    command <- paste("aov(update(object,~INT+.),dn)")
  }
  if (any(class(object)=="lme")){
    M <- model.matrix(update(formula(object), ~INT+.), dn)
    command<-paste("mainEffects.lme(object, M)")
  }
  if (any(class(object)=="mer")){
    print("Not implimented yet - there is an issue with update and lmer")
  }
  eval(parse(text=command))
}

## calculateInteraction <- function(object, at1){
##   ## helper function for mainEffects
##   ## Basically, this function gets the various interactions from the model
##   ## not to be called
##   at1<-strsplit(at1," *== *")[[1]]
##   at1[2]<-gsub("\"","",at1[2])
##   if(!is.null(object$call)) data_name<-object$call[[3]] #aov
##    else data_name<-attr(object, "call")[[3]]            #aovlist
##   dn<-eval(data_name)
##   #Determine what the interaction involving the at term is
##   if(!is.null(object$terms)) nms<-names(which(attr(object$terms,"factors")[at1[1],]==1))  #aov
##   else nms<-names(which(attr(attr(object,"terms"),"factors")[at1[1],]==1))                #aovlist
##   nms<-nms[length(nms)]
##   intt <- paste("dn$INT<-interaction(dn$",nms,")",sep="")
##   intt <- gsub(":",",dn$",intt)
##   eval(parse(text=intt))
##   #create a list to hold the base level of each factor
##   nms <- strsplit(nms,":")[[1]]
##   eval(parse(text=paste("base<-lapply(nms, function(x) levels(dn[,x])[1])", sep="")))
##   names(base) <- nms
##   base[at1[1]]<-at1[2]
##   intt2 <- paste("dn$INT[dn$",at1[1],"=='",at1[2],"'] <- '",paste(base, collapse="."),"'", sep="")
##   eval(parse(text=intt2))
##   dn$INT
## }

calculateInteraction <- function(object, at1){
  ## helper function for mainEffects
  ## Basically, this function gets the various interactions from the model
  ## not to be called
  at1<-strsplit(at1," *== *")[[1]]
  at1[2]<-gsub("\"","",at1[2])
  if(class(object)=="mer") data_name<-attr(object, "call")[[3]]            #lmer
  else {
    if(!is.null(object$call)) data_name<-object$call[[3]] #aov
    else data_name<-attr(object, "call")[[3]]            #aovlist
  }
  dn<-eval(data_name)
  #Determine what the interaction involving the at term is
  if(class(object)=="mer") nms<-names(which(attr(terms(object),"factors")[at1[1],]==1))                #lmer
  else {
    if(!is.null(object$terms)) nms<-names(which(attr(object$terms,"factors")[at1[1],]==1))  #aov
    else nms<-names(which(attr(attr(object,"terms"),"factors")[at1[1],]==1))                #aovlist
  }
  nms<-nms[length(nms)]
  intt <- paste("dn$INT<-interaction(dn$",nms,")",sep="")
  intt <- gsub(":",",dn$",intt)
  eval(parse(text=intt))
  #create a list to hold the base level of each factor
  nms <- strsplit(nms,":")[[1]]
  eval(parse(text=paste("base<-lapply(nms, function(x) levels(dn[,x])[1])", sep="")))
  names(base) <- nms
  base[at1[1]]<-at1[2]
  intt2 <- paste("dn$INT[dn$",at1[1],"=='",at1[2],"'] <- '",paste(base, collapse="."),"'", sep="")
  eval(parse(text=intt2))
  dn$INT
}

mainEffects.lme<-function(object, M, ...){
  ## This is a version of the mainEffects function that specifically handles lme
  if(!is.null(object$call)) data_name<-object$call[[3]] #aov and lme
  else data_name<-attr(object, "call")[[3]]            #aovlist
  dn<-eval(data_name)
                                        #this step is necessary to make the object subsequently
                                        #update in this environment - otherwise it will update in the
                                        #aov creation environment and this is not the environment
                                        #that contains Ms
  object<-update(object)
  QR <- qr(M, tol = 1e-07, LAPACK = FALSE)
  ass<-attr(QR$qr,"assign")[QR$pivot[seq_len(QR$rank)]]
  M2 <- stats:::Thin.col(M)
  rm(Ms)
  mmm<-""
  tms <-paste(attr(terms(object), "term.labels"),"", sep="")
  for (i in unique(ass)[-1]){
    eval(parse(text=paste("dn$M",i,"<-M2[,ass==i]",sep="")))
    mmm <- paste(mmm,"+","M",i, sep="")
  }
  mm <- paste("anova(update(object, ~",mmm, ", data=dn))", sep="")
  mm <- paste("update(object, ~",mmm, ", data=dn)", sep="")
  tbl<-eval(parse(text=mm))
                                        #rownames(tbl) <- c("Intercept",tms[unique(ass)[-1:-2]-1:-2])
                                        #  rownames(tbl) <- c("(Intercept)","",tms[unique(ass)[-1]-1])
                                        #rownames(tbl) <- c("(Intercept","",tms[unique(ass)[-1]-1])
  tbl
}

mainEffects.lmer<-function(object, M, ...){
  ## This is a version of the mainEffects function that specifically handles lme
  data_name<-attr(object, "call")[[3]]            #aovlist
  dn<-eval(data_name)
                                        #this step is necessary to make the object subsequently
                                        #update in this environment - otherwise it will update in the
                                        #aov creation environment and this is not the environment
                                        #that contains Ms
  object<-update(object)
  QR <- qr(M, tol = 1e-07, LAPACK = FALSE)
  ass<-attr(QR$qr,"assign")[QR$pivot[seq_len(QR$rank)]]
  M2 <- stats:::Thin.col(M)
  rm(Ms)
  mmm<-""
  tms <-paste(attr(terms(object), "term.labels"),"", sep="")
  for (i in unique(ass)[-1]){
    eval(parse(text=paste("dn$M",i,"<-M2[,ass==i]",sep="")))
    mmm <- paste(mmm,"+","M",i, sep="")
  }
  mm <- paste("anova(update(object, ~",mmm, ", data=dn))", sep="")
  mm <- paste("update(object, ~",mmm, ", data=dn)", sep="")
  tbl<-eval(parse(text=mm))
                                        #rownames(tbl) <- c("Intercept",tms[unique(ass)[-1:-2]-1:-2])
                                        #  rownames(tbl) <- c("(Intercept)","",tms[unique(ass)[-1]-1])
                                        #rownames(tbl) <- c("(Intercept","",tms[unique(ass)[-1]-1])
  tbl
}



GGHF <- function (object) {
  ## this function performs the grunt of the work for calculating
  ## Greenhouse-Geisser and HuhnFeldt epsilons
  nterms <- length(object$terms)
  GG <- function(SSPE, P) {
    p <- nrow(SSPE)
    if (p < 2) 
      return(NA)
    lambda <- eigen(SSPE %*% solve(t(P) %*% P))$values
    lambda <- lambda[lambda > 0]
    ((sum(lambda)/p)^2)/(sum(lambda^2)/p)
  }
  HF <- function(gg, error.df, p) {
    ((error.df + 1) * p * gg - 2)/(p * (error.df - p * gg))
  }
  gg.hf <- NULL
  for (term in 2:nterms) {
    SSP <- object$SSP[[term]]
    SSPE <- object$SSPE[[term]]
    P <- object$P[[term]]
    p <- ncol(P)
    gg <- GG(SSPE, P)
    error.df <- object$error.df
    hf <- HF(gg, error.df, p)
    gg.hf <- rbind(gg.hf, cbind(gg, hf))
  }
  gg.hf
}





## epsi.GG.HF <- function (object) 
## {
##     require(car)
##     gg.hf <- NULL
##     call <- attr(object, "call")
##     block <- names(object)[2]
##     resp <- as.character(call[[2]][[2]])
##     rept <- row.names(summary(object[[3]])[[1]])
##     rept_lab <- trim.blanks(rept[-grep("Residuals", rept)])
##     rept <- trim.blanks(rept[-grep(":|Residuals", rept)])
##     rept <- sub("C\\(([A-Za-z_\\.]*).*\\)", "\\1", rept)
##     data <- eval(call[[3]])
##     between <- row.names(summary(object[[2]])[[1]])
##     between_lab <- trim.blanks(between[-grep("Residuals", between)])
##     between <- trim.blanks(between[-grep(":|Residuals", between)])
##     data.rm <- reshape(data, direction = "wide", timevar = rept, 
##         idvar = block, v.names = resp)
##     if (length(between_lab) == 0) {
##         data.manova <- lm(as.matrix(data.rm[grep(resp, colnames(data.rm))]) ~ 
##             1)
##     }
##     else {
##         between <- eval(parse(text = paste("data.rm$", between, 
##             sep = "")))
##         data.manova <- lm(as.matrix(data.rm[grep(resp, colnames(data.rm))]) ~ 
##             between)
##     }
##     idata <- data.frame(within = as.factor(levels(eval(parse(text = paste(as.character(call[[3]]), 
##         "$", rept, sep = ""))))))
##     data.ok <- Anova(data.manova, idata = idata, idesign = ~within, 
##         type = "III")
##     gg.hf <- GGHF(data.ok)
##     colnames(gg.hf) <- c("Greenhouse-Geisser", "Hunyh-Feldt")
##     rownames(gg.hf) <- c(between_lab, rept_lab)
##     cat("Epsilon values\n")
##     gg.hf
## }


AnovaM <- function (mod, type = "I", denoms = "Resid", split = list(), 
    ...) 
{
    dots <- list(...)
    if (type == "IIII") {
        if (suppressWarnings(!mod == "aovlist" && is.null(split))) {
            mod <- aov(mod)
            mod <- do.call("anova", c(list(mod, dots)))
            Anovam(mod, denoms = denoms)
        }
        else {
            mod <- aov(mod)
            mod <- do.call("summary", c(list(mod, intercept = FALSE, 
                split = split, expand.split = TRUE, keep.zero.df = TRUE, 
                dots)))
            Anovam(mod, denoms = denoms)
        }
    }
    else UseMethod("AnovaM", mod)
}



## AnovaM.aov <- function (object, type = "I", error, intercept = FALSE, split = NULL, 
##     expand.split = TRUE, keep.zero.df = TRUE, denoms = "Resid", 
##     RM = F, ...) 
## {
##     require(car)
##     dots <- list(...)
##     op <- options()
##     options(warn = -1)
##     if (!is.null(split)) {
##         s <- do.call("summary", c(list(object, intercept = intercept, 
##             split = split, expand.split = expand.split, keep.zero.df = keep.zero.df, 
##             dots)))
##     }
##     else {
##         s <- do.call("summary", c(list(object, intercept = intercept, 
##             expand.split = expand.split, keep.zero.df = keep.zero.df)))
##     }
##     names <- character(0)
##     if (length(s) > 1) {
##         for (i in 1:length(s)) {
##             nms <- attr(s[[i]][[1]], "row.names")
##             nms <- sub("^ +", "", nms)
##             nms <- sub(" +$", "", nms)
##             nms <- nms[-length(nms)]
##             names <- c(names, nms)
##         }
##     }
##     else {
##         if (length(s[[1]]) == 1) 
##             nms <- attr(s[[1]][[1]], "row.names")
##         else nms <- attr(s[[1]], "row.names")
##         nms <- sub("^ +", "", nms)
##         nms <- sub(" +$", "", nms)
##         nms <- nms[-length(nms)]
##         names <- c(names, nms)
##         if (!missing(error)) {
##             sumry <- summary(error, corr = FALSE)
##             s2 <- sumry$sigma^2
##             error.df <- error$df.residual
##             nms <- row.names(s[[1]])
##             nms <- sub("^ +", "", nms)
##             nms <- sub(" +$", "", nms)
##             ind <- match("Residuals", nms, nomatch = 0)
##             s[[1]][ind, "Df"] <- error.df
##             s[[1]][ind, "Sum Sq"] <- s2 * error.df
##             s[[1]][ind, "Mean Sq"] <- s2
##         }
##     }
##     if (length(s) > 1) 
##         form <- attr(object, "call")
##     else {
##         form <- object$call
##         if (is.null(form)) 
##             form <- attr(object, "call")
##     }
##     original_form <- form
##     form <- deparse(form)
##     form <- sub("Error\\(([A-Za-z0-9_\\.\\/\\,\\ \\(\\)]*)\\)", 
##         "\\1", form)
##     if (RM == T) {
##         gg.hf <- epsi.GG.HF(object)
##         print(gg.hf)
##     }
##     form <- eval(parse(text = form))
##     x <- model.matrix(form)
##     offset <- model.offset(model.frame(form))
##     asgn <- attr(x, "assign")
##     tl <- attr(form$terms, "term.labels")
##     fac <- attr(form$terms, "factors")
##     split.fact <- vector("list", length(tl))
##     names(split.fact) <- tl
##     for (k in 1:length(tl)) {
##         scope <- attr(terms(update.formula(form, paste("~", tl[k], 
##             sep = ""))), "term.labels")
##         ndrop <- match(scope, tl)
##         ii <- seq(along = asgn)[asgn == ndrop]
##         split.fact[k] <- vector("list", length(ii))
##         split.fact[k] <- list(ii)
##     }
##     names.scope <- vector("list", length(names))
##     names(names.scope) <- names
##     main.scope <- names.scope
##     for (ll in 1:length(names)) {
##         for (mm in 1:length(split.fact)) {
##             if (names(names.scope[ll]) == names(split.fact[mm])) {
##                 names.scope[[ll]] <- split.fact[[mm]]
##                 main.scope[[ll]] <- names(split.fact[mm])
##             }
##         }
##     }

##     if (!is.null(split)) {
##         split.cont <- split
##         split.cont <- splitInteractions(split.cont, factors = attr(form$terms, 
##             "factors"), c("(Intercept)", tl), asgn, names(form$coef))
##         for (ll in 1:length(names)) {
##             for (mm in 1:length(split.cont)) {
##                 for (nn in 1:length(split.cont[[mm]])) {
##                   if (names(names.scope[ll]) == paste(names(split.cont[mm]), 
##                     ": ", names(split.cont[[mm]][nn]), sep = "")) {
##                     names.scope[[ll]] <- split.fact[[names(split.cont[mm])]][split.cont[[mm]][[nn]]]
##                     main.scope[[ll]] <- names(split.cont[mm])
##                   }
##                 }
##             }
##         }
##     }
##     chisq <- deviance(form)
##     y <- form$residuals + predict(form)
##     tol <- .Machine$double.eps^0.5
##     for (ll in 1:length(names.scope)) {
##         if (type == "II") {
##             if (ncol(fac) > 1) {
##                 ndrop <- relatives(main.scope[[ll]], main.scope, 
##                   fac)
##                 ii <- seq(along = asgn)[asgn == ndrop]
##                 jj <- setdiff(seq(ncol(x)), ii)
##             }
##             else jj <- seq(ncol(x))
##             z <- lm.fit(x[, jj, drop = FALSE], y, offset = offset)
##             oldClass(z) <- "lm"
##             chisq <- deviance(z)
##             jj <- setdiff(jj, names.scope[[ll]])
##             z <- lm.fit(x[, jj, drop = FALSE], y, offset = offset)
##             oldClass(z) <- "lm"
##             RSS <- deviance(z)
##         }
##         if (type == "III") {
##             jj <- setdiff(seq(ncol(x)), names.scope[[ll]])
##             z <- lm.fit(x[, jj, drop = FALSE], y, offset = offset)
##             oldClass(z) <- "lm"
##             RSS <- deviance(z)
##         }
##         if (type != "I") 
##             SS <- RSS - chisq
##         if (length(s) > 1 & type != "I") {
##             for (mm in 1:length(s)) {
##                 nms <- row.names(s[[mm]][[1]])
##                 nms <- sub("^ +", "", nms)
##                 nms <- sub(" +$", "", nms)
##                 ind <- match(names(names.scope[ll]), nms, nomatch = 0)
##                 if (ind && SS > tol) {
##                   s[[mm]][[1]][ind, "Sum Sq"] <- SS
##                   s[[mm]][[1]][ind, "Mean Sq"] <- SS/s[[mm]][[1]][ind, 
##                     "Df"]
##                   s[[mm]][[1]][ind, "F value"] <- s[[mm]][[1]][ind, 
##                     "Mean Sq"]/s[[mm]][[1]]["Residuals", "Mean Sq"]
##                   s[[mm]][[1]][ind, "Pr(>F)"] <- pf(s[[mm]][[1]][ind, 
##                     "F value"], s[[mm]][[1]][ind, "Df"], s[[mm]][[1]]["Residuals", 
##                     "Df"], lower.tail = FALSE)
##                   if (RM == T && 1 == 2) {
##                     if (!is.na(match(trim.blanks(rownames(s[[mm]][[1]][ind, 
##                       ])), rownames(gg.hf)))) 
##                       mtch <- match(trim.blanks(rownames(s[[mm]][[1]][ind, 
##                         ])), rownames(gg.hf))
##                     if (is.na(mtch)) 
##                       mtch <- 1
##                     if (RM == T && gg.hf[mm][1] != "NA") {
##                       eps.GG <- gg.hf[mtch][1]
##                       dfsGR_df <- s[[mm]][[1]]["Residuals", "Df"]
##                       dfsG <- s[[mm]][[1]][ind, "Df"] * eps.GG
##                       dfsH <- s[[mm]][[1]][ind, "Df"] * gg.hf[mtch, 
##                         2]
##                       s[[mm]][[1]][ind, "GG.P"] <- 1 - pf(s[[mm]][[1]][ind, 
##                         "F value"], dfsG, dfsGR_df)
##                       s[[mm]][[1]][ind, "HF.P"] <- 1 - pf(s[[mm]][[1]][ind, 
##                         "F value"], dfsH, dfsGR_df)
##                     }
##                   }
##                 }
##                 else next
##             }
##         }
##         else {
##             if (type != "I") {
##                 nms <- row.names(s[[1]])
##                 nms <- sub("^ +", "", nms)
##                 nms <- sub(" +$", "", nms)
##                 ind <- match(names(names.scope[ll]), nms, nomatch = 0)
##                 if (ind && SS > tol) {
##                   s[[1]][ind, "Sum Sq"] <- SS
##                   s[[1]][ind, "Mean Sq"] <- SS/s[[1]][ind, "Df"]
##                   s[[1]][ind, "F value"] <- s[[1]][ind, "Mean Sq"]/s[[1]]["Residuals", 
##                     "Mean Sq"]
##                   s[[1]][ind, "Pr(>F)"] <- pf(s[[1]][ind, "F value"], 
##                     s[[1]][ind, "Df"], s[[1]]["Residuals", "Df"], 
##                     lower.tail = FALSE)
##                 }
##             }
##         }
##     }
##     if (length(s) > 1) {
##         ss <- NULL
##         for (i in 1:length(s)) {
##             if (denoms != "Resid") 
##                 s[[i]] <- Anovam(s[[i]], denoms = denoms[[i]])
##         }
##     }
##     else if (denoms != "Resid") 
##         s <- Anovam(s, denoms = denoms)
##     if (RM == T) {
##         if (length(s) > 1) {
##             for (ll in 1:length(names.scope)) {
##                 for (mm in 1:length(s)) {
##                   nms <- row.names(s[[mm]][[1]])
##                   nms <- sub("^ +", "", nms)
##                   nms <- sub(" +$", "", nms)
##                   ind <- match(names(names.scope[ll]), nms, nomatch = 0)
##                   if (!ind) 
##                     next
##                   if (!is.na(match(trim.blanks(rownames(s[[mm]][[1]][ind, 
##                     ])), rownames(gg.hf)))) 
##                     mtch <- match(trim.blanks(rownames(s[[mm]][[1]][ind, 
##                       ])), rownames(gg.hf))
##                   if (is.na(mtch)) 
##                     mtch <- 1
##                   if (RM == T && !is.na(gg.hf[mm][1])) {
##                     eps.GG <- gg.hf[mtch][1]
##                     dfsGR_df <- s[[mm]][[1]]["Residuals", "Df"]
##                     dfsG <- s[[mm]][[1]][ind, "Df"] * eps.GG
##                     dfsH <- s[[mm]][[1]][ind, "Df"] * gg.hf[mtch, 
##                       2]
##                     s[[mm]][[1]][ind, "GG.P"] <- 1 - pf(s[[mm]][[1]][ind, 
##                       "F value"], dfsG, dfsGR_df)
##                     s[[mm]][[1]][ind, "HF.P"] <- 1 - pf(s[[mm]][[1]][ind, 
##                       "F value"], dfsH, dfsGR_df)
##                   }
##                 }
##             }
##         }
##     }
##     options(op)
##     attr(s, "heading") <- c("Anova Table (Type III tests)\n", 
##         paste("Response:", responseName(object)))
##     attr(s, "call") <- original_form
##     s
##     s
## }

AnovaM.aov<-function (object, type = "I", error, intercept = FALSE, split = NULL, 
                      expand.split = TRUE, keep.zero.df = TRUE, denoms = "Resid", 
                      RM = F, epsi=NULL,...) 
{
  require(car)
  dots <- list(...)
  op <- options()
  options(warn = -1)
  if (!is.null(split)) {
    s <- do.call("summary", c(list(object, intercept = intercept, 
                                   split = split, expand.split = expand.split, keep.zero.df = keep.zero.df, 
                                   dots)))
  }
  else {
    s <- do.call("summary", c(list(object, intercept = intercept, 
                                   expand.split = expand.split, keep.zero.df = keep.zero.df)))
  }
  names <- character(0)
  if (length(s) > 1) {
    for (i in 1:length(s)) {
      nms <- attr(s[[i]][[1]], "row.names")
      nms <- sub("^ +", "", nms)
      nms <- sub(" +$", "", nms)
      nms <- nms[-length(nms)]
      names <- c(names, nms)
    }
  }
  else {
    if (length(s[[1]]) == 1) 
      nms <- attr(s[[1]][[1]], "row.names")
    else nms <- attr(s[[1]], "row.names")
    nms <- sub("^ +", "", nms)
    nms <- sub(" +$", "", nms)
    nms <- nms[-length(nms)]
        names <- c(names, nms)
    if (!missing(error)) {
      sumry <- summary(error, corr = FALSE)
      s2 <- sumry$sigma^2
      error.df <- error$df.residual
      nms <- row.names(s[[1]])
      nms <- sub("^ +", "", nms)
      nms <- sub(" +$", "", nms)
      ind <- match("Residuals", nms, nomatch = 0)
      s[[1]][ind, "Df"] <- error.df
      s[[1]][ind, "Sum Sq"] <- s2 * error.df
      s[[1]][ind, "Mean Sq"] <- s2
    }
  }
  if (length(s) > 1) 
    form <- attr(object, "call")
  else {
    form <- object$call
    if (is.null(form)) 
      form <- attr(object, "call")
  }
  original_form <- form
  form <- deparse(form,width.cutoff = 500, 
                  backtick = TRUE)
                                        #form <- sub("Error\\(([A-Za-z0-9_\\.\\/\\,\\ \\(\\)\\+]*)\\)", "\\1", form)
  form <- sub("Error\\(([A-Za-z0-9_\\.\\,\\$]*)[A-Za-z0-9_\\,\\.()\\$\\/\\ \\+\\*\\-]*\\)", "\\1", form)
  if (RM == T) {
    if (is.null(epsi)) {
      gg.hf <- epsi.GG.HF(object)
      cat("\n   Sphericity Epsilon Values \n-------------------------------\n")
      spher.table<-data.frame("Greenhouse-Geisser"=gg.hf[[1]], "Huynh-Feldt"=gg.hf[[2]])
      rownames(spher.table) <- ""
      print(spher.table)
      
    }
    else {
      gg.hf<-epsi
    }
  }
  form <- eval(parse(text = form))
  x <- model.matrix(form)
  offset <- model.offset(model.frame(form))
  asgn <- attr(x, "assign")
  tl <- attr(form$terms, "term.labels")
  fac <- attr(form$terms, "factors")
  split.fact <- vector("list", length(tl))
  names(split.fact) <- tl
  for (k in 1:length(tl)) {
    scope <- attr(terms(update.formula(form, paste("~", tl[k], 
                                                   sep = ""))), "term.labels")
    ndrop <- match(scope, tl)
    ii <- seq(along = asgn)[asgn == ndrop]
    split.fact[k] <- vector("list", length(ii))
    split.fact[k] <- list(ii)
  }
  names.scope <- vector("list", length(names))
  names(names.scope) <- names
  main.scope <- names.scope
  for (ll in 1:length(names)) {
    for (mm in 1:length(split.fact)) {
      if (names(names.scope[ll]) == names(split.fact[mm])) {
        names.scope[[ll]] <- split.fact[[mm]]
        main.scope[[ll]] <- names(split.fact[mm])
      }
    }
  }
  if (!is.null(split)) {
    split.cont <- split
    split.cont <- splitInteractions(split.cont, factors = attr(form$terms, 
                                                  "factors"), c("(Intercept)", tl), asgn, names(form$coef))
    for (ll in 1:length(names)) {
      for (mm in 1:length(split.cont)) {
        for (nn in 1:length(split.cont[[mm]])) {
          if (names(names.scope[ll]) == paste(names(split.cont[mm]), 
                     ": ", names(split.cont[[mm]][nn]), sep = "")) {
            names.scope[[ll]] <- split.fact[[names(split.cont[mm])]][split.cont[[mm]][[nn]]]
            main.scope[[ll]] <- names(split.cont[mm])
          }
        }
      }
    }
  }
  chisq <- deviance(form)
  y <- form$residuals + predict(form)
  tol <- .Machine$double.eps^0.5
  for (ll in 1:length(names.scope)) {
    if (type == "II") {
      if (ncol(fac) > 1) {
        ndrop <- relatives(main.scope[[ll]], main.scope, 
                           fac)
        ii <- seq(along = asgn)[asgn == ndrop]
        jj <- setdiff(seq(ncol(x)), ii)
      }
      else jj <- seq(ncol(x))
      z <- lm.fit(x[, jj, drop = FALSE], y, offset = offset)
      oldClass(z) <- "lm"
      chisq <- deviance(z)
      jj <- setdiff(jj, names.scope[[ll]])
      z <- lm.fit(x[, jj, drop = FALSE], y, offset = offset)
      oldClass(z) <- "lm"
      RSS <- deviance(z)
    }
    if (type == "III") {
      jj <- setdiff(seq(ncol(x)), names.scope[[ll]])
      z <- lm.fit(x[, jj, drop = FALSE], y, offset = offset)
      oldClass(z) <- "lm"
      RSS <- deviance(z)
    }
    if (type != "I") 
      SS <- RSS - chisq
    if (length(s) > 1 & type != "I") {
      for (mm in 1:length(s)) {
        nms <- row.names(s[[mm]][[1]])
        nms <- sub("^ +", "", nms)
        nms <- sub(" +$", "", nms)
        ind <- match(names(names.scope[ll]), nms, nomatch = 0)
        if (ind && SS > tol) {
          s[[mm]][[1]][ind, "Sum Sq"] <- SS
          s[[mm]][[1]][ind, "Mean Sq"] <- SS/s[[mm]][[1]][ind, 
                                                          "Df"]
          s[[mm]][[1]][ind, "F value"] <- s[[mm]][[1]][ind, 
                                                       "Mean Sq"]/s[[mm]][[1]]["Residuals", "Mean Sq"]
          s[[mm]][[1]][ind, "Pr(>F)"] <- pf(s[[mm]][[1]][ind, 
                                                         "F value"], s[[mm]][[1]][ind, "Df"], s[[mm]][[1]]["Residuals", 
                                                                                                           "Df"], lower.tail = FALSE)
          if (RM == T && 1 == 2) {
            if (!is.na(match(trim.blanks(rownames(s[[mm]][[1]][ind, 
                                                               ])), rownames(gg.hf)))) 
              mtch <- match(trim.blanks(rownames(s[[mm]][[1]][ind, 
                                                              ])), rownames(gg.hf))
            if (is.na(mtch)) 
              mtch <- 1
            if (RM == T && gg.hf[mm][1] != "NA") {
              eps.GG <- gg.hf[mtch][1]
              dfsGR_df <- s[[mm]][[1]]["Residuals", "Df"]
              dfsG <- s[[mm]][[1]][ind, "Df"] * eps.GG
              dfsH <- s[[mm]][[1]][ind, "Df"] * gg.hf[mtch, 
                                                      2]
              s[[mm]][[1]][ind, "GG.P"] <- 1 - pf(s[[mm]][[1]][ind, 
                                                               "F value"], dfsG, dfsGR_df)
              s[[mm]][[1]][ind, "HF.P"] <- 1 - pf(s[[mm]][[1]][ind, 
                                                               "F value"], dfsH, dfsGR_df)
            }
          }
        }
        else next
      }
    }
    else {
      if (type != "I") {
        nms <- row.names(s[[1]])
        nms <- sub("^ +", "", nms)
        nms <- sub(" +$", "", nms)
        ind <- match(names(names.scope[ll]), nms, nomatch = 0)
        if (ind && SS > tol) {
          s[[1]][ind, "Sum Sq"] <- SS
          s[[1]][ind, "Mean Sq"] <- SS/s[[1]][ind, "Df"]
          s[[1]][ind, "F value"] <- s[[1]][ind, "Mean Sq"]/s[[1]]["Residuals", 
                                                                  "Mean Sq"]
          s[[1]][ind, "Pr(>F)"] <- pf(s[[1]][ind, "F value"], 
                                      s[[1]][ind, "Df"], s[[1]]["Residuals", "Df"], 
                                      lower.tail = FALSE)
        }
      }
    }
  }
  if (length(s) > 1) {
    ss <- NULL
    for (i in 1:length(s)) {
      if (denoms != "Resid")
        s[[i]] <- Anovam(s[[i]], denoms = denoms[[i]])
    }
  }
  else if (denoms != "Resid")
    s <- Anovam(s, denoms = denoms)
  sHF<- sGG <- s
  #Made a change to remove responseName() function.  It is now an internal of car that is not accessable!
  #alternatively could have changed to
  # car:::responseNames.default()
  respName <- deparse(attr(terms(object), "variables")[[2]])
  attr(s, "heading") <- c(paste("Anova Table (Type",type, "tests)\n"), 
                          paste("Response:", respName))
  
  attr(s, "call") <- original_form
  
  if (RM == T) {
    if (length(s) > 1) {
      for (ll in 1:length(names.scope)) {
        for (mm in 1:length(s)) {
          nms <- row.names(s[[mm]][[1]])
          nms <- sub("^ +", "", nms)
          nms <- sub(" +$", "", nms)
          ind <- match(names(names.scope[ll]), nms, nomatch = 0)
          if (!ind) 
            next
          #gh<-gg.hf[match(names(names.scope[ll]),row.names(gg.hf)),]
          if (RM == T) {# && !is.na(gh[1])) {
            eps.GG <- gg.hf[[1]] #gg.hf[mtch][1]
            dfsGR_df <- s[[mm]][[1]]["Residuals", "Df"]
            dfsG <- s[[mm]][[1]][ind, "Df"] * eps.GG
            dfsH <- s[[mm]][[1]][ind, "Df"] * gg.hf[[2]]#gg.hf[mtch,2]
            sGG[[mm]][[1]][ind, "Pr(>F)"] <- 1 - pf(s[[mm]][[1]][ind, 
                                                                 "F value"], dfsG, dfsGR_df)
            sGG[[mm]][[1]][ind,"Df"] <-dfsG
            sHF[[mm]][[1]][ind, "Pr(>F)"] <- 1 - pf(s[[mm]][[1]][ind, 
                                                                 "F value"], dfsH, dfsGR_df)
            sHF[[mm]][[1]][ind,"Df"] <-dfsH
          }
        }
      }
    }
    attr(sGG, "heading") <- c("Greenhouse-Geisser corrected ANOVA table\n", 
                                      paste("Response:", respName))
    attr(sHF, "heading") <- c("Huynh-Feldt corrected ANOVA table\n", 
                                      paste("Response:", respName))
    s<-list(s,"Greenhouse-Geisser"=sGG,sHF)
    attr(s, "class") <- c("RM",class(s))
  }
  options(op)

  #s
  
  s
}



print.RM<-function(ob)
  {
    for (i in 1:length(ob)){
      cat("\n\n",attr(ob[[i]],"heading"),sep="")
     print(ob[[i]]) 
    }
  }

AnovaM.aovlist <- function (object, ...) 
{
    AnovaM.aov(object, ...)
}

AnovaM.lm <- function (object, ...) 
{
    AnovaM.aov(aov(object), ...)
}


Anovam <- function (object, denoms = "Resid") 
{
    if (length(object[[1]]) > 1) 
        a <- object[[1]]
    else {
        a <- object[[1]][[1]]
        if (is.null(a)) 
            a <- object[[1]]
    }
    var_names <- rownames(a)
    terms <- grep("^[^ ]", rownames(a))
    var_length <- length(var_names) - 1
    if (match("summary.aov", class(denoms), nomatch = 0)) {
        a[length(rownames(a)), "Sum Sq"] <- denoms[[1]]["Residuals", 
            "Sum Sq"]
        a[length(rownames(a)), "Df"] <- denoms[[1]]["Residuals", 
            "Df"]
        a[length(rownames(a)), "Mean Sq"] <- denoms[[1]]["Residuals", 
            "Sum Sq"]/denoms[[1]]["Residuals", "Df"]
        b <- denoms[[1]]["Resid", ]
    }
    if (match("list", class(denoms), nomatch = 0)) {
        b <- NULL
        for (i in 1:length(denoms)) b <- rbind(b, denoms[[i]])
        a[length(rownames(a)), "Sum Sq"] <- b[length(denoms), 
            "Sum Sq"]
        a[length(rownames(a)), "Df"] <- b[length(denoms), "Df"]
        a[length(rownames(a)), "Mean Sq"] <- b[length(denoms), 
            "Sum Sq"]/b[length(denoms), "Df"]
    }
    if (match("anova", class(denoms), nomatch = 0)) {
        if (dim(denoms)[1] != 1) 
            denoms <- denoms[rownames(denoms)[length(rownames(denoms))], 
                ]
        b <- denoms
        bb <- b
        if (dim(b)[1] == 1) 
            for (i in 1:var_length) bb <- rbind(bb, b)
        b <- bb
        a[length(rownames(a)), "Sum Sq"] <- b[var_length, "Sum Sq"]
        a[length(rownames(a)), "Df"] <- b[var_length, "Df"]
        a[length(rownames(a)), "Mean Sq"] <- b[var_length, "Sum Sq"]/b[var_length, 
            "Df"]
    }
    if (match("character", class(denoms), nomatch = 0)) {
        b <- a
        if (length(denoms) == 1 | length(denoms) == length(terms) - 
            1 | length(denoms) == var_length) {
        }
        else stop("The number of denominators is not a function of the number of tests")
        if (length(denoms) == 1) 
            denoms <- rep(denoms, var_length)
        if (length(denoms) > length(terms)) 
            terms <- grep("", rownames(a))
        b <- a[denoms, ]
    }
    for (i in 1:var_length) {
        if (!is.na(match(i, terms))) 
            j <- match(i, terms)
        SS <- a[i, "Sum Sq"]
        MS <- SS/a[i, "Df"]
        a[i, "F value"] <- f <- MS/b[j, "Mean Sq"]
        a[i, "Pr(>F)"] <- 1 - pf(f, a[i, "Df"], b[j, "Df"])
    }
    a
}




Hist <- function (x, scale = c("frequency", "percent", "density"), ...) 
{
    xlab <- deparse(substitute(x))
    xlab <- (strsplit(xlab, "\\$")[[1]][2])
    x <- na.omit(x)
    scale <- match.arg(scale)
    if (scale == "frequency") 
        hist(x, xlab = xlab, main = "", ...)
    else if (scale == "density") 
        hist(x, freq = FALSE, xlab = xlab, main = "", ...)
    else {
        n <- length(x)
        hist(x, axes = FALSE, xlab = xlab, ylab = "Percent", 
            main = "", ...)
        axis(1)
        max <- ceiling(10 * par("usr")[4]/n)
        at <- if (max <= 3) 
            (0:(2 * max))/20
        else (0:max)/10
        axis(2, at = at * n, labels = at * 100)
    }
    box()
    abline(h = 0, col = "gray")
    invisible(NULL)
}




Mbargraph <- function (dv, factorlist, ylabel = paste("Mean", deparse(substitute(dv))), 
    xlabel = "", type = "b", ymin = NULL, ymax = NULL, errorbars = "se", 
    symbols = "") 
{
     means <- tapply(dv, factorlist, mean, na.rm = T)
    sds <- tapply(dv, factorlist, sd, na.rm = T)
    ns <- tapply(dv, factorlist, function(x) length(na.omit(x)))
    ses <- sds/sqrt(ns)
    ys <- pretty(c(means - ses, means + (2 * ses)))
    if (!is.null(ymin)) 
        ys <- c(ys, ymin)
    if (!is.null(ymax)) 
        ys <- c(ys, ymax)
    if (type == "b") {
        if (is.list(factorlist)) {
            xs <- barplot(means, ylim = c(min(ys), max(ys)), 
                beside = T, axes = F, xpd = F, axisnames = F, ylabel="", xlabel="",
                axis.lty = 2, legend.text = T, args.legend=list(title=gsub(".*\\$(.*)\\(\\)","\\1",deparse(substitute(factorlist)[2]))))
            #print(deparse(substitute(factorlist)))
            if (xlabel == "") {
                xlabel = deparse(substitute(factorlist)[[3]])
            }
            xAT<-apply(xs,2,mean)
            xLab<-colnames(means)
        }
        else {
            xs <- barplot(means, ylim = c(min(ys), max(ys)), 
                beside = T, axes = F, xpd = F, axis.lty = 1, 
                legend.text = F,ann=F)
            xAT<-xs[,1]
            xLab<-names(means)
            #if (xlabel == "") {
            #    xlabel = deparse(substitute(factorlist))
            #}
        }
    }
    else {
        if (is.list(factorlist)) {
            xs <- barplot(means, ylim = c(min(ys), max(ys)), 
                beside = T, axes = F, xpd = F, axis.lty = 1, 
                legend.text = F, plot = FALSE)
            interaction.plot(factorlist[[1]], factorlist[[2]], 
                dv, ylim = c(min(ys), max(ys)))
            xAT<-apply(xs,2,mean)
            xLab<-colnames(means)
        }
        else {
            xs <- plotMeans(means, factorlist[[1]], factorlist[[2]],axes = F)
            xAT=xs[,1]
            xLab<-names(means)
        }
    }
    #library(Hmisc)
    #errbar(xs, means, means + ses, means - ses, add = T)
     print(xAT)
     print(xLab)
    axis(1,at=xAT, lab=xLab)
    arrows(xs, means+ses, xs, means-ses, ang=90, length=0.1, code=3)
    axis(2, las = 1)
    mtext(2, text = ylabel, line = 3, cex = 1.5)
    mtext(1, text = xlabel, line = 3, cex = 1.5)
    lines(c(0, max(xs) + 0.5), c(min(ys), min(ys)))
    Mgraph <- data.frame(xs, means = means + ses + ((max(ys) - 
        min(ys))/20))
    if (exists("symbols")) {
        text(xs, means + ses, labels = symbols, pos = 3)
    }
}

Mbargraph1<-function (x.factor, trace.factor = NULL, response, legend = T, 
    names.arg = NULL, xlabel = deparse(substitute(x.factor)), 
    axes = T, xlabs = levels(x.factor), ylabel = paste("Mean", 
        deparse(substitute(response))), ylim = NULL, errorbars = "se", 
    col = c(0, 1)) 
{
    if (is.null(trace.factor)) {
        cells <- tapply(response, x.factor, mean, na.rm = T)
        sds <- tapply(response, x.factor, sd, na.rm = T)
        ns <- tapply(response, x.factor, function(x) length(na.omit(x)))
        ses <- sds/sqrt(ns)
        col <- col[2]
    }
    else {
        cells <- t(tapply(response, list(x.factor, trace.factor), 
            mean, na.rm = T))
        sds <- t(tapply(response, list(x.factor, trace.factor), 
            sd, na.rm = T))
        ns <- t(tapply(response, list(x.factor, trace.factor), 
            function(x) length(na.omit(x))))
        ses <- sds/sqrt(ns)
        leg.text <- levels(trace.factor)
    }
    ys <- pretty(c(cells - ses, cells + (2 * ses)))
    if (!is.null(ylim)) 
        ys <- pretty(ylim)
    if (!is.null(trace.factor)) {
        xs <- barplot(cells, ylim = c(min(ys), max(ys)), beside = T, 
            axes = F, xpd = F, axisnames = F, axis.lty = 2, legend.text = F, 
            col = col)
    }
    else {
        xs <- barplot(cells, ylim = c(min(ys), max(ys)), beside = T, 
            axes = F, xpd = F, axis.lty = 1, legend.text = F, 
            space = 1, col = col)
    }
    if (axes == T) {
        if (is.null(trace.factor)) 
            xm <- xs
        else xm <- apply(xs, 2, median)
        axis(2, las = 1)
        axis(1, at = xm, lab = xlabs, padj = 1, mgp = c(0, 0, 
            0))
        mtext(2, text = ylabel, line = 3, cex = 1.5)
        mtext(1, text = xlabel, line = 3, cex = 1.5)
        box(bty = "l")
    }
    arrows(xs, cells, xs, cells + ses, code = 3, angle = 90, 
        len = 0.05)
    if (legend != F & !is.null(trace.factor)) 
        legend("topright", leg = legend, fill = col, col = col, 
            bty = "n", cex = 1.5)
}



## Model.selection <- function (object, data) 
## {
##     library(hier.part)
##     library(nlme)
##     library(R.utils)
##     predictors <- dimnames(attr(object$terms, "factors"))[[2]]
##     comb <- combos(length(predictors))[[1]]
##     ll <- length(comb[, 1])
##     var.names <- vector(mode = "character", length = ll)
##     d <- data.frame(r.squared = 1:ll, Adj.r.sq = 1:ll, P = 1:ll, 
##         Cp = 1:ll, AIC = 1:ll, BIC = 1:ll, Select = factor(rep(" ", 
##             ll), levels = c(" ", "*")))
##     form <- object$call$formula
##     minAIC = 1e+05
##     minBIC = 1e+05
##     minCp <- 1e+05
##     for (i in 1:length(comb[, 1])) {
##         pred <- predictors[comb[i, ]]
##         var.names[i] <- paste(pred, collapse = "+")
##         ii <- paste(pred, collapse = "+", sep = "")
##         new_fit <- update.formula(form, paste("~  +", ii))
##         fit <- lm(new_fit, data = data)
##         fit_sum <- summary(lm(fit))
##         d$r.squared[i] <- fit_sum$r.squared
##         d$Adj.r.sq[i] <- fit_sum$adj.r.squared
##         d$P[i] <- fit_sum$fstatistic[["numdf"]]
##         d$Cp[i] <- Cp(fit, object)
##         deltaCp <- abs(d$Cp[i] - d$P[i])
##         d$AIC[i] <- extractAIC(fit)[2]
##         if (minCp > d$Cp[i]) 
##             minCp <- d$Cp[i]
##         if (minAIC > d$AIC[i]) 
##             minAIC <- d$AIC[i]
##         d$BIC[i] <- extractAIC(fit, k = log(nrow(data)))[2]
##         if (minBIC > d$BIC[i]) 
##             minBIC <- d$BIC[i]
##     }
##     rownames(d) <- var.names
##     d$Select[d$BIC == minBIC] <- "*"
##     d$Select[d$AIC == minAIC] <- "*"
##     d$Select[d$Cp == minCp] <- "*"
##     d
## }

#qAIC is used as an alternative to AIC when overdispersal
qAIC<-function(object, based="deviance"){
  qaic <- vector(mode="numeric", length=2)
  k<-attr(logLik(object),"df")
  chat <-object$deviance/object$df.residual
  qaic[1]<- -(logLik(object)[1]/chat)+(2*k)
  chat <-sum(resid(object, type="pearson")^2)/object$df.residual
  qaic[2]<- -(logLik(object)[1]/chat)+(2*k)
  names(qaic) <- c("Deviance based", "Residuals based")
  ifelse(based=="residuals", qaic[2],qaic[1])
}

qAICc<-function(object,based="deviance"){
  qaic <- vector(mode="numeric", length=2)
  k<-attr(logLik(object),"df")
  chat <-object$deviance/object$df.residual
  qaic[1]<- -(logLik(object)[1]/chat)+(2*k)+(((2*k)*(k+1))/(length(object$data[,1])-k-1))
  chat <-sum(resid(object, type="pearson")^2)/object$df.residual
  -(logLik(object)[1]/chat)+(2*k)+(((2*k)*(k+1))/(length(object$data[,1])-k-1))
  qaic[2]<- -(logLik(object)[1]/chat)+(2*k)+(((2*k)*(k+1))/(length(object$data[,1])-k-1))
  names(qaic) <- c("Deviance based", "Residuals bases")
  ifelse(based=="residuals", qaic[2],qaic[1])
}



AICc<-function(object, aic){
	k<-attr(logLik(object),"df")
	n<-attr(logLik(object),"nobs")
    #Note, the above line was different in the AICc function definition in glm.R (Book)
	AICc<-aic+(2*k*(k+1))/(n-k-1)
	AICc	
}

AIC.glm<-function(object){
  (object$deviance)+(2*object$rank)
}

Cp <- function(model, full){
  SSres.red<-anova(model)["Residuals","Sum Sq"]
  MSres.full<-anova(full)["Residuals","Mean Sq"]
  n<-nrow(model$model)
  p<-summary(model)$fstatistic[["numdf"]]
  ((SSres.red/MSres.full)-(n-2*(p+1)))
                                
}

## Model.selection <- function(object,data){
## 	library(hier.part)
## 	library(nlme)
##         #data<-object$model
## 	predictors <- dimnames(attr(object$terms, "factors"))[[2]]
## 	comb <-combos(length(predictors))[[1]]
## 	ll<-length(comb[,1])
## 	var.names <- vector(mode="character",length=ll)
## 	d <- data.frame("r.squared"=1:ll,"Adj.r.sq"=1:ll, "P"=1:ll,"Cp"=1:ll,"AIC"=1:ll, "AICc"=1:ll, deltaAIC = 1:ll, wAIC = 1:ll, "BIC"=1:ll, "Select"=factor(rep(" ",ll),levels=c(" ","*")))
## 	form <- object$call$formula
## 	minAIC=100000
## 	minBIC =100000
## 	minCp<-100000
##         fits <- list()
## 	e <- matrix(NA,nrow=ll,ncol=length(predictors), dimnames=list(1:ll, predictors))
##         ese <- matrix(NA,nrow=ll,ncol=length(predictors), dimnames=list(1:ll, predictors))
## 	for (i in 1:length(comb[,1]))
## 	{
## 		pred <- predictors[comb[i,]]
## 		var.names[i]<-paste(pred, collapse="+") 
## 		#ii <- paste("`",pred,"`", collapse="+", sep="")
## 	      ii <- paste(pred, collapse="+", sep="")
                
## 		new_fit <-update.formula(form, paste("~  +",  ii))
## 		fit<-lm(new_fit, data=data)
##                 fits[i]<-fit
##                 fit_sum <-summary(lm(fit))
##                 e[i,pred]<-fit[[1]][pred]
##                 ese[i,pred] <- fit_sum$coef[pred,2]
## 		d$"r.squared"[i] <-fit_sum$"r.squared"
## 		d$"Adj.r.sq"[i] <-fit_sum$"adj.r.squared"
##                 d$"P"[i] <- fit_sum$fstatistic[["numdf"]]
##                 d$"Cp"[i] <- Cp(fit,object)
## 			deltaCp<-abs(d$"Cp"[i]-d$"P"[i])
##             		d$"AIC"[i]<-extractAIC(fit)[2]
## 			if (minCp>d$"Cp"[i]) minCp <- d$"Cp"[i]
## 		if (minAIC>d$"AIC"[i]) minAIC <- d$"AIC"[i]
## 		AICk<-attr(logLik(object),"df")
## 		AICn<-attr(logLik(object),"nobs")
## 		d$"AICc"[i]<-d$"AIC"[i]+(2*AICk*(AICk+1))/(AICn-AICk-1)
## 		d$"BIC"[i] <-extractAIC(fit,k=log(nrow(data)))[2]
## 		if (minBIC>d$"BIC"[i]) minBIC <- d$"BIC"[i]

## 	}
##         d$deltaAIC <- d$"AICc"-min(d$"AICc")
##         d$wAIC <- round(exp(-d$deltaAIC/2)/sum(exp(-d$deltaAIC/2)),3)
##         f<-e*d$wAIC
##     g <-apply(f,2,sum, na.rm=T)
##     print("Model-averaged estimates")
##     print(g)
##     gg <- matrix(g,ncol=length(g), nrow=ll, byrow=T)
##     #print(e[e>0])
##     av_e <-apply(e,2,mean, na.rm=T)
##     unconditional_se <- apply(d$wAIC*sqrt((ese^2)+(e-gg)^2),2,sum,na.rm=T)
##     print("Unconditional SEs")
##     print(unconditional_se)
##     print("Confidence intervals")
##     print(unconditional_se*1.96)

##     model_averaging <- cbind(Estimate=g,Unconditional_SE= unconditional_se, Lower95CI=g-unconditional_se*1.96, Upper95CI=g+unconditional_se*1.96)
##         print(model_averaging)
## 	rownames(d) <- var.names
## 	d$"Select"[d$BIC==minBIC]<-"*"
## 	d$"Select"[d$AIC==minAIC]<-"*"
## 	d$"Select"[d$Cp==minCp]<-"*"
##       d
## }

Model.selection <- function(object){
	library(hier.part)
	library(nlme)
        #data<-object$model
	predictors <- dimnames(attr(object$terms, "factors"))[[2]]
	comb <-combos(length(predictors))[[1]]
	ll<-length(comb[,1])
	var.names <- vector(mode="character",length=ll)
	d <- data.frame("r.squared"=1:ll,"Adj.r.sq"=1:ll, "P"=1:ll,"Cp"=1:ll,"AIC"=1:ll, "AICc"=1:ll, deltaAIC = 1:ll, wAIC = 1:ll, "BIC"=1:ll, "Select"=factor(rep(" ",ll),levels=c(" ","*")))
	form <- object$call$formula
	minAIC=100000
	minBIC =100000
	minCp<-100000
        fits <- list()
	e <- matrix(NA,nrow=ll,ncol=length(predictors), dimnames=list(1:ll, predictors))
        ese <- matrix(NA,nrow=ll,ncol=length(predictors), dimnames=list(1:ll, predictors))
	for (i in 1:length(comb[,1]))
	{
      #print(i)
		pred <- predictors[comb[i,]]
		var.names[i]<-paste(pred, collapse="+") 
		#ii <- paste("`",pred,"`", collapse="+", sep="")
	      ii <- paste(pred, collapse="+", sep="")
                
		#new_fit <-update.formula(form, paste("~  +",  ii))
		#fit<-lm(new_fit, data=data)
        fit<-update(object,paste("~  +",  ii))
                fits[i]<-fit
                fit_sum <-summary(lm(fit))
                e[i,pred]<-fit[[1]][pred]
                ese[i,pred] <- fit_sum$coef[pred,2]
		d$"r.squared"[i] <-fit_sum$"r.squared"
		d$"Adj.r.sq"[i] <-fit_sum$"adj.r.squared"
                d$"P"[i] <- fit_sum$fstatistic[["numdf"]]
                d$"Cp"[i] <- Cp(fit,object)
			deltaCp<-abs(d$"Cp"[i]-d$"P"[i])
            		d$"AIC"[i]<-extractAIC(fit)[2]
			if (minCp>d$"Cp"[i]) minCp <- d$"Cp"[i]
		if (minAIC>d$"AIC"[i]) minAIC <- d$"AIC"[i]
		AICk<-attr(logLik(object),"df")
		AICn<-attr(logLik(object),"nobs")
		d$"AICc"[i]<-d$"AIC"[i]+(2*AICk*(AICk+1))/(AICn-AICk-1)
		d$"BIC"[i] <-extractAIC(fit,k=log(nrow(fit$model)))[2] #nrow(data)
		if (minBIC>d$"BIC"[i]) minBIC <- d$"BIC"[i]

	}
        d$deltaAIC <- d$"AICc"-min(d$"AICc")
        d$wAIC <- round(exp(-d$deltaAIC/2)/sum(exp(-d$deltaAIC/2)),3)
        f<-e*d$wAIC
    g <-apply(f,2,sum, na.rm=T)
    #print("Model-averaged estimates")
    #print(g)
    gg <- matrix(g,ncol=length(g), nrow=ll, byrow=T)
    #print(e[e>0])
    av_e <-apply(e,2,mean, na.rm=T)
    unconditional_se <- apply(d$wAIC*sqrt((ese^2)+(e-gg)^2),2,sum,na.rm=T)
    #print("Unconditional SEs")
    #print(unconditional_se)
    #print("Confidence intervals")
    #print(unconditional_se*1.96)

    model_averaging <- cbind(Estimate=g,Unconditional_SE= unconditional_se, Lower95CI=g-unconditional_se*1.96, Upper95CI=g+unconditional_se*1.96)
        #print(model_averaging)
    attr(model_averaging,"heading") <- c("Model averaging\n", 
                                      paste("Response:", car:::responseName.default(object),"\n"))
	rownames(d) <- paste(1:length(var.names),". ", var.names, sep="")
	d$"Select"[d$BIC==minBIC]<-"*"
	d$"Select"[d$AIC==minAIC]<-"*"
	d$"Select"[d$Cp==minCp]<-"*"
    attr(d,"heading") <- c("Model selection\n", 
                           paste("Response:", car:::responseName.default(object),"\n"))
    d<-list(d,model_averaging)
    attr(d, "class") <- c("modelSelection",class(d))
    d
}

print.modelSelection<-function(ob)
  {
    for (i in 1:length(ob)){
      cat("\n\n",attr(ob[[i]],"heading"),sep="")
     print(ob[[i]]) 
    }
  }

## See the changes made in Rat analysis/Model.Selection.R
## particularly for Model.selection.glm
## put these in the one below.
## Model.selection.glm <- function(object, ab=2){
##   predictors <- dimnames(attr(object$terms, "factors"))[[2]]
##   prdcts<-lapply(lapply(strsplit(predictors,":"),abbreviate, ab),paste, collapse=":")
##   comb <- combos(length(predictors))[[1]]
##   ll <- length(comb[, 1])
##   var.names <- vector(mode = "character", length = ll)
##   v.names <- vector(mode = "character", length = ll)
##   d <- data.frame(Deviance=1:ll,AIC = 1:ll, Select = factor(rep(" ", 
##             ll), levels = c(" ", "*")))
##   form <- object$call$formula
##   minAIC = 1e+05
##   minDev = 1e+05
##   for (i in 1:length(comb[, 1])) {
##     pred <- predictors[comb[i, ]]
##     prd<-prdcts[comb[i, ]]
##     var.names[i] <- paste(prd, collapse = "+")
##     ii <- paste(pred, collapse = "+", sep = "")
##     new_fit <- update(object, paste("~  +", ii))
##     d$Deviance[i] <-(new_fit$deviance)
##     #d$AIC[i] <-AIC(new_fit)
##     d$AIC[i] <-AIC.glm(new_fit) #(new_fit$deviance)+(2*new_fit$rank)
##     #d$BIC[i] <- AIC(new_fit, k = log(nrow(object$data)))
##     if (minAIC > d$AIC[i]) 
##       minAIC <- d$AIC[i]
##     #if (minBIC > d$BIC[i]) 
##     #  minBIC <- d$BIC[i]
##     if (minDev > d$Deviance[i]) 
##       minDev <- d$Deviance[i]
##   }
##   rownames(d) <- var.names
##   d$Select[d$Deviance == minDev] <- "*"
##   d$Select[d$AIC == minAIC] <- "*"
##   cat("Predictor names abbreviated\n")
##   d
## }

Model.selection.glm <- function(object, ab=2,use="AICc"){
require(hier.part)
  #predictors <- dimnames(attr(object$terms, "term.labels"))[[2]]
predictors <- attr(object$terms, "term.labels")
  prdcts<-lapply(lapply(strsplit(predictors,":"),abbreviate, ab),paste, collapse=":")
  comb <- combos(length(predictors))[[1]]
  ll <- length(comb[, 1])
  var.names <- vector(mode = "character", length = ll)
  v.names <- vector(mode = "character", length = ll)
  d <- data.frame(Deviance=1:ll,AIC = 1:ll, AICc=1:ll, deltaAIC = 1:ll, wAIC = 1:ll, qAIC=1:ll,qAICc=1:ll,Select = factor(rep(" ", 
            ll), levels = c(" ", "*")))
  form <- object$call$formula
  minAIC = 1e+05
  minDev = 1e+05
  e <- matrix(NA,nrow=ll,ncol=length(predictors), dimnames=list(1:ll, predictors))
  eW <- matrix(NA,nrow=ll,ncol=length(predictors), dimnames=list(1:ll, predictors))
  ese <- matrix(NA,nrow=ll,ncol=length(predictors), dimnames=list(1:ll, predictors))
  for (i in 1:length(comb[, 1])) {
    pred <- predictors[comb[i, ]]
    prd<-prdcts[comb[i, ]]
    var.names[i] <- paste(prd, collapse = "+")
    ii <- paste(pred, collapse = "+", sep = "")
    new_fit <- update(object, paste("~  +", ii))
    e[i,pred]<-new_fit[[1]][pred]
    eW[i,pred]<-1
    fit_sum <-summary(glm(new_fit))
    print(i)
    ese[i,pred] <- fit_sum$coef[-1,2]
    d$Deviance[i] <-(new_fit$deviance)
    #d$AIC[i] <-AIC(new_fit)
    d$AIC[i] <-AIC.glm(new_fit) #(new_fit$deviance)+(2*new_fit$rank)
    #d$BIC[i] <- AIC(new_fit, k = log(nrow(object$data)))
    if (minAIC > d$AIC[i]) 
      minAIC <- d$AIC[i]
    #if (minBIC > d$BIC[i]) 
    #  minBIC <- d$BIC[i]
    if (minDev > d$Deviance[i]) 
      minDev <- d$Deviance[i]
    #AICk<-attr(logLik(object),"df")
     AICk<-attr(logLik(new_fit),"df")
    #AICn<-ifelse(!is.null(attr(logLik(object),"nobs")),attr(logLik(object),"nobs"),length(object$data[,1]))
    AICn<-ifelse(!is.null(attr(logLik(new_fit),"nobs")),attr(logLik(new_fit),"nobs"),length(new_fit$data[,1]))
    d$"AICc"[i]<-d$"AIC"[i]+(2*AICk*(AICk+1))/(AICn-AICk-1)
    #chat <-new_fit$deviance/new_fit$df.residual
    ##d$"qAICc"[i] <-(logLik(object)[1]/chat)+(2*AICk)
    #d$"qAIC"[i] <- -(logLik(new_fit)[1]/chat)+(2*AICk)
    #d$"qAICc"[i] <- -(logLik(new_fit)[1]/chat)+(2*AICk)+(((2*AICk)*(AICk+1))/(length(new_fit$data[,1])-AICk-1))

    d$"qAIC"[i] <- qAIC(new_fit)
    d$"qAICc"[i] <-qAICc(new_fit)
  }

  #NOTE FOR TERRY 22 JUNE 2009
  #run the following line for correct dispersion
  if (use=="AICc") d$deltaAIC <- d$"AICc"-min(d$"AICc")
  else d$deltaAIC <- d$"qAICc"-min(d$"qAICc")
  #run the following if you have overdispersion
  #d$deltaAIC <- d$"qAICc"-min(d$"qAICc")
  #print("HERE")

  #not sure why I rounded it!!!
 # d$wAIC <- round(exp(-d$deltaAIC/2)/sum(exp(-d$deltaAIC/2)),3)
  d$wAIC <- exp(-d$deltaAIC/2)/sum(exp(-d$deltaAIC/2))
  wwAIC<-exp(-d$deltaAIC/2)/sum(exp(-d$deltaAIC/2))
  
  
  f<-e*d$wAIC
  g <-apply(f,2,sum, na.rm=T)
  gg <- matrix(g,ncol=length(g), nrow=ll, byrow=T)
  av_e <-apply(e,2,mean, na.rm=T)
  unconditional_se <- apply(d$wAIC*sqrt((ese^2)+(e-gg)^2),2,sum,na.rm=T)
  TerrySum <-apply(eW*wwAIC,2,sum, na.rm=T)
  
  #print(e)
  model_averaging <- cbind(Sum=TerrySum,Estimate=g,Unconditional_SE= unconditional_se, Lower95CI=g-unconditional_se*1.96, Upper95CI=g+unconditional_se*1.96)
        #print(model_averaging)
  attr(model_averaging,"heading") <- c("Model averaging\n", 
                                      paste("Response:", car:::responseName.default(object),"\n"))
  rownames(d) <- paste(1:length(var.names),". ", var.names, sep="")
  #d$"Select"[d$BIC==minBIC]<-"*"
  d$Select[d$Deviance == minDev] <- "*"
  d$"Select"[d$AIC==minAIC]<-"*"
	#d$"Select"[d$Cp==minCp]<-"*"
  attr(d,"heading") <- c("Model selection\n", 
                         paste("Response:", car:::responseName.default(object),"\n"))
  d<-list(d,model_averaging)
  attr(d, "class") <- c("modelSelection",class(d))
  
  #rownames(d) <- var.names
  #d$Select[d$Deviance == minDev] <- "*"
  #d$Select[d$AIC == minAIC] <- "*"
  #cat("Predictor names abbreviated\n")
  d
  #d$wAIC
  
}




allocations <- function (object) 
{
    r <- sqrt((c.2 * SS.3)/(c.3 * SS.2))
    r
}


bargraph <- function (dv, iv, grp = NULL, width = 1, space = 0.5, ymin, xmin = 0, 
    ymax, xlabs = c("N"), xmax, xlab, ylab1, ylab2, ylabs = c("N"), 
    xlabelsfontsize = 1, ylabelsfontsize = 1, ytitlefontsize, 
    ylabpos = 4, sublabel = c(""), plotting = T, legnd, errorbars = "None", 
    ...) 
{
    old.par <- par(no.readonly = T)
    on.exit(par(old.par))
    if (missing(ylab1)) 
        ylab1 <- deparse(substitute(dv))
    if (missing(ylab2)) 
        ylab2 <- ""
    if (missing(grp)) {
        means <- tapply(dv, iv, mean)
        stdev <- tapply(dv, iv, sd)
        n <- tapply(dv, iv, length)
        sterr <- stdev/sqrt(n)
        if (missing(xlab)) 
            xlab <- deparse(substitute(iv))
        legnd <- NULL
    }
    else {
        means <- tapply(dv, list(grp, iv), mean)
        stdev <- tapply(dv, list(grp, iv), sd)
        n <- tapply(dv, list(grp, iv), length)
        sterr <- stdev/sqrt(n)
        if (missing(xlab)) 
            xlab <- deparse(substitute(iv))
        if (missing(legnd)) 
            legnd <- c(levels(grp))
    }
    if (ymax == "NULL") {
        if (errorbars == "se") {
            y.heights <- means + sterr + ((means + sterr) * 0.05)
        }
        else if (errorbars == "sd") {
            y.heights = means + stdev + ((means + stdev) * 0.05)
        }
        else y.heights = means
        ymax <- max(y.heights)
    }
    if (ymin == "NULL") {
        ymin <- 0
    }
    y <- seq(ymin, ymax, by = ymax/5)
    yy <- seq(ymin, ymax, by = ymax/25)
    if (ylabs == "N") 
        ylabs <- seq(ymin, ymax, by = ymax/5)
    if (plotting == F) 
        barcolor <- c("white")
    if (plotting == F) 
        means <- c(0, 0)
    par(lwd = 2)
    x <- barplot(means, ylab = "", axes = FALSE, ylim = c(ymin, 
        ymax), width = width, axisnames = FALSE, beside = T, 
        xpd = F, legend.text = legnd, )
    if (errorbars == "se") {
        if (plotting == T) 
            arrows(x, means, x, means - sterr, angle = 90, length = 0.05, 
                lwd = 2)
        if (plotting == T) 
            arrows(x, means, x, means + sterr, angle = 90, length = 0.05, 
                lwd = 2)
    }
    else if (errorbars == "sd") {
        if (plotting == T) 
            arrows(x, means, x, means - stdev, angle = 90, length = 0.05, 
                lwd = 2)
        if (plotting == T) 
            arrows(x, means, x, means + stdev, angle = 90, length = 0.05, 
                lwd = 2)
    }
    axis(2, las = 1, lwd = 2, cex.axis = ylabelsfontsize, font = 2)
    if (missing(grp)) {
        x <- c(0.02, x)
        xlabs <- levels(iv)
    }
    else {
        xlabs <- levels(iv)
        xleg <- max(x)
        xx <- split(x, gl(length(levels(iv)), length(levels(grp))))
        xx <- sapply(xx, mean)
        x <- c(0.02, xx)
    }
    print(xlabs)
    xlabs <- c("", xlabs)
    axis(1, at = x, pos = ymin, labels = xlabs, lwd = 2, cex.axis = xlabelsfontsize, 
        tick = TRUE, font = 2)
    mtext(ylab1, 2, line = 3, font = 2, cex = 1.5)
    mtext(ylab2, 2, line = 2.2, font = 2, cex = 1.5)
    mtext(xlab, 1, line = 3, font = 2, cex = 1.5)
    par(old.par)
    x
}


dist.which <- function (dist, FUN = min) 
{
    FUN <- match.fun(FUN)
    command <- paste("a<-apply(as.matrix(dist),1,function(y) y==FUN(dist))", 
        sep = "")
    eval(base:::parse(text = command))
    colnames(a[, colSums(a, dims = 1) == 1])
}

in.boundary <- function (b.x, b.y, asc.x, asc.y) 
{
    if (b.x[1] < asc.x[1]) 
        return(FALSE)
    if (b.x[2] > asc.x[2]) 
        return(FALSE)
    if (b.y[1] < asc.y[1]) 
        return(FALSE)
    if (b.y[2] > asc.y[2]) 
        return(FALSE)
    return(TRUE)
}

 kver.to.shapefile <- function (ver, path = paste(getwd(), "/", sep = "", collapse = NULL)) 
{
    if (class(ver) != "kver") 
        stop("The first parameter must be a kver class object (output of getverticeshr)")
    a <- lapply(ver, data.frame)
    b <- lapply(1:length(ver), function(i) data.frame(Id = c(1:nlevels(ver[[i]][[1]])), 
        name = c(1:nlevels(ver[[i]][[1]]))))
    shape <- lapply(1:length(ver), function(i) convert.to.shapefile(a[[i]], 
        b[[i]], "Id", 5))
    out.files <- lapply(1:length(ver), function(i) write.shapefile(shape[[i]], 
        paste(path, names(ver[i]), sep = "", collapse = NULL), 
        arcgis = T))
    ifelse(length(ver) == 1, message <- (paste(length(ver), " shapefile has been saved in ", 
        path, sep = "", collapse = NULL)), message <- (paste(length(ver), 
        " shapefiles have been saved in ", path, sep = "", collapse = NULL)))
    message
}

lm.II <- function (formula, data, print = "all", type = "RMA", plot = F, 
    zero = F, ...) 
{
    cl <- match.call()
    mf <- form.call <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("lm")
    mf <- eval(mf, parent.frame())
    terms <- mf$terms
    x <- model.matrix(mf)
    var.names <- colnames(x)
    asgn <- attr(x, "assign")
    tl <- attr(mf$terms, "term.labels")
    response <- model.response(model.frame(mf))
    responseName <- car:::responseName.default(mf)
    tab <- NULL
    stats <- NULL
    if (is.null(weights(mf))) 
        wt <- rep(1, length(response))
    else wt <- weights(mf)
    if (length(var.names) > 2) {
        if (plot == T) {
            save.mfrow <- par(mfrow = mfrow(length(var.names)))
            on.exit(par(mfrow = save.mfrow))
        }
        for (i in var.names) {
            var <- which(i == var.names)
            res <- lsfit(x[, -var], cbind(x[, var], response), 
                wt = wt, intercept = FALSE)$residuals
            z <- m2regress(res[, 1], res[, 2], type = type, zero = zero)$coef[2, 
                ]
            zz <- matrix(z$stats[2, ], nrow = 1)
            rownames(zz) <- i
            z <- matrix(z$coef[2, ], nrow = 1)
            a <- mean(res[, 2]) - (z[1] * mean(res[, 1]))
            rownames(z) <- i
            tab <- rbind(tab, i = z)
            stats <- rbind(stats, i = zz)
            if (plot == T) {
                plot(res[, 1], res[, 2], xlab = paste(i, "| others"), 
                  ylab = paste(responseName, " | others"), las = 1, 
                  col = "black", pch = 16)
                abline(a, z[[1]], col = "red", lwd = 1.5)
            }
        }
    }
    else {
        var <- which(var.names[2] == var.names)
        z <- m2regress(x[, var], response, type = type, zero = zero)
        zz <- matrix(z$stats, nrow = 2)
        rownames(zz) <- c("(Intercept)", var.names[2])
        z <- matrix(z$coef, nrow = 2)
        rownames(z) <- c("(Intercept)", var.names[2])
        tab <- rbind(tab, z)
        stats <- rbind(stats, zz)
        if (plot == T) {
            plot(x[, var], response, xlab = paste(var.names[2]), 
                ylab = paste(responseName), las = 1, col = "black", 
                pch = 16)
            abline(z[1], z[2], col = "red", lwd = 1.5)
        }
    }
    #colnames(stats) <- c("St. Error", "t value", "Pr(>|t|)", 
    #    "CI 2.5%", "CI 97.5%")
    colnames(stats) <- c("Lower 95% CI", "Upper 95% CI")
    colnames(tab) <- c("Estimate")
    l <- list(coefficients = tab, statistics = stats, call = form.call, 
        terms = terms, model = mf$model, xlevels = mf$xlevels)
    class(l) <- "lm2"
    l
}


m2regress <- function (X, Y, type = "RMA", zero = F) 
{
    IV <- X
    DV <- Y
    ssx <- sum((IV - mean(IV))^2)
    ssy <- sum((DV - mean(DV))^2)
    sxy <- sum(((IV - mean(IV))) * ((DV - mean(DV))))
    r <- sxy/sqrt(ssx * ssy)
    r2 <- sxy/(ssx * ssy)
    df.DV <- length(DV) - 2
    n <- length(IV)
    s2xy <- abs((ssy - (sxy^2)/ssx)/(length(IV) - 2))
    mod <- lm(DV ~ IV)
    if (type == "OLS") {
        b <- sxy/ssx
        a <- mean(DV) - b * mean(IV)
        pp <- length(IV)
        sb <- sqrt((s2xy/ssx))
        sa <- sqrt(s2xy * (1/(length(IV) - 0) + (mean(IV)^2)/ssx))
        ci <- qt(0.975, length(IV) - 2) * sb
        bLower <- b - ci
        bUpper <- b + ci
        aCI <- qt(0.975, length(IV) - 2) * sa
        aLower <- a - aCI
        aUpper <- a + aCI
    }
    if (type == "MA") {
        b_ols <- sxy/ssx
        r2 <- (sxy/sqrt(ssx * ssy))^2
        d <- (b_ols^2 - r2)/(r2 * b_ols)
        b <- (d + sqrt(d^2 + 4))/2
        a <- mean(DV) - b * mean(IV)
        sb <- sqrt((s2xy/ssx))
        sa <- sqrt(s2xy * (1/(length(IV) - 0) + (mean(IV)^2)/ssx))
        B <- ((qt(0.975, n - 2)^2)) * ((1 - r2))/(n - 2)
        bUpper <- b * (sqrt(B + 1) + sqrt(B))
        bLower <- b * (sqrt(B + 1) - sqrt(B))
        aLower <- mean(DV) - (bUpper * mean(IV))
        aUpper <- mean(DV) - (bLower * mean(IV))
    }
    if (type == "rMA") {
        if (zero == F) {
            rIV <- (X - min(X))/(max(X) - min(X))
            rDV <- (Y - min(Y))/(max(Y) - min(Y))
        }
        if (zero == T) {
            rIV <- (X)/(max(X))
            rDV <- (Y)/(max(Y))
        }
        rssx <- sum((rIV - mean(rIV))^2)
        rssy <- sum((rDV - mean(rDV))^2)
        rsxy <- sum(((rIV - mean(rIV))) * ((rDV - mean(rDV))))
        rr <- rsxy/sqrt(rssx * rssy)
        rr2 <- rr^2
        rb_ols <- rsxy/rssx
        d <- (rb_ols^2 - rr2)/(rr2 * rb_ols)
        b <- (d + sqrt(d^2 + 4))/2
        rn <- length(rIV)
        rs2xy <- abs((rssy - (rsxy^2)/rssx)/(length(rIV) - 2))
        if (zero == F) 
            b <- b * ((max(Y) - min(Y))/(max(X) - min(X)))
        if (zero == T) 
            b <- b * ((max(Y))/(max(X)))
        a <- mean(DV) - b * mean(IV)
        sb <- sqrt((rs2xy/rssx))
        sa <- sqrt(rs2xy * (1/(length(rIV) - 0) + (mean(rIV)^2)/rssx))
        B <- ((qt(0.975, rn - 2)^2)) * ((1 - rr2))/(n - 2)
        bUpper <- b * (sqrt(B + 1) + sqrt(B))
        bLower <- b * (sqrt(B + 1) - sqrt(B))
        aLower <- mean(DV) - (bUpper * mean(IV))
        aUpper <- mean(DV) - (bLower * mean(IV))
    }
    if (type == "RMA") {
        r2 <- (sxy/sqrt(ssx * ssy))^2
        b <- sqrt(ssy/ssx)
        if (sxy < 0) 
            b <- b * -1
        a <- mean(DV) - (b * mean(IV))
        sb <- sqrt((s2xy/ssx))
        sa <- sqrt(s2xy * (1/(length(IV) - 0) + (mean(IV)^2)/ssx))
        n <- length(IV)
        ci <- qt(0.975, length(IV) - 2) * sb
        bLower <- b - ci
        bUpper <- b + ci
        B <- ((qt(0.975, n - 2)^2)) * ((1 - r2))/(n - 2)
        bUpper <- b * (sqrt(B + 1) + sqrt(B))
        bLower <- b * (sqrt(B + 1) - sqrt(B))
        aCI <- qt(0.975, length(IV) - 2) * sa
        aLower <- a - aCI
        aUpper <- a + aCI
        aLower <- mean(DV) - (bUpper * mean(IV))
        aUpper <- mean(DV) - (bLower * mean(IV))
    }
    coef <- matrix(NA, 2, 1, dimnames = list(c("(Intercept)", 
        "X"), c("Coef")))
   # stats <- matrix(NA, 2, 5, dimnames = list(c("(Intercept)", 
   #     "X"), c("SE", "t", "P|>t|", "2.5", "97.5")))
     stats <- matrix(NA, 2, 2, dimnames = list(c("(Intercept)", 
        "X"), c("2.5", "97.5")))
    coef[1, 1] <- a
    coef[2, 1] <- b
    #stats[1, 1] <- sa
    #stats[2, 1] <- sb
    #stats[1, 2] <- a/sa
    #stats[2, 2] <- b/sb
    #stats[1, 3] <- 2 * pt(abs(abs(a/sa)), length(DV) - 2, lower.tail = F)
    #stats[2, 3] <- 2 * pt(abs(b/sb), length(DV) - 2, lower.tail = F)
    #stats[1, 4] <- aLower
    #stats[1, 5] <- aUpper
    #stats[2, 4] <- bLower
    #stats[2, 5] <- bUpper
    stats[1, 1] <- aLower
    stats[1, 2] <- aUpper
    stats[2, 1] <- bLower
    stats[2, 2] <- bUpper
    z <- list(coefficients = coef, stats = stats)
    z
}




 mbuffer <- function (asc, sites, dist, fun) 
{
    require(adehabitat)
    require(sp)
    nm <- deparse(substitute(asc))
    bpu <- .buffer.point.unic(asc, dist)
    attr(bpu, "cellsize") <- attr(asc, "cellsize")
    ss <- coordinates(sites)
    cs <- attr(bpu, "cellsize")
    b <- list()
    cc <- data.frame(sites = NULL)
    cc.data <- NULL
    cc.sites <- NULL
    asc.x <- c(attr(asc, "xll"), attr(asc, "xll") + (attr(asc, 
        "cellsize") * (attr(asc, "dim")[1])))
    asc.y <- c(attr(asc, "yll"), attr(asc, "yll") + (attr(asc, 
        "cellsize") * (attr(asc, "dim")[2])))
    j <- 0
    for (i in 1:length(ss[, 1])) {
        xl <- c(ss[i, 1] - ((ncol(bpu) - 1) * cs)/2, ss[i, 1] + 
            ((ncol(bpu) - 1) * cs)/2)
        yl <- c(ss[i, 2] - ((nrow(bpu) - 1) * cs)/2, ss[i, 2] + 
            ((nrow(bpu) - 1) * cs)/2)
        if (in.boundary(xl + c(-dist, dist), yl + c(-dist, dist), 
            asc.x, asc.y)) {
            bb <- try(subsetmap(asc, xlim = xl, ylim = yl), TRUE)
            j <- j + 1
            b[[j]] <- bb * bpu
            if (fun == "mean") 
                cc.d <- mean(b[[j]], na.rm = T)
            if (fun == "sum") 
                cc.d <- sum(b[[j]], na.rm = T)
            if (fun == "percent") 
                cc.d <- sum(b[[j]], na.rm = T)/sum(bpu, na.rm = T)
            cc.data <- rbind(cc.data, DATA = cc.d)
            cc.sites <- rbind(cc.sites, SITES = as.character(sites[i, 
                ]$SITE))
        }
        else {
            DEM <- NA
            cc.data <- rbind(cc.data, DEM)
            cc.sites <- rbind(cc.sites, SITES = as.character(sites[i, 
                ]$SITE))
        }
    }
    cc <- data.frame(SITES = cc.sites, DATA = cc.data)
    dimnames(cc)[[2]] <- c("SITES", nm)
    return(cc)
}

 mcmcpvalue <- function (samp, fac) 
{
    s <- grep(fac, colnames(samp))
    samp <- as.matrix(samp)[, s]
    std <- backsolve(chol(var(samp)), cbind(0, t(samp)) - colMeans(samp), 
        transpose = TRUE)
    sqdist <- colSums(std * std)
    sum(sqdist[-1] > sqdist[1])/nrow(samp)
}

mdes.plot<-function (power = 0.8, sd, n) 
{
    library(R.utils)
    delta <- numeric(length(n))
    j <- 0
    if (length(power) == 1) 
        constant <- "power"
    else constant <- "sd"
    if (constant == "sd") 
        for (i in n) {
            j <- j + 1
            delta[j] <- power.t.test(n = i, sd = sd, power = power[1])$delta
        }
    else for (i in n) {
        j <- j + 1
        delta[j] <- power.t.test(n = i, sd = sd[1], power = power)$delta
    }
    if (getBuiltinOs.GString() == "windows") 
        windows(height = 4, width = 4)
    else X11(height = 4, width = 4)
    old_par <- par(mar = c(4, 4, 2, 1))
    if (constant == "sd") 
        plot(n, delta, axes = F, xlab = "", ylab = "", lwd = 2, 
            type = "n", main = substitute(paste("sd = ", sd), 
                list(sd = sd)))
    else plot(n, delta, axes = F, xlab = "", ylab = "", lwd = 2, 
        type = "n", main = substitute(paste("power = ", power), 
            list(power = power)))
    axis(1, tick = T, tcl = -0.2, lab = T, cex.axis = 1, lwd = 2)
    mtext("n", 1, font = 2, cex = 1.5, line = 3)
    axis(2, tick = T, tcl = -0.2, lab = T, cex.axis = 1, lwd = 2, 
        las = 1)
    mtext("Minimum detectable effect size", 2, font = 2, cex = 1.5, 
        line = 3)
    lines(c(0, max(n)), c(0.8, 0.8), lty = 2)
    box(bty = "l", lwd = 2)
    if (constant == "sd") {
        if (length(power) > 0) {
            kk <- 0
            for (k in power) {
                kk <- kk + 1
                delta <- numeric(length(n))
                j <- 0
                for (i in n) {
                  j <- j + 1
                  delta[j] <- power.t.test(n = i, sd = sd, power = power[kk])$delta
                }
                points(n, delta, col = "black", pch = 16, lwd = 2, 
                  type = "l", lty = kk)
            }
        }
        legend("topright", title = "Power", legend = power, lty = seq(1, 
            length(power), by = 1))
    }
    else {
        if (length(sd) > 0) {
            kk <- 0
            for (k in sd) {
                kk <- kk + 1
                delta <- numeric(length(n))
                j <- 0
                for (i in n) {
                  j <- j + 1
                  delta[j] <- power.t.test(n = i, sd = sd[kk], 
                    power = power)$delta
                }
                points(n, delta, col = "black", pch = 16, lwd = 2, 
                  type = "l", lty = kk)
            }
        }
        legend("topright", title = "sd", legend = sd, lty = seq(1, 
            length(sd), by = 1))
    }
}

 north.arrow <- function (s_x, s_y = NULL, height = 0.1, lwd = 1, cex = 1) 
{
    require(Hmisc)
    if (is.null(s_y)) {
        s_y <- s_x$y
        s_x <- s_x$x
    }
    dims <- par()$usr
    ht <- dims[2] - dims[1]
    wt <- dims[4] - dims[3]
    fac_y <- ht * height
    fac_x <- wt * (height * 1)
    b1 <- bezier(c(0, -0.01, -0.04), c(1, 0.55, 0.5))
    b2 <- bezier(c(-0.04, -0.02, 0), c(0.5, 0.484, 0.48))
    xs <- (c(b1$x, b2$x, 0, 0) * fac_x) + s_x
    ys <- (c(b1$y, b2$y, 0, 1) * fac_y) + s_y
    polygon(xs, ys, xlim = c(0, 1), ylim = c(0, 1), col = "black")
    lines((c(-0.04, 0.04) * fac_x) + s_x, (c(0.38, 0.38) * fac_y) + 
        s_y, lwd = 1, lend = 2)
    text((0 * fac_x) + s_x, (1 * fac_y) + s_y, "N", font = 2, 
        pos = 3, cex = cex)
}


 plotMeans <- function (response, factor1, factor2, error.bars = c("se", "sd", 
    "conf.int", "none"), level = 0.95, xlab = deparse(substitute(factor1)), 
    ylab = paste("mean of", deparse(substitute(response))), legend.lab = deparse(substitute(factor2)), 
    main = "Plot of Means", pch = 1:n.levs.2, lty = 1:n.levs.2, 
    col = palette()) 
{
    if (!is.numeric(response)) 
        stop("Argument response must be numeric.")
    xlab
    ylab
    legend.lab
    error.bars <- match.arg(error.bars)
    if (missing(factor2)) {
        if (!is.factor(factor1)) 
            stop("Argument factor1 must be a factor.")
        valid <- !(is.na(factor1) | is.na(response))
        factor1 <- factor1[valid]
        response <- response[valid]
        means <- tapply(response, factor1, mean)
        sds <- tapply(response, factor1, sd)
        ns <- tapply(response, factor1, length)
        if (error.bars == "se") 
            sds <- sds/sqrt(ns)
        if (error.bars == "conf.int") 
            sds <- qt((1 - level)/2, df = ns - 1, lower.tail = FALSE) * 
                sds/sqrt(ns)
        yrange <- if (error.bars != "none") 
            c(min(means - sds), max(means + sds))
        else range(means)
        levs <- levels(factor1)
        n.levs <- length(levs)
        plot(c(1, n.levs), yrange, type = "n", xlab = xlab, ylab = ylab, 
            axes = FALSE, main = main)
        points(1:n.levs, means, type = "b", pch = 16, cex = 2)
        box()
        axis(2)
        axis(1, at = 1:n.levs, labels = levs)
        if (error.bars != "none") 
            arrows(1:n.levs, means - sds, 1:n.levs, means + sds, 
                angle = 90, lty = 2, code = 3, length = 0.125)
    }
    else {
        if (!(is.factor(factor1) | is.factor(factor2))) 
            stop("Arguments factor1 and factor2 must be factors.")
        valid <- !(is.na(factor1) | is.na(factor2) | is.na(response))
        factor1 <- factor1[valid]
        factor2 <- factor2[valid]
        response <- response[valid]
        means <- tapply(response, list(factor1, factor2), mean)
        sds <- tapply(response, list(factor1, factor2), sd)
        ns <- tapply(response, list(factor1, factor2), length)
        if (error.bars == "se") 
            sds <- sds/sqrt(ns)
        if (error.bars == "conf.int") 
            sds <- qt((1 - level)/2, df = ns - 1, lower.tail = FALSE) * 
                sds/sqrt(ns)
        yrange <- if (error.bars != "none") 
            c(min(means - sds), max(means + sds))
        else range(means)
        levs.1 <- levels(factor1)
        levs.2 <- levels(factor2)
        n.levs.1 <- length(levs.1)
        n.levs.2 <- length(levs.2)
        if (n.levs.2 > length(col)) 
            stop(paste("Number of groups for factor2, ", n.levs.2, 
                ", exceeds number of distinct colours, ", length(col), 
                ".", sep = ""))
        plot(c(1, n.levs.1 + 1), yrange, type = "n", xlab = xlab, 
            ylab = ylab, axes = FALSE, main = main)
        box()
        axis(2)
        axis(1, at = 1:n.levs.1, labels = levs.1)
        for (i in 1:n.levs.2) {
            points(1:n.levs.1, means[, i], type = "b", pch = pch[i], 
                cex = 2, col = col[i], lty = lty[i])
            if (error.bars != "none") 
                arrows(1:n.levs.1, means[, i] - sds[, i], 1:n.levs.1, 
                  means[, i] + sds[, i], angle = 90, code = 3, 
                  col = col[i], lty = lty[i], length = 0.125)
        }
        x.posn <- n.levs.1 + 0.25
        y.posn <- sum(c(0.1, 0.9) * par("usr")[c(3, 4)])
        text(x.posn, y.posn, legend.lab, adj = c(0, -0.5))
        legend(x.posn, y.posn, levs.2, pch = pch, col = col, 
            lty = lty)
    }
    invisible(NULL)
}


power.plot <- function (delta, sd, n) 
{
    library(R.utils)
    power <- numeric(length(n))
    j <- 0
    if (length(delta) == 1) 
        constant <- "delta"
    else constant <- "sd"
    if (constant == "sd") 
        for (i in n) {
            j <- j + 1
            power[j] <- power.t.test(n = i, sd = sd, delta = delta[1])$power
        }
    else for (i in n) {
        j <- j + 1
        power[j] <- power.t.test(n = i, sd = sd[1], delta = delta)$power
    }
    if (getBuiltinOs.GString() == "windows") 
        windows(height = 4, width = 4)
    else X11(height = 4, width = 4)
    old_par <- par(mar = c(4, 4, 2, 1))
    if (constant == "sd") 
        plot(n, power, axes = F, xlab = "", ylab = "", lwd = 2, 
            type = "n", main = substitute(paste("sd = ", sd), 
                list(sd = sd)))
    else plot(n, power, axes = F, xlab = "", ylab = "", lwd = 2, 
        type = "n", main = substitute(paste("delta = ", delta), 
            list(delta = delta)))
    axis(1, tick = T, tcl = -0.2, lab = T, cex.axis = 1, lwd = 2)
    mtext("n", 1, font = 2, cex = 1.5, line = 3)
    axis(2, tick = T, tcl = -0.2, lab = T, cex.axis = 1, lwd = 2, 
        las = 1)
    mtext("power", 2, font = 2, cex = 1.5, line = 3)
    lines(c(0, max(n)), c(0.8, 0.8), lty = 2)
    box(bty = "l", lwd = 2)
    if (constant == "sd") {
        if (length(delta) > 0) {
            kk <- 0
            for (k in delta) {
                kk <- kk + 1
                power <- numeric(length(n))
                j <- 0
                for (i in n) {
                  j <- j + 1
                  power[j] <- power.t.test(n = i, sd = sd, delta = delta[kk])$power
                }
                points(n, power, col = "black", pch = 16, lwd = 2, 
                  type = "l", lty = kk)
            }
        }
        legend("bottomright", title = "Delta", legend = delta, 
            lty = seq(1, length(delta), by = 1))
    }
    else {
        if (length(sd) > 0) {
            kk <- 0
            for (k in sd) {
                kk <- kk + 1
                power <- numeric(length(n))
                j <- 0
                for (i in n) {
                  j <- j + 1
                  power[j] <- power.t.test(n = i, sd = sd[kk], 
                    delta = delta)$power
                }
                points(n, power, col = "black", pch = 16, lwd = 2, 
                  type = "l", lty = kk)
            }
        }
        legend("bottomright", title = "sd", legend = sd, lty = seq(1, 
            length(sd), by = 1))
    }
}


power.r.plot <- function (r = NULL, power = NULL, n = NULL) 
{
    library(pwr)
    library(R.utils)
    if (getBuiltinOs.GString() == "windows") 
        windows(height = 4, width = 4)
    else X11(height = 4, width = 4)
    old_par <- par(mar = c(4, 4, 2, 1))
    if (length(n) > 0) {
        power_val <- numeric(length(n))
        j <- 0
        for (i in n) {
            j <- j + 1
            power_val[j] <- pwr.r.test(n = i, r = r)$power
        }
        plot(n, power_val, axes = F, xlab = "", ylab = "", lwd = 2, 
            type = "l", main = substitute(paste("r = ", r), list(r = r)))
    }
    else {
        n_val <- numeric(length(power))
        j <- 0
        for (i in power) {
            j <- j + 1
            n_val[j] <- pwr.r.test(power = i, r = r)$n
        }
        plot(n_val, power, axes = F, xlab = "", ylab = "", lwd = 2, 
            type = "l", main = substitute(paste("r = ", r), list(r = r)))
    }
    axis(1, tick = T, tcl = -0.2, lab = T, cex.axis = 1, lwd = 2)
    mtext("n", 1, font = 2, cex = 1.5, line = 3)
    axis(2, tick = T, tcl = -0.2, lab = T, cex.axis = 1, lwd = 2, 
        las = 1)
    mtext("power", 2, font = 2, cex = 1.5, line = 3)
    lines(c(0, max(n)), c(0.8, 0.8), lty = 2)
    box(bty = "l", lwd = 2)
}


print.lm2 <- function (object) 
{
    object$coef
}


relatives <- function (term, names, factors) 
{
    is.relative <- function(term1, term2) {
        all(!(factors[, term1] & (!factors[, term2])))
    }
    if (length(names) == 1) 
        return(NULL)
    which.term <- which(term == names)
    (1:length(names))[-which.term][sapply(names[-which.term], 
        function(term2) is.relative(term, term2))]
}

 rmMatrix <- function (data, v.names = NULL, idvar, repeated, drop = NULL, 
    split = NULL) 
{
    r <- reshape(data, v.names = v.names, idvar = idvar, timevar = repeated, 
        drop = drop, direction = "wide")
    rnames <- names(r)
    i <- match(idvar, rnames)
    r <- r[, -i]
    if (!is.null(split)) {
        r <- split(r, eval(parse(text = paste("r$", split, sep = ""))))
        for (i in 1:length(r)) {
            ii <- match(split, rnames)
            r[[i]] <- r[[i]][, -ii]
            r[[i]] <- as.data.frame(r[[i]])
        }
    }
    else r <- as.data.frame(r)
    r
}


scale.bar <- function (s_x, s_y = NULL, len, ht, q = 2, cex = 1) 
{
    require(sp)
    require(rgdal)
    if (is.null(s_y)) {
        scale_start <- data.frame(x = s_x$x, y = s_x$y)
    }
    else {
        scale_start <- data.frame(x = 145.9692, y = -34.47851)
        scale_start <- data.frame(x = s_x, y = s_y)
    }
    coordinates(scale_start) = ~x + y
    proj4string(scale_start) <- CRS("+proj=longlat +ellps=WGS84")
    scale_start.utm <- spTransform(scale_start, CRS("+proj=utm +zone=55 +south"))
    s_xx <- attr(scale_start.utm, "coords")[, "x"]
    s_yy <- attr(scale_start.utm, "coords")[, "y"]
    x_start <- s_xx
    x_end <- s_xx + len
    y_start <- s_yy
    y_end <- s_yy + ht
    scale_end <- data.frame(x = x_end, y = y_end)
    coordinates(scale_end) = ~x + y
    proj4string(scale_end) <- CRS("+proj=utm +zone=55 +south")
    scale_end.latlong <- spTransform(scale_end, CRS("+proj=longlat +ellps=WGS84"))
    s_x2 <- attr(scale_end.latlong, "coords")[, "x"]
    s_y2 <- attr(scale_end.latlong, "coords")[, "y"]
    horiz <- s_x2 - attr(scale_start, "coords")[, "x"]
    vert <- s_y2 - attr(scale_start, "coords")[, "y"]
    m <- NULL
    x <- c(0, 1/q, 1/q, 0, 0)
    for (i in 1:q) {
        m[[i]] <- matrix(0, nrow = 5, ncol = 2)
        m[[i]][, 1] <- attr(scale_start, "coords")[, "x"] + (x * 
            (horiz))
        m[[i]][, 2] <- c(attr(scale_start, "coords")[, "y"], 
            attr(scale_start, "coords")[, "y"], attr(scale_start, 
                "coords")[, "y"] + vert, attr(scale_start, "coords")[, 
                "y"] + vert, attr(scale_start, "coords")[, "y"])
        x <- x + (1/q)
        if (i == (as.integer(i/2) * 2)) 
            polygon(m[[i]], col = "black")
        else polygon(m[[i]], col = "white")
        lx <- m[[i]][1, 1]
        ly <- m[[i]][1, 2]
        lines(c(lx, lx), c(ly, ly + (2 * vert)))
        lines(c(m[[i]][3, 1], m[[i]][3, 1]), c(ly, ly + (2 * 
            vert)))
        text(lx, ly + (cex * vert), (i - 1) * ((len/1000)/q), 
            pos = 3, cex = cex)
        text(m[[i]][3, 1], ly + (cex * vert), (i) * ((len/1000)/q), 
            pos = 3, cex = cex)
    }
    text(attr(scale_start, "coords")[, "x"] + (0.5 * (horiz)), 
        attr(scale_start, "coords")[, "y"] + (cex * vert), "Kilometres", 
        pos = 1, cex = cex)
}


splitInteractions <- function (split, factors, names, asgn, df.names) 
{
    ns <- names(split)
    for (i in unique(asgn)) {
        if (i == 0 || names[i + 1] %in% ns) 
            next
        f <- rownames(factors)[factors[, i] > 0]
        sp <- f %in% ns
        if (any(sp)) {
            if (sum(sp) > 1) {
                old <- split[f[sp]]
                nn <- f[sp]
                names(nn) <- nn
                marg <- lapply(nn, function(x) df.names[asgn == 
                  (match(x, names) - 1)])
                term.coefs <- strsplit(df.names[asgn == i], ":", 
                  fixed = TRUE)
                ttc <- sapply(term.coefs, function(x) x[sp])
                rownames(ttc) <- nn
                splitnames <- apply(expand.grid(lapply(old, names)), 
                  1, function(x) paste(x, collapse = "."))
                names(splitnames) <- splitnames
                tmp <- sapply(nn, function(i) names(old[[i]])[match(ttc[i, 
                  ], marg[[i]])])
                tmp <- apply(tmp, 1, function(x) paste(x, collapse = "."))
                new <- lapply(splitnames, function(x) match(x, 
                  tmp))
                split[[names[i + 1]]] <- new[sapply(new, function(x) length(x) > 
                  0)]
            }
            else {
                old <- split[[f[sp]]]
                marg.coefs <- df.names[asgn == (match(f[sp], 
                  names) - 1)]
                term.coefs <- strsplit(df.names[asgn == i], ":", 
                  fixed = TRUE)
                ttc <- sapply(term.coefs, function(x) x[sp])
                new <- lapply(old, function(x) seq(along = ttc)[ttc %in% 
                  marg.coefs[x]])
                split[[names[i + 1]]] <- new
            }
        }
    }
    split
}


summary.lm2 <- function (object) 
{
    tab <- cbind(object$coef, object$stat)
    ls <- list(Call = object$call, Coefficients = tab)
    ls
}


summaryRM <- function (object, rm, ...) 
{
    sum <- object
    eps <- epsi.GG.HF_old(rm)
    eps.GG <- eps[1]
    eps.HF <- eps[2]
    term.length <- length(sum[[2]][[1]]$Df)
    dfsG <- sum[[2]][[1]]$Df * eps.GG
    dfsH <- sum[[2]][[1]]$Df * eps.HF
    GG.p <- 1 - pf(sum[[2]][[1]]$"F value", dfsG, dfsG[term.length])
    HF.p <- 1 - pf(sum[[2]][[1]]$"F value", dfsH, dfsH[term.length])
    sum[[2]][[1]]$GG.P <- GG.p
    sum[[2]][[1]]$HF.P <- HF.p
    print(sum, digits = 5)
    cat(paste("\n", "GG epsilon: ", format(eps.GG, digits = 7), 
        "\n", sep = ""))
    cat(paste("HF epsilon: ", format(eps.HF, digits = 7), "\n"))
}

trim.blanks <- function (text) 
{
    gsub("^ ", "", gsub(" *$", "", text))
  }

legend.vfont<-
function (x, y = NULL, legend, fill = NULL, col = par("col"), 
    lty, lwd, pch, angle = 45, density = NULL, bty = "o", bg = par("bg"), 
    box.lwd = par("lwd"), box.lty = par("lty"), box.col = par("fg"), 
    pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd, xjust = 0, 
    yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5), 
    text.width = NULL, text.col = par("col"), merge = do.lines && 
        has.pch, trace = FALSE, plot = TRUE, ncol = 1, horiz = FALSE, 
    title = NULL, inset = 0, xpd, title.col = text.col, ...) 
{
    if (missing(legend) && !missing(y) && (is.character(y) || 
        is.expression(y))) {
        legend <- y
        y <- NULL
    }
    mfill <- !missing(fill) || !missing(density)
    if (!missing(xpd)) {
        op <- par("xpd")
        on.exit(par(xpd = op))
        par(xpd = xpd)
    }
    title <- as.graphicsAnnot(title)
    if (length(title) > 1) 
        stop("invalid title")
    legend <- as.graphicsAnnot(legend)
    n.leg <- if (is.call(legend)) 
        1
    else length(legend)
    if (n.leg == 0) 
        stop("'legend' is of length 0")
    auto <- if (is.character(x)) 
        match.arg(x, c("bottomright", "bottom", "bottomleft", 
            "left", "topleft", "top", "topright", "right", "center"))
    else NA
    if (is.na(auto)) {
        xy <- xy.coords(x, y)
        x <- xy$x
        y <- xy$y
        nx <- length(x)
        if (nx < 1 || nx > 2) 
            stop("invalid coordinate lengths")
    }
    else nx <- 0
    xlog <- par("xlog")
    ylog <- par("ylog")
    rect2 <- function(left, top, dx, dy, density = NULL, angle, 
        ...) {
        r <- left + dx
        if (xlog) {
            left <- 10^left
            r <- 10^r
        }
        b <- top - dy
        if (ylog) {
            top <- 10^top
            b <- 10^b
        }
        rect(left, top, r, b, angle = angle, density = density, 
            ...)
    }
    segments2 <- function(x1, y1, dx, dy, ...) {
        x2 <- x1 + dx
        if (xlog) {
            x1 <- 10^x1
            x2 <- 10^x2
        }
        y2 <- y1 + dy
        if (ylog) {
            y1 <- 10^y1
            y2 <- 10^y2
        }
        segments(x1, y1, x2, y2, ...)
    }
    points2 <- function(x, y, ...) {
        if (xlog) 
            x <- 10^x
        if (ylog) 
            y <- 10^y
        points(x, y, ...)
    }
    text2 <- function(x, y, ...) {
        if (xlog) 
            x <- 10^x
        if (ylog) 
            y <- 10^y
        text(x, y, ...)
    }
    if (trace) 
        catn <- function(...) do.call("cat", c(lapply(list(...), 
            formatC), list("\n")))
    cin <- par("cin")
    Cex <- cex * par("cex")
    if (is.null(text.width)) 
        text.width <- max(abs(strwidth(legend, units = "user", 
            cex = cex, ...)))
    else if (!is.numeric(text.width) || text.width < 0) 
        stop("'text.width' must be numeric, >= 0")
    xc <- Cex * xinch(cin[1L], warn.log = FALSE)
    yc <- Cex * yinch(cin[2L], warn.log = FALSE)
    if (xc < 0) 
        text.width <- -text.width
    xchar <- xc
    xextra <- 0
    yextra <- yc * (y.intersp - 1)
    ymax <- yc * max(1, strheight(legend, units = "user", cex = cex, ...)/yc)
    ychar <- yextra + ymax
    if (trace) 
        catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra, 
            ychar))
    if (mfill) {
        xbox <- xc * 0.8
        ybox <- yc * 0.5
        dx.fill <- xbox
    }
    do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 
        0))) || !missing(lwd)
    n.legpercol <- if (horiz) {
        if (ncol != 1) 
            warning("horizontal specification overrides: Number of columns := ", 
                n.leg)
        ncol <- n.leg
        1
    }
    else ceiling(n.leg/ncol)
    has.pch <- !missing(pch) && length(pch) > 0
    if (do.lines) {
        x.off <- if (merge) 
            -0.7
        else 0
    }
    else if (merge) 
        warning("'merge = TRUE' has no effect when no line segments are drawn")
    if (has.pch) {
        if (is.character(pch) && !is.na(pch[1L]) && nchar(pch[1L], 
            type = "c") > 1) {
            if (length(pch) > 1) 
                warning("not using pch[2..] since pch[1L] has multiple chars")
            np <- nchar(pch[1L], type = "c")
            pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
        }
    }
    if (is.na(auto)) {
        if (xlog) 
            x <- log10(x)
        if (ylog) 
            y <- log10(y)
    }
    if (nx == 2) {
        x <- sort(x)
        y <- sort(y)
        left <- x[1L]
        top <- y[2L]
        w <- diff(x)
        h <- diff(y)
        w0 <- w/ncol
        x <- mean(x)
        y <- mean(y)
        if (missing(xjust)) 
            xjust <- 0.5
        if (missing(yjust)) 
            yjust <- 0.5
    }
    else {
        h <- (n.legpercol + (!is.null(title))) * ychar + yc
        w0 <- text.width + (x.intersp + 1) * xchar
        if (mfill) 
            w0 <- w0 + dx.fill
        if (do.lines) 
            w0 <- w0 + (2 + x.off) * xchar
        w <- ncol * w0 + 0.5 * xchar
        if (!is.null(title) && (abs(tw <- strwidth(title, units = "user", 
            cex = cex, ...) + 0.5 * xchar)) > abs(w)) {
            xextra <- (tw - w)/2
            w <- tw
        }
        if (is.na(auto)) {
            left <- x - xjust * w
            top <- y + (1 - yjust) * h
        }
        else {
            usr <- par("usr")
            inset <- rep(inset, length.out = 2)
            insetx <- inset[1L] * (usr[2L] - usr[1L])
            left <- switch(auto, bottomright = , topright = , 
                right = usr[2L] - w - insetx, bottomleft = , 
                left = , topleft = usr[1L] + insetx, bottom = , 
                top = , center = (usr[1L] + usr[2L] - w)/2)
            insety <- inset[2L] * (usr[4L] - usr[3L])
            top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3L] + 
                h + insety, topleft = , top = , topright = usr[4L] - 
                insety, left = , right = , center = (usr[3L] + 
                usr[4L] + h)/2)
        }
    }
    if (plot && bty != "n") {
        if (trace) 
            catn("  rect2(", left, ",", top, ", w=", w, ", h=", 
                h, ", ...)", sep = "")
        rect2(left, top, dx = w, dy = h, col = bg, density = NULL, 
            lwd = box.lwd, lty = box.lty, border = box.col)
    }
    xt <- left + xchar + xextra + (w0 * rep.int(0:(ncol - 1), 
        rep.int(n.legpercol, ncol)))[1L:n.leg]
    yt <- top - 0.5 * yextra - ymax - (rep.int(1L:n.legpercol, 
        ncol)[1L:n.leg] - 1 + (!is.null(title))) * ychar
    if (mfill) {
        if (plot) {
            fill <- rep(fill, length.out = n.leg)
            rect2(left = xt, top = yt + ybox/2, dx = xbox, dy = ybox, 
                col = fill, density = density, angle = angle, 
                border = "black")
        }
        xt <- xt + dx.fill
    }
    if (plot && (has.pch || do.lines)) 
        col <- rep(col, length.out = n.leg)
    if (missing(lwd)) 
        lwd <- par("lwd")
    if (do.lines) {
        seg.len <- 2
        if (missing(lty)) 
            lty <- 1
        lty <- rep(lty, length.out = n.leg)
        lwd <- rep(lwd, length.out = n.leg)
        ok.l <- !is.na(lty) & (is.character(lty) | lty > 0)
        if (trace) 
            catn("  segments2(", xt[ok.l] + x.off * xchar, ",", 
                yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
        if (plot) 
            segments2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len * 
                xchar, dy = 0, lty = lty[ok.l], lwd = lwd[ok.l], 
                col = col[ok.l])
        xt <- xt + (seg.len + x.off) * xchar
    }
    if (has.pch) {
        pch <- rep(pch, length.out = n.leg)
        pt.bg <- rep(pt.bg, length.out = n.leg)
        pt.cex <- rep(pt.cex, length.out = n.leg)
        pt.lwd <- rep(pt.lwd, length.out = n.leg)
        ok <- !is.na(pch) & (is.character(pch) | pch >= 0)
        x1 <- (if (merge && do.lines) 
            xt - (seg.len/2) * xchar
        else xt)[ok]
        y1 <- yt[ok]
        if (trace) 
            catn("  points2(", x1, ",", y1, ", pch=", pch[ok], 
                ", ...)")
        if (plot) 
            points2(x1, y1, pch = pch[ok], col = col[ok], cex = pt.cex[ok], 
                bg = pt.bg[ok], lwd = pt.lwd[ok])
    }
    xt <- xt + x.intersp * xchar
    if (plot) {
        if (!is.null(title)) 
            text2(left + w/2, top - ymax, labels = title, adj = c(0.5, 
                0), cex = cex, col = title.col)
        text2(xt, yt, labels = legend, adj = adj, cex = cex, 
            col = text.col)
    }
    invisible(list(rect = list(w = w, h = h, left = left, top = top), 
        text = list(x = xt, y = yt)))
}


