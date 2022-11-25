library(ggplot2)

# check KM plot
run_kmplot.wt.group <- function(sdf,group.col,time.col,event.col)
{
  require(ggplot2)
  require(ggkm)
  require(survival)
  
  # check best split of gene expression groups
  # if (is.null(qtiles)) qtiles = seq(0.1,0.9,0.1)
  tcol <- match(time.col,colnames(sdf))
  ecol <- match(event.col,colnames(sdf))
  gcol = match(group.col,colnames(sdf))
  
  lr.test <- survdiff(Surv(sdf[[tcol]],sdf[[ecol]]) ~ sdf[[gcol]])
  logrank.p = 1-pchisq(lr.test$chisq,1);
  
  str <- paste("Logrank P=",format(signif(logrank.p,3),scientific = TRUE),sep = "")
  
  kmobj <- ggplot(data = sdf,
                  aes_string(time = time.col,status = event.col,colour = group.col)) + 
    geom_km() + geom_kmticks() + labs(x = "Time",y = "Survival") + 
    labs(title = str) + 
    theme_bw() 
  output = list(kmplot = kmobj,logrank.p = logrank.p,survdiff.obj = lr.test)
  return(output)
}

get_clusterwise_dev <- function(x)
{
  if (is.matrix(x$obs)){
    otmp <- apply(x$obs,1,sum)
    etmp <- apply(x$exp,1,sum)
  }else{
    otmp <- x$obs
    etmp <- x$exp
  }
  vec = ((otmp-etmp)^2)/ diag(x$var)
  names(vec) = names(x$n)
  vec = vec*sign(otmp-etmp)
  return(vec)
}

model_ESTIMATE_score <- function(x,mu.o,lambda.o = c(0.5,0.5),k = 2,title.name = "ESTIMATE Score Distribution")
{
  sw.out = shapiro.test(x)
  
  plot.new()
  par(mfrow = c(1,2))
  qqnorm(x, pch = 1, frame = FALSE)
  qqline(x, col = "steelblue", lwd = 2)
  text(x = 1,y = 5000,label = paste0("Shapiro-Wilk p=",format(signif(sw.out$p.value,3),scientific = T)),cex = 1.5)
  
  # fit bimodal
  library(mixtools)
  set.seed(99)  
  x.mixmdl <- normalmixEM(x, lambda=lambda.o, mu=mu.o, k=k)
  #cutval = c(x.mixmdl$mu[1] + 3*x.mixmdl$sigma[1],x.mixmdl$mu[2])
  cutval = c(x.mixmdl$mu)
  
  plot(x.mixmdl,density = TRUE,whichplots = 2,main2 = title.name,xlab2 = "ESTIMATE Score")
  abline(v = cutval,col = "black",lty = "dotted")
  
  group.list = list("TP" = names(x)[x <= cutval[1]],"TE" = names(x)[x > cutval[1] & x <= cutval[2]],"HTME" = names(x)[x > cutval[2]])
  
  output = list(MixEM.result = x.mixmdl,cutoff.value = cutval,group.labels = group.list)
  return(output)
}