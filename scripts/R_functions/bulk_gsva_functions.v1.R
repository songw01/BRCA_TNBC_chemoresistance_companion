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
  output = list(kmplot = kmobj,logrank.p = logrank.p)
  return(output)
}
