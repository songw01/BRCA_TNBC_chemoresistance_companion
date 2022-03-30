rm(list = ls())

library(scran)
library(scater)
library(slingshot)

source("scripts/R_functions/lineage_inference_functions.v1.R")
###########

sce = readRDS("C:/Users/songw01/Documents/BRCA_TNBC_scRNAseq_chemoresistance_v2/Data/SCE/SCE_object.final.RDS")

sce$exp.design = rename_func(paste(sce$clonal.persistence,sce$treatment,sep = "_"))
sce$exp.design = factor(sce$exp.design,levels = c("SB","SM","SA","RB","RM","RA"))

table(sce$exp.design,sce$identifier)
table(sce$cluster,sce$exp.design)

global.sres = run_slingshot(sce.subset = sce,start.clus = c("CLS5"))

lineage = slingLineages(SlingshotDataSet(global.sres[[1]]))

#### run Spearman correlation analysis
if (TRUE)
{
  pseudo.df = colData(global.sres$slingshot.res)[grepl("slingPseudotime",colnames(colData(global.sres$slingshot.res)))]
  
  lineage.cor = data.frame()
  for (i in 1:ncol(pseudo.df))
  {
    # load expression matrix and lineage specific pseudo-time data
    ii = !is.na(pseudo.df[[i]])
    cells = rownames(pseudo.df)[ii]
    mat = assay(global.sres$slingshot.res[,match(cells,colnames(global.sres$slingshot.res))],"corrected_logcounts")
    gene.annot = rowData(global.sres$slingshot.res)
    pseudo.time = pseudo.df[[i]][ii];
    names(pseudo.time) = cells
    
    # filter out genes dominated by zeros
    n.zero = rowSums(mat < 1E-320,na.rm = TRUE)
    ii = n.zero <= 0.8*(ncol(mat))
    mat.filter = mat[ii,]
    gene.annot = gene.annot[ii,]
    dim(mat.filter)
    
    # correlation analysis
    res = apply(mat.filter,1,function(x,y) do.call('c',cor.test(x =x,y = y,method = "spearman")[c("estimate","p.value")]),y = pseudo.time)
    res.df = data.frame(gene.annot,t(res),adj.p.value = p.adjust(res[2,],"bonferroni"))
    res.df = res.df[order(res.df$p.value),]
    res.df = subset(res.df,!is_feature_control_Mt)
    res.df$lineage = rep(colnames(pseudo.df)[i],nrow(res.df))
    lineage.cor = rbind.data.frame(lineage.cor,res.df);
    rm(mat,cells,gene.annot,pseudo.time,res.df)
  }
  lineage.cor$lineage = gsub("slingPseudotime_","lineage.",lineage.cor$lineage)
  
  ### mark top 500 correlated genes: stored into lin.data 
  lineage.cor$COR.ID = rep(NA,nrow(lineage.cor))
  lin.data = data.frame()
  for (i in c("lineage.1","lineage.2","lineage.3"))
  {
    df = subset(lineage.cor,lineage == i)
    df = df[order(df$estimate.rho,decreasing = T),]
    df$COR.ID[1:500] = "POS"
    lin.data = rbind.data.frame(lin.data,df)
    rm(df)
  }
  table(lin.data$COR.ID,lin.data$lineage)
  
  
}
