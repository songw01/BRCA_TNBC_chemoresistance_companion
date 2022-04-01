rm(list = ls())

library(fgsea)
library(ggplot2)
rename_signature <- function(x)
{
  x = gsub("\\.vs\\.","_vs_",gsub("pre","B",gsub("post","A",gsub("extinct","S",gsub("persistent","R",gsub("poor.prognosis.CoxPH","PoorProg",gsub("TumorEnriched","METABRIC:TE",gsub("RemoveHighESTIMATE","METABRIC:TP",gsub("GetHighESTIMATE","METABRIC:HTME",x)))))))))
  x = gsub("^B_","B:",gsub("^A_","A:",gsub("^R_","R:",x)))
  x = gsub("_vs_SUM159|_vs_SUM149","",gsub("PTX_resistance","TNBCLine:PTX_R_",x))
  return(x)
}

###### load resistance signatures
if (TRUE)
{
  # load pre-collected sigantures 
  surv.lst = lapply(read.geneSet("Data/brca_signature_library/brca_survival_compendium.GMT"),function(x) x[!is.na(x) & x != "NA"])
  metabric.lst = lapply(read.geneSet("Data/brca_signature_library/metabric_survival_compendium.GMT"),function(x) x[!is.na(x) & x != "NA"])
  chemo.lst = lapply(read.geneSet("Data/brca_signature_library/chemo_resistance_compendium.GMT"),function(x) x[!is.na(x) & x != "NA"])
  schemo.lst = lapply(read.geneSet("Data/brca_signature_library/intersected_chemo_signature.txt"),function(x) x[!is.na(x) & x != "NA"])
  schemores.lst = lapply(read.geneSet("Data/brca_signature_library/sc_chemoresistance_signature.txt"),function(x) x[!is.na(x) & x != "NA"])
  resist.lst = c(surv.lst,chemo.lst,schemo.lst,schemores.lst,metabric.lst)
  resist.lst <- lapply(resist.lst,function(x) unique(gsub("\\|(.*)$","",x)))
  
  # clusterwise resistance signature, then add to resist.lst
  pcut = 0.05;
  fc.cut = 2;
  
  cluster.resist = read.delim(file = "Data/CellCluster_DEG/Resistant_vs_Sensitive.limma.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  cluster.resist$DEG.ID = rep(NA,nrow(cluster.resist))
  cluster.resist$DEG.ID[cluster.resist$logFC > 0] = "UP"
  cluster.resist$DEG.ID[cluster.resist$logFC < 0] = "DN"
  cluster.resist$full.name = paste(cluster.resist$cluster,".R_vs_S",sep = "")
  cluster.resist.sig = subset(cluster.resist,adj.P.Val < pcut & abs(logFC) > log2(fc.cut))[,-c(1,2)]
  cls.res.sig = lapply(split(cluster.resist.sig[[1]],factor(paste(cluster.resist.sig$full.name,cluster.resist.sig$DEG.ID,sep = "_"))),function(x) unique(gsub("\\|(.*)$","",x)))
  
  resist.lst = c(resist.lst,cls.res.sig)
  
  # get marker list
  marker.res = readRDS("Data/Marker_Results.RDS")
  marker.lst = lapply(marker.res,function(x) unique(subset(x,FDR < 0.05 & summary.logFC > log2(1.2))$hgnc_symbol))
  
  ###
  ranking.sigs = resist.lst[names(resist.lst) %in% Reduce("union",ranking.signatures)]
  names(ranking.sigs) = rename_signature(names(ranking.sigs))
}

##### load combo treatment signatures
source("scripts/R_functions/enrichment_functions.R")
combo.deg = read.delim(file = "Data/PTX_Afatinib/DEG.limma.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
bg = union(Reduce("union",lapply(input.tables,function(x) as.character(x[[1]]))),combo.deg$gene.symbol)

#### GSEA on consensus signatures: Figure 5A, top left
if (TRUE)
{
  conc.sigs = readLines("Consensus_signature.txt")
  comp.id = unique(combo.deg$comparison);
  comp.id = comp.id[grep("DMSO",comp.id)]
  gsea.conc.results = data.frame()
  for (i in 1:length(comp.id))
  {
    tbl = subset(combo.deg,comparison == comp.id[i])
    tbl = subset(tbl,!(duplicated(gene.symbol) | adj.P.Val > 0.05))
    tstat = tbl$t;names(tstat) = tbl$gene.symbol
    gsea.conc = fgsea(pathways = list("Consensus" = conc.sigs), stats = tstat, nperm = 1000, minSize = 1, maxSize = Inf, nproc = 0,
                      gseaParam = 1, BPPARAM = NULL)
    gsea.conc$comparison = rep(comp.id[i],nrow(gsea.conc))
    gsea.conc.results = rbind.data.frame(gsea.conc.results,gsea.conc)
    rm(gsea.conc,tbl,tstat)
  }
  
  library(ggplot2)
  bnes = ggplot(data = gsea.conc.results,aes(x = comparison,y = NES,fill = NES)) + geom_bar(stat = "identity") + 
    scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + 
    theme_bw() + 
    theme(legend.text = element_text(size = 15),legend.title = element_text(size = 18),
          axis.title = element_text(size = 19),axis.text = element_text(size = 15))
  print(bnes)
}

#### GSEA on resistance signatures: Figure 5A
if (TRUE)
{
  pathways = c(ranking.sigs,marker.lst)
  pathways = lapply(pathways,function(x,y) intersect(x,y),y = unique(combo.deg$gene.symbol))
  pathways = pathways[sapply(pathways,length) > 10]
  sapply(pathways,length)
  
  comp.id = unique(combo.deg$comparison)
  gsea.results = data.frame()
  for (i in 1:length(comp.id))
  {
    tbl = subset(combo.deg,comparison == comp.id[i])
    tbl = subset(tbl,!(duplicated(gene.symbol) | adj.P.Val > 0.05))
    tstat = tbl$t;names(tstat) = tbl$gene.symbol
    gsea.res = fgsea(pathways = pathways, stats = tstat, nperm = 1000, minSize = 1, maxSize = Inf, nproc = 0,
                     gseaParam = 1, BPPARAM = NULL)
    gsea.res$comparison = rep(comp.id[i],nrow(gsea.res))
    gsea.results = rbind.data.frame(gsea.results,gsea.res)
    rm(tbl,tstat,gsea.res)
  }
  vec = rep("Resistance Signatures",nrow(gsea.results))
  vec[gsea.results$pathway %in% names(marker.lst)] = "Cluster Markers"
  gsea.results$sig.type = vec
  
  ## make heatmap
  library(reshape2)
  #pdata = subset(gsea.results,padj < 0.05 & !(pathway %in% names(conc.sigs)))
  pdata = subset(gsea.results,grepl(" - DMSO$",comparison))
  #pdata = gsea.results
  #pdata$set2_Name = rename_signature(x = as.character(pdata$set2_Name))
  pmat = acast(data = pdata,formula = pathway ~ comparison,value.var = "NES",fun.aggregate = function(x) min(x,na.rm = TRUE)) 
  pmat[is.infinite(pmat)] = 1
  hobj = hclust(dist(pmat),"complete")
  siglabs = hobj$labels[hobj$order]
  
  pdata$sig.class = rep("> 0.2",nrow(pdata))
  pdata$sig.class[pdata$padj < 0.2] = "0.05 > & < 0.2"
  pdata$sig.class[pdata$padj < 0.05] = "< 0.05"
  pdata$sig.class = factor(pdata$sig.class,levels = c("< 0.05","0.05 > & < 0.2","> 0.2"))
  pdata$pathway = gsub("^CLS","C",pdata$pathway)
  
  # get the significantly enriched signatures
  library(ggplot2)
  pobj = ggplot() + 
    #geom_tile(data = subset(combo.res,corrected.FET.pvalue < 0.05),aes(x = set1_Name,y = set2_Name,fill = -log10(corrected.FET.pvalue))) + 
    #scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = -log10(0.05)) + 
    geom_point(data = pdata,
               aes(x = comparison,y = pathway,size = -log10(padj),fill = NES,colour = sig.class),
               shape = 21,stroke = 1) + 
    scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + 
    scale_colour_manual(values = c("> 0.2" = "grey","0.05 > & < 0.2" = "magenta","< 0.05" = "red")) + 
    scale_size_continuous(range = c(1,7)) + 
    #scale_y_discrete(limits = siglabs) + 
    facet_grid( sig.type ~ .,scale = "free_y",space = "free_y") + 
    #scale_x_discrete(limits = drglabs) + 
    guides(fill = guide_colorbar(title = "NES"),size = guide_legend(title = "-log10(adj. P)",ncol = 2),
           colour = guide_legend(title = "Significance")) + 
    #facet_grid(. ~ set1_Name) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 14),axis.text.y = element_text(size = 15),
                       axis.title = element_blank(),strip.text = element_text(size = 15),
                       legend.position = "right",legend.direction = "vertical")
  print(pobj)
}

#### GSEA on hallmark signatures: Figure 5B
if (TRUE)
{
  gs = lapply(read.geneSet("Data/MSigDB_FrequentSets/h.all.v6.2.symbols.gmt"),function(x) x[-1])
  comp.id = unique(combo.deg$comparison);
  comp.id = comp.id[grep("DMSO",comp.id)]
  gsea.hrk.results = data.frame()
  for (i in 1:length(comp.id))
  {
    tbl = subset(combo.deg,comparison == comp.id[i])
    tbl = subset(tbl,!(duplicated(gene.symbol) | adj.P.Val > 0.05))
    tstat = tbl$t;names(tstat) = tbl$gene.symbol
    gsea.hrk = fgsea(pathways = gs, stats = tstat, nperm = 1000, minSize = 1, maxSize = Inf, nproc = 0,
                     gseaParam = 1, BPPARAM = NULL)
    gsea.hrk$comparison = rep(comp.id[i],nrow(gsea.hrk))
    gsea.hrk.results = rbind.data.frame(gsea.hrk.results,gsea.hrk)
    rm(gsea.conc,tbl,tstat)
  }
  
  ### make heatmap
  library(ggplot2)
  library(reshape2)
  sig.terms = unique(subset(gsea.hrk.results,padj < 0.05)$pathway)
  tbl = subset(gsea.hrk.results,pathway %in% sig.terms)
  tbl$pathway = gsub("^HALLMARK_","",tbl$pathway)
  
  # get the pathway limits by NES matrix
  pmat = acast(data = tbl,formula = pathway ~ comparison,value.var = "NES",fun.aggregate = function(x) mean(x,na.rm = T))
  hcls = hclust(dist(pmat),"complete")
  
  # mark significance class
  tbl$sig.class = rep("> 0.2",nrow(tbl))
  tbl$sig.class[tbl$padj < 0.2] = "0.05 > & < 0.2"
  tbl$sig.class[tbl$padj < 0.05] = "< 0.05"
  tbl$sig.class = factor(tbl$sig.class,levels = c("< 0.05","0.05 > & < 0.2","> 0.2"))
  tbl$class = "Hallmark Signatures"
  # generate plot with ggplot2
  hobj = ggplot() + 
    facet_grid(class ~ .) +
    #geom_tile(data = tbl,aes(x = comparison,y = pathway,fill = NES,alpha = -log10(padj))) + 
    #geom_point(data = subset(tbl,padj < 0.05),aes(x = comparison,y = pathway),colour = "black",shape = 19) + 
    geom_point(data = tbl,
               aes(x = comparison,y = pathway,size = -log10(padj),fill = NES,colour = sig.class),
               shape = 21,stroke = 1.5) + 
    scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0)  + 
    scale_colour_manual(values = c("> 0.2" = "grey","0.05 > & < 0.2" = "magenta","< 0.05" = "red")) + 
    scale_size_continuous(range = c(1,7)) + 
    scale_y_discrete(limits = hcls$labels[hcls$order]) + 
    #guides(size = guide_legend(title = "-log10(adj. P)"),colour = guide_legend(title = "Significance\n(adj. P)")) + 
    guides(size = FALSE,colour = FALSE,fill = FALSE) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45,vjust= 1,hjust = 1,size = 15),
          axis.text.y = element_text(size = 15),axis.title.y = element_blank(),
          axis.title.x = element_text(size = 18),
          strip.text = element_text(size = 16),
          legend.title = element_text(size = 17),legend.text = element_text(size = 14))
  print(hobj)
}
