library(scran)
library(scater)
library(slingshot)

rename_func <- function(x) 
{
  out = gsub("cellCluster","C",gsub("pre","B",gsub("mid","M",gsub("post","A",gsub("extinct_","S",gsub("persistent_","R",x))))))
  out =gsub("^c1_","M",out)
  return(out)
}
run_slingshot <- function(sce.subset,start.clus = NULL)
{
  exp.colmap = c("SB" = "darkseagreen","SM" = "darkolivegreen1","SA" = "green",
                 "RB" = "gold","RM" = "darkgoldenrod","RA" = "red")
  
  require(slingshot)
  require(RColorBrewer);
  require(ggrepel)
  
  # run slingshot
  if (!is.null(start.clus)) 
  {
    slingres = slingshot(data = sce.subset, clusterLabels = 'cluster', reducedDim = 'PCA',start.clus = start.clus)
  }else{
    slingres = slingshot(data = sce.subset, clusterLabels = 'cluster', reducedDim = 'PCA')
  }
  lineage = slingLineages(SlingshotDataSet(slingres))
  print(names(lineage))
  print(lineage)
  # add lineage information to slingres
  lineage.mat = matrix(0,nrow = ncol(slingres),ncol = length(lineage))
  colnames(lineage.mat) = names(lineage);
  rownames(lineage.mat) = colnames(slingres)
  for (i in 1:length(lineage))
  {
    lineage.mat[,i] = match(slingres$cluster,lineage[[i]])
  }
  
  lineage.df = as.data.frame(lineage.mat)
  colData(slingres) = cbind(colData(slingres),lineage.df)
  
  #### create a ggplot object
  require(ggplot2)
  require(grid)
  
  # PCA data
  pca.data = reducedDim(slingres,"PCA");
  cdat = as.data.frame(colData(slingres)[,(colnames(colData(slingres)) %in% c("cluster","exp.design")) | grepl("^slingPseudotime|^Lineage",colnames(colData(slingres)))])
  pca.data = cbind.data.frame(as.data.frame(pca.data),cdat)
  
  # lineage data
  lin.index = as.integer(gsub("^Lineage","",names(lineage)))
  #lin.col = c("black","dimgray","darkslategrey","indianred4")
  lin.col = rep("black",4)
  plst = vector("list",length(lin.index))
  names(plst) = names(lineage)
  for (i in 1:length(lin.index))
  {
    lineage.id = paste("Lineage",lin.index[i],sep = "")
    pt.id = paste("slingPseudotime_",lin.index[i],sep = "")
    
    # get cell plot
    pobj = ggplot() + 
      geom_point(data = pca.data,
                 aes_string(x = "PC1",y = "PC2",fill = pt.id,colour="exp.design"),
                 shape = 21,size = 3,alpha = 0.7,stroke = 1.4) + 
      scale_fill_gradient(low = "green",high = "red") + scale_colour_manual(values = exp.colmap)
    
    # get lineage data
    lineage.data = split(1:nrow(pca.data),factor(pca.data[[which(colnames(pca.data) == lineage.id)]]))
    lineage.data = lineage.data[order(as.integer(names(lineage.data)))]
    lineage.data = do.call('rbind',lapply(lineage.data,function(rr,m) colMeans(m[rr,]),m = as.matrix(pca.data[,c("PC1","PC2")])))
    colnames(lineage.data) = c("lineage.x","lineage.y")
    lineage.data = as.data.frame(lineage.data);
    lineage.data$cluster.id = lineage[[i]]
    
    # segment data
    seg.data = data.frame(x1 = lineage.data$lineage.x[1:(nrow(lineage.data)-1)],
                          y1 = lineage.data$lineage.y[1:(nrow(lineage.data)-1)],
                          x2 = lineage.data$lineage.x[2:(nrow(lineage.data))],
                          y2 = lineage.data$lineage.y[2:(nrow(lineage.data))])
    pobj = pobj + 
      geom_point(data = lineage.data,
                 aes_string(x = "lineage.x",y = "lineage.y"),colour = "black",size = 4,alpha = 0.3) + 
      geom_text_repel(data = lineage.data,
                      aes_string(x = "lineage.x",y = "lineage.y",label = "cluster.id")) + 
      geom_segment(data = seg.data,aes(x = x1,y = y1,xend = x2,yend = y2),colour = lin.col[i],alpha = 0.6,
                   size = 1,
                   arrow = arrow(length = unit(0.03, "npc"))) + 
      guides(fill = guide_colorbar(title = "Pseudo-time"),colour = guide_legend(title = "Study\nDesign")) + 
      #labs(title = paste(pid,":",lineage.id,sep = "")) + 
      theme_bw() + theme(legend.position = "bottom",legend.direction = "horizontal")
    plst[[i]] = pobj
    
  }
  
  output = list(slingshot.res = slingres,pca.res = pca.data,plot.obj = plst)
  
  return(output)
}

get_expression_violins <- function(sce,genes,cluster.order = NULL)
{
  genes = genes[genes %in% rowData(sce)$hgnc_symbol]
  plst = vector("list",length(genes))
  for (i in 1:length(genes))
  {
    gene_df = subset(rowData(sce),hgnc_symbol == genes[i])
    gene.id = subset(rowData(sce),hgnc_symbol == genes[i])$Geneid
    if (nrow(gene_df) > 1) gene.id = gene_df$Geneid[which.max(gene_df$total_counts)]
    if (length(gene.id) == 1)
    {
      pobj = plotExpression(sce, features=gene.id,x = "cluster",exprs_values = "corrected_logcounts")
      pobj = pobj + labs(y = paste(genes[i],sep = "")) + theme_bw() + 
        theme(strip.background = element_blank(),strip.text.x = element_blank(),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 13),
              axis.title.x = element_blank(),axis.title.y = element_text(size = 15))
      if (!is.null(cluster.order)) pobj = pobj + scale_x_discrete(limits = cluster.order,labels = gsub("^CLS","C",cluster.order))
      plst[[i]] = pobj
      names(plst)[i] = paste(genes[i],"\nexpression",sep = "")
    }
  }
  #plst = plst[!sapply(plst,is.null)]
  return(plst)
}

read.geneSet <- function(geneSet.file)
{
  gene.list <- readLines(geneSet.file)
  gene.list <- strsplit(gene.list,"\t")
  
  names(gene.list) <- sapply(gene.list,function(x) x[1])
  gene.list <- lapply(gene.list,function(x) x[2:length(x)])
  return(gene.list)
}

output.geneSet.file <- function(geneSet,outputfname)
{
  if (!is.list(geneSet)) stop("geneSet is not a list.")
  if (is.null(names(geneSet))) stop("names(geneSet) is not defined properly.")
  
  sink(outputfname)
  cat(paste(paste(names(geneSet),"\t",sapply(geneSet,function(x) paste(x,collapse = "\t")),sep = ""),collapse = "\n"))
  sink()
  
  return(0)
}
