rm(list = ls())

library(scran)
library(scater)
library(slingshot)

source("scripts/R_functions/lineage_inference_functions.v1.R")

########### load processed data
sce = readRDS("Data/SCE_object.final.RDS")

sce$exp.design = rename_func(paste(sce$clonal.persistence,sce$treatment,sep = "_"))
sce$exp.design = factor(sce$exp.design,levels = c("SB","SM","SA","RB","RM","RA"))

table(sce$exp.design,sce$identifier)
table(sce$cluster,sce$exp.design)

############# summarize cell composition
if (TRUE)
{
  cls.tick = paste("CLS",c(1,6,10,5,8,7,3,4,9,2),sep = "")
  # summarize clonal persistence
  require(reshape)
  tbl = table(rename_func(paste(sce$clonal.persistence,sce$treatment,sep = "_")),sce$cluster)
  cell.cnt = colSums(tbl)
  df = melt(tbl)
  colnames(df) = c("exp.design","cluster","count")
  df$count = df$count/cell.cnt[match(df$cluster,names(cell.cnt))]
  df$cell.group = gsub("[A-Z]$","",df$exp.design)
  # plot out
  require(ggplot2)
  df$exp.design = factor(df$exp.design,levels = c("SB","SM","SA","RB","RM","RA"))
  pobj.group = ggplot(data = df,aes(x = cluster,y = count,fill = exp.design,colour = cell.group)) + 
    geom_bar(stat = "identity",position = "stack") + 
    scale_fill_manual(values = c("SB" = "darkseagreen","SM" = "darkolivegreen1","SA" = "green",
                                 "RB" = "gold","RM" = "darkgoldenrod","RA" = "red"),
                      limits = c("SB","SM","SA","RB","RM","RA")) + 
    scale_colour_manual(values = c("S" = "green","R" = "red")) + 
    labs(x = "cell cluster",y = "Proportion") + 
    #scale_x_discrete(limits = paste("CLS",1:10,sep = "")) +
    scale_x_discrete(limits = cls.tick) + 
    theme_bw() +
    guides(fill = guide_legend(title = "Exp. Design"),colour = guide_legend(title = "Chemo-resistance\nGroup")) + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust =1 ),
          legend.position = "top",legend.direction = "horizontal")
  print(pobj.group)
}

############# Get expression markers
if (TRUE)
{
  # expression marker thresholds
  fc.cut = 1.2
  pval.cut = 0.05
  
  # now, get markers
  marker.out = findMarkers(x = sce,assay.type = "corrected_logcounts",groups = sce$cluster,row.data = rowData(sce)[,c("Geneid","ensembl_transcript_id","hgnc_symbol","composite.id","is_feature_control_Mt")])
  saveRDS(marker.out,file = "Data/Marker_Results.RDS")
  
  # update with CRISPRi hits
  chemo.sigs = read.geneSet("Data/brca_signature_library/chemo_resistance_compendium.gmt")
  chemo.sigs = chemo.sigs[grep("PTX|SUM|ZPOS|ZNEG",names(chemo.sigs))]
  zneg = Reduce("union",chemo.sigs[c("DOX_ZNEG","PTX_ZNEG")])
  zpos = Reduce("union",chemo.sigs[c("DOX_ZPOS","PTX_ZPOS")])
  intc = intersect(zneg,zpos)
  zneg = setdiff(zneg,intc)
  zpos = setdiff(zpos,intc)
  
  marker.out = lapply(marker.out,function(x,y) {out = x;out$is.ZNEG = out$jhgnc_symbol %in% y;return(out)},y = zneg)
  marker.out = lapply(marker.out,function(x,y) {out = x;out$is.ZPOS = out$hgnc_symbol %in% y;return(out)},y = zpos)
  
  #  get markers
  marker.lst = lapply(marker.out,function(x) {
    out = subset(x,(FDR < 0.05 & !is_feature_control_Mt & summary.logFC > log2(1.2)))
    out = out[order(out$p.value),]
    #out = out[1:min(c(100,nrow(out))),]
    out$composite.id
  })
  sapply(marker.lst,length)
  
  
}

############# trajectory analysis via slingshot
if (TRUE)
{
  global.sres = run_slingshot(sce.subset = sce,start.clus = c("CLS5"))
  
  lineage = slingLineages(SlingshotDataSet(global.sres[[1]]))
  
  # plot results in PCA space
  multiplot(plotlist = global.sres[[3]])
  
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
}

############# run PAM50, plot cell portions per cluster, then plot HER family expressions
if (TRUE)
{
  library(genefu)
  PAM50.col = c("Basal" = "red","Her2" = "cadetblue" ,"LumB" = "burlywood1","LumA" = "azure2","Normal" = "chartreuse2")
  
  data(vdxs)
  anno = as.data.frame(rowData(sce));
  anno_sub = subset(anno,hgnc_symbol %in% annot.vdxs$NCBI.gene.symbol)
  anno_sub = do.call('rbind.data.frame',lapply(split(1:nrow(anno_sub),factor(anno_sub$hgnc_symbol)),function(x,y) {ou = y[x,];ou[which.max(ou$ave.count),]},y = anno_sub))
  rownames(anno_sub) = anno_sub$Geneid
  anno_sub$EntrezGene.ID = annot.vdxs$EntrezGene.ID[match(anno_sub$hgnc_symbol,annot.vdxs$NCBI.gene.symbol)]
  colnames(anno_sub)[8] = "Gene.Symbol"
  colnames(anno_sub)[1] = "probe"
  mat = assay(sce,"corrected_logcounts");
  mat = mat[match(rownames(anno_sub),rownames(mat)),]
  
  pam50.res = molecular.subtyping(sbt.model="pam50", data=t(mat), 
                                  annot=anno_sub, do.mapping=TRUE)
  pam50.res$subtype[apply(pam50.res$subtype.proba,1,function(x) max(x) < 0.5)] = NA
  sce$PAM50 = pam50.res$subtype[match(names(pam50.res$subtype),colnames(sce))]
  
  tbl = table(sce$PAM50,sce$cluster)
  cell.cnt = colSums(tbl)
  df = melt(tbl)
  colnames(df) = c("PAM50","cluster","count")
  df$count = df$count/cell.cnt[match(df$cluster,names(cell.cnt))]
  
  # plot out
  require(ggplot2)
  
  pobj.PAM50 = ggplot(data = df,aes(x = cluster,y = count,fill = PAM50)) + 
    geom_bar(stat = "identity",position = "stack") + 
    scale_fill_manual(values = PAM50.col) + 
    labs(x = "cell cluster",y = "Proportion") + 
    #scale_x_discrete(limits = paste("CLS",1:10,sep = "")) +
    scale_x_discrete(limits = cls.tick) + 
    theme_bw() +
    guides(fill = guide_legend(title = "PAM50"),colour = guide_legend(title = "Chemo-resistance\nGroup")) + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust =1 ),
          legend.position = "top",legend.direction = "horizontal")
  print(pobj.PAM50)
  
  # get EGFR expressions
  pres = get_expression_violins(sce,genes = c("EGFR","ERBB2","ERBB3","ERBB4"),cluster.order = cls.tick)
  for (i in 1:length(pres)) pres[[i]] = pres[[i]] + scale_y_continuous(limits = c(0,4)) + theme(axis.text.y = element_text(size = 14))
  multiplot(plotlist = pres,cols = 1)
  
}

############# Generate marker heatmap: Figure 2B
if (TRUE)
{
  # get top markers to label on right side
  top.marker.lst <- lapply(marker.out,function(x) {
    out = subset(x,FDR < 0.05 & !is_feature_control_Mt & (is.ZNEG | is.ZPOS))
    
    tops = out$composite.id
    return(tops)
  })
  top.markers = Reduce("union",top.marker.lst)
  length(top.markers)
  
  # generate cluster tagged marker labels to improve the viz.
  top.marker.lab = rep(NA,length(top.markers))
  for (i in 1:length(top.markers))
  {
    cls.vec = gsub("CLS","",names(top.marker.lst)[sapply(top.marker.lst,function(x,y) any(x == y),y = top.markers[i])])
    top.marker.lab[i] = paste(gsub("\\|(.*)","",top.markers[i]),"(",paste(sort(as.integer(cls.vec)),collapse = ","),")",sep = "")
  }
  print(top.marker.lab)
  
  # mark them according to CRISPRi screening results
  top.marker.type = rep(NA,length(top.markers))
  top.marker.type[gsub("\\|(.*)$","",top.markers) %in% zneg] = "Antagonist"
  top.marker.type[gsub("\\|(.*)$","",top.markers) %in% zpos] = "Agonist"
  
  #### draw heatmap
  library(ComplexHeatmap)
  library(circlize)
  
  # marker gene space
  marker.genes = Reduce("union",marker.lst)
  
  # extract matrix
  expr.mat = assay(sce[match(gsub("^(.*)\\|","",marker.genes),rownames(sce)),],"corrected_logcounts")
  
  cls.order = paste("CLS",c(1,6,10,5,8,7,3,4,9,2),sep = "")
  cell.labs = c()
  for (i in 1:length(cls.order))
  {
    sce.cls = sce[match(gsub("^(.*)\\|","",marker.genes),rownames(sce)),sce$cluster == cls.order[i]]
    mat = assay(sce.cls,"corrected_logcounts")
    hobj = hclust(as.dist(sqrt(2*(1-cor(mat,method = "spearman")))),"complete")
    cell.labs = c(cell.labs,hobj$labels[hobj$order])
  }
  
  ii = match(cell.labs,colnames(expr.mat))
  
  # create cell annotation for columns
  cls.col = c("lightblue1","deepskyblue1","plum1","chartreuse3","tan3","coral2","mediumpurple","chocolate","palegreen4","firebrick")
  names(cls.col) = paste("CLS",1:10,sep = "")
  
  trt.col = c("SB" = "darkseagreen","SM" = "darkolivegreen1","SA" = "green",
              "RB" = "gold","RM" = "darkgoldenrod","RA" = "red")
  
  PAM50.col = c("Basal" = "red","Her2" = "cadetblue" ,"LumB" = "burlywood1","LumA" = "azure2","Normal" = "chartreuse2")
  
  df = colData(sce)[,c("cluster","exp.design","identifier","PAM50")]
  colnames(df) = c("Cluster","Treatment","Patient","PAM50")
  cell.anno = HeatmapAnnotation(df = df[ii,],name = "cell.annot",
                                col = list(Treatment = trt.col,
                                           Cluster = cls.col,
                                           PAM50 = PAM50.col),
                                show_legend = FALSE,na_col = "white",
                                annotation_legend_param = list(Cluster = gpar(ncol = 3),Treatment = gpar(ncol = 4),Patient = gpar(ncol = 3)))
  
  # add marker annotation for top 5 genes per cluster
  mark.col = rep("darkgreen",length(top.markers))
  mark.col[top.marker.type == "Agonist"] = "red"
  gene.anno = rowAnnotation(top.marker = anno_mark(at = match(gsub("^(.*)\\|","",top.markers),rownames(expr.mat)),
                                                   labels = top.marker.lab,side = "right",link_width = unit(10, "mm"),
                                                   labels_gp = list(col = mark.col)),show_legend = FALSE)
  
  # add celltype markers
  ct.anno = HeatmapAnnotation("EPCAM" = anno_barplot(x = assay(sce,"corrected_logcounts")["ENST00000263735.8",ii],ylim = c(0,5)),
                              "CD45" = anno_barplot(x = assay(sce,"corrected_logcounts")["ENST00000442510.7",ii],ylim = c(0,5)),
                              "ESR1" = anno_barplot(assay(sce,"corrected_logcounts")["ENST00000427531.6",ii],ylim = c(0,5)),
                              "PGR" = anno_barplot(assay(sce,"corrected_logcounts")["ENST00000325455.9",ii],ylim = c(0,5)),
                              "ERBB2" = anno_barplot(assay(sce,"corrected_logcounts")["ENST00000578199.5",ii],ylim = c(0,5)),
                              "Pseudotime1" = colData(global.sres[[1]])$slingPseudotime_1[ii],
                              "Pseudotime2" = colData(global.sres[[1]])$slingPseudotime_2[ii],
                              "Pseudotime3" = colData(global.sres[[1]])$slingPseudotime_3[ii],
                              annotation_legend_param = list("Pseudotime1" = list(direction = "horizontal"),
                                                             "Pseudotime2" = list(direction = "horizontal"),
                                                             "Pseudotime3" = list(direction = "horizontal")),
                              show_legend = TRUE)
  
  # heatmap color
  heat.col = colorRamp2(c(0, 1,11), c("white", "yellow", "red"))
  
  ### get main heatmap
  ht = Heatmap(matrix = expr.mat[,match(cell.labs,colnames(expr.mat))],col = heat.col,name = "log2(TPM+1)",
               cluster_columns = FALSE,
               clustering_method_rows = "complete",
               clustering_distance_rows = "spearman",
               show_column_names = FALSE,
               show_row_names = FALSE,
               top_annotation = cell.anno,
               right_annotation = gene.anno,
               bottom_annotation = ct.anno)
  
  ###### draw heatmap
  # legend
  heat.lgd = Legend(col_fun = heat.col,title = "log2(TPM+1)",direction = "horizontal",legend_width = unit(4, "cm"),title_position = "topcenter")
  trt.lgd = Legend(labels = names(trt.col)[c(1,4,2,5,3,6)],title = "Treatment",direction = "horizontal",ncol = 3,legend_gp = gpar(fill = trt.col[c(1,4,2,5,3,6)]),title_position = "topcenter")
  PAM50.lgd = Legend(labels = names(PAM50.col),title = "PAM50",ncol = 3,direction = "horizontal",legend_gp = gpar(fill = PAM50.col))
  pd = packLegend(heat.lgd,trt.lgd,PAM50.lgd,direction = "horizontal")
  
  # main heatmap
  draw(ht,show_heatmap_legend = FALSE, annotation_legend_list = pd,annotation_legend_side = "bottom")
  
  # add cluster numbers
  cls.id = paste("CLS",1:10,sep = "")
  #cell.names = as.hclust(cluster.dendro)$labels[as.hclust(cluster.dendro)$order]
  cell.names = cell.labs
  cls.names = sce$cluster[match(cell.names,colnames(sce))]
  for (i in 1:length(cls.id))
  {
    decorate_annotation(annotation = "Cluster",
                        code = {
                          grid.text(gsub("CLS","",cls.id[i]),y = 0.5,x = (mean(which(cls.names == cls.id[i])) - 0.5)/length(cls.names),
                                    gp = gpar(col = "black"),just = "center")
                        })
    
  }
  
}

