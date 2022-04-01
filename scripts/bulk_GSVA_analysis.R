rm(list = ls())

library(GSVA)

source("scripts/R_functions/bulk_gsva_functions.v1.R")

############ load markers
if (TRUE)
{
  marker.res = readRDS("Data/Marker_Results.RDS")
  #marker.res = readRDS("Data/marker_res.RDS")
  marker.lst = lapply(marker.res,function(x) subset(x,FDR < 0.05 & summary.logFC > log2(1.2) & !is_feature_control_Mt))
  
  marker.lst = marker.lst[paste0("CLS",c(5,8,7,3,4,9,2))]
  names(marker.lst) = gsub("CLS","C",names(marker.lst))
  
  print(sapply(marker.lst,nrow))
}

################################## Analyze METABRIC data
if (TRUE)
{
  ############ load METABRIC data
  if (TRUE)
  {
    ## load data
    metabric.data = read.delim(file = "Data/bulk_data/metabric_tnbc.expr.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
    metabric.data = subset(metabric.data,!duplicated(Symbol))
    metabric.mat = as.matrix(metabric.data[,-c(1,2)])
    rownames(metabric.mat) = metabric.data$Symbol
    colnames(metabric.mat) = gsub("\\.","-",colnames(metabric.mat))
    
    metabric.cif = read.delim(file = "Data/bulk_data/metabric_clinical.survgroupUpdate.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
    metabric.cif = metabric.cif[match(colnames(metabric.mat),metabric.cif[[1]]),]
    metabric.cif$'chemo received' = grepl("CT",metabric.cif$Treatment)
    
    #metabric.cif = subset(metabric.cif,ER.Expr == "-" & PR.Expr == "-" & Her2.Expr == "-")
    
    common.sid = intersect(metabric.cif[[1]],colnames(metabric.mat))
    metabric.mat = metabric.mat[,common.sid]
    metabric.cif = metabric.cif[match(common.sid,metabric.cif[[1]]),]
    
    # add ESTIMATE score
    sample.list = readRDS("Data/bulk_data/METABRIC_SampleIDs_by_ESTIMATE_score_distribution.RDS")
    sample.list = sample.list[c(3,2,1)]
    sapply(sample.list,function(x,y) sum(x %in% y),y = colnames(metabric.mat))
    vec = rep(NA,nrow(metabric.cif));
    for (i in 1:length(sample.list)) vec[colnames(metabric.mat) %in% sample.list[[i]]] = names(sample.list)[i]
    metabric.cif$ESTIMATE.group = vec
    
    # add ESTIMATE data
    edf = read.delim(file = "Data/bulk_data/METABRIC.ESTIMATE.txt",sep = "\t",skip = 2,header = TRUE,stringsAsFactors = FALSE)
    edf = edf[,-2];
    edf.mat = as.matrix(edf[,-1]);rownames(edf.mat) = edf[[1]];
    colnames(edf.mat) = gsub("\\.","-",colnames(edf.mat))
    metabric.cif$ESTIMATE.score = edf.mat[3,match(metabric.cif[[1]],colnames(edf.mat))]
    
  }
  
  ############ Run GSVA and plot heatmap
  if (TRUE)
  {
    
    library(GSVA)
    
    ## run GSVA
    gsva.metabric = gsva(expr = metabric.mat,gset.idx.list= lapply(marker.lst,function(x) unique(x$hgnc_symbol)),method = "zscore")
    hout = hclust(dist(t(gsva.metabric)),"complete")
    
    cls = cutree(hout,k = 3) 
    table(cls)
    metabric.cif$GSVA.cluster = factor(cls[metabric.cif[[1]]])
    
    #### make heatmap
    library(ComplexHeatmap)
    library(circlize)
    col.df = metabric.cif[,c("lymph_nodes_positive","tnbc_chemo_disease_response",
                             "NOT_IN_OSLOVAL_P53_mutation_status","ESTIMATE.group","chemo received","GSVA.cluster")]
    colnames(col.df)[1:6] =c("Lymph Node","Disease Survival","P53 Mutation","ESTIMATE Group","chemo received","GSVA cluster")
    
    col.df[[2]][col.df[[2]] == "response"] = "Good"
    col.df[[2]][col.df[[2]] == "no_response"] = "Poor"
    
    hc = HeatmapAnnotation(df= col.df,col = list("Lymph Node" = colorRamp2(c(0,4,22),c("chartreuse","salmon","red")),
                                                 "Disease Survival" = c("Good" = "green","Poor" = "red"),
                                                 "ESTIMATE Group" = c("TP" = "yellow","HMTE" = "red","TE" = "pink"),
                                                 "chemo received" = c("TRUE" = "black","FALSE" = "grey"),
                                                 #"GSVA cluster" = c("1" = "yellow","2" = "red","3" = "chartreuse"),
                                                 "P53 Mutation" = c("MUT" = "cyan","WT" = "aquamarine3")))
    ht = Heatmap(matrix = gsva.metabric,name = "METABRIC\nGSVA",
                 col = colorRamp2(c(-24, 0, 12), c("green", "white", "red")),
                 #col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")),
                 show_column_names = FALSE,
                 #cluster_columns = hout,
                 column_split = factor(clus.mat[,kf-1]),
                 cluster_rows = FALSE,top_annotation = hc)
    
    draw(ht,heatmap_legend_side = "right",annotation_legend_side = "bottom")
  }
  
  ##### show KM curve 
  if (TRUE)
  {

    out = run_kmplot.wt.group(sdf = metabric.cif,group.col = "GSVA.cluster",
                              time.col = "DiseaseSpecificSurvival.time",
                              event.col = "DiseaseSpecificSurvival.status")
    out$kmplot = out$kmplot + 
      scale_colour_manual(values = c("1" = "yellow","2" = "red","3" = "chartreuse")) + 
      guides(colour = guide_legend(title = "GSVA\nCluster")) + 
      labs(x = "Time (Days)",y = "Disease-specific survival") + 
      theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18),
            legend.text = element_text(size = 14),legend.title = element_text(size = 16),
            plot.title = element_text(hjust = 0.5,size = 18),
            legend.position = "bottom",legend.direction = "horizontal")
    
    print(out$kmplot)
  
  }
}

################################## Analyze TCGA data
if (TRUE)
{
  library(GSVA)
  library(ComplexHeatmap)
  library(circlize)
  ### load TCGA data
  if (TRUE)
  {
    ## load data
    if (TRUE)
    {
      data.file = "Data/bulk_data/tcga_basal.expr.txt";
      cif.file = "Data/bulk_data/clinical.cart.2016-07-27T23_37_55.958352.biotab.survgroupUpdate.txt"
      annot.col = 3;
      symbol.col = 2
      time.col = "survival.overall.followup"
      event.col = "survival.overall.event"
      estf = "Data/bulk_data/TCGA.ESTIMATE.txt"
      
      
      data.df = read.delim(file = data.file,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
      data.mat = as.matrix(data.df[,-c(1:annot.col)]);
      rownames(data.mat) = data.df[[symbol.col]];
      colnames(data.mat) = sapply(strsplit(colnames(data.mat),"\\."),function(x) paste(x[1:3],collapse = "-"))
      cif.df = read.delim(file = cif.file,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
      cif.df = cif.df[match(colnames(data.mat),cif.df[[1]]),]
      
      
      # load ESTIMATE results
      estimate.res = read.delim(file = estf,sep = "\t",header = TRUE,stringsAsFactors = FALSE,skip = 2)
      emat = as.matrix(estimate.res[,-c(1:2)]);rownames(emat) = estimate.res[[1]]
      emat = t(emat)
      rownames(emat) = gsub("\\.","-",rownames(emat))
      emat = emat[,c("StromalScore","ImmuneScore","ESTIMATEScore")]
      est.cov = data.frame(sid= rownames(emat),as.data.frame(emat))
    }
    
    ## draw heatmap
    if (TRUE)
    {
      ## run GSVA
      gsva.tcga = gsva(expr = data.mat,gset.idx.list= lapply(marker.lst,function(x) unique(x$hgnc_symbol)),method = "zscore")
      #gsva.tcga = gsva(expr = data.mat,gset.idx.list= marker.lst,method = "zscore")
      hout = hclust(dist(t(gsva.tcga)),"complete")
      
      # get cluster
      cls = cutree(hout,k = 4) 
      table(cls)
      cif.df$'GSVA.cluster' = factor(cls)
      ## make heatmap
      col.df = cif.df[,c("tumor_stage","tnbc_chemo_overall_response")]
      colnames(col.df) =c("Stage","Overall Survival")
      
      col.df[[2]][col.df[[2]] == "response"] = "Good"
      col.df[[2]][col.df[[2]] == "no_response"] = "Poor"
      col.df$'GSVA cluster' = factor(cls)
      
      # make clinical annotation 
      stage.col = c("not reported" = "grey","stage i" = "cadetblue1","stage ia" = "cyan","stage ii" = "darksalmon",
                    "stage iia" = "darkorange1","stage iib" = "coral2",
                    "stage iiia" = "deeppink","stage iiib" = "brown1","stage iiic" = "brown",
                    "stage iv" = "red",
                    "stage x" = "cornsilk3")
      
      
      hc = HeatmapAnnotation(df = col.df,col = list(Stage = stage.col,'Overall Survival' = c("Good" = "green","Poor" = "red"),
                                                    'GSVA cluster' = c("1" = "red","3" = "green","2" = "palegreen1","4" = "pink")),
                             show_legend = c(FALSE,TRUE,TRUE))
      hest = HeatmapAnnotation(ESTIMATE = est.cov$ESTIMATEScore[match(colnames(gsva.tcga),est.cov[[1]])],
                               col = list(ESTIMATE = colorRamp2(c(0,5000,10000,20000),c("white","yellow","pink","red"))),
                               show_legend = FALSE)
      #est.anno = as.matrix(est.cov[,-c(1,4)])
      #colnames(est.anno) = c("Stromal","Immune")
      
      #hest = HeatmapAnnotation(ESTIMATE = anno_barplot(est.anno,gp = gpar(fill = c("deepskyblue","chocolate1"),col = c("deepskyblue","chocolate1"))),
      #                        show_legend = FALSE)
      
      ht = Heatmap(matrix = gsva.tcga,name = "TCGA\nGSVA",
                   col = colorRamp2(c(-10, 0, 10), c("green", "white", "red")),
                   #col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")),
                   show_column_names = FALSE,cluster_columns = hout,
                   cluster_rows = FALSE,top_annotation = hc,
                   bottom_annotation = hest)
      
      lgd.est = Legend(col_fun = colorRamp2(c(0,5000,10000,20000),c("white","yellow","pink","red")),
                       title = 'ESTIMATE\nScore',direction = "horizontal") 
      lgd.stage = Legend(labels = names(stage.col),legend_gp = gpar(fill = stage.col),title = "Stage",ncol = 3)
      
      draw(ht,annotation_legend_list = list(lgd.est,lgd.stage),heatmap_legend_side = "right",annotation_legend_side = "bottom")
    }
    
    ## draw KM curve
    if (TRUE)
    {
      out = run_kmplot.wt.group(sdf = cif.df,group.col = "GSVA.cluster",
                                time.col = time.col,
                                event.col = event.col)
      
      out$kmplot = out$kmplot + 
        #scale_colour_manual(values = c("1" = "red","3" = "green","2" = "palegreen1","4" = "pink")) +
        scale_colour_manual(values = c("1" = "red","3" = "yellow","2" = "green","4" = "pink")) + 
        guides(colour = guide_legend(title = "GSVA\nCluster")) + 
        labs(x = "Time (Days)",y = "Overall survival") + 
        theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18),
              legend.text = element_text(size = 14),legend.title = element_text(size = 16),
              plot.title = element_text(hjust = 0.5,size = 18),
              legend.position = "bottom",legend.direction = "horizontal")
      
      print(out$kmplot)
    }
    
  }
}