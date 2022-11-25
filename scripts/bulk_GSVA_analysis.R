rm(list = ls())

library(GSVA)
library(patchwork)
source("scripts/R_functions/bulk_gsva_functions.v2.R")
source("scripts/R_functions/enrichment_functions.v2.R")

############ load markers
if (TRUE)
{
  nf = 500
  marker.res = readRDS("Data/Marker_Results.RDS")
  marker.res = lapply(marker.res,function(x) x[order(x$summary.logFC,decreasing = T),])
  
  marker.lst = lapply(marker.res,function(x) 
  {
    out = subset(x,FDR < 0.05 & summary.logFC > log2(1) & !is_feature_control_Mt)
    out = out[1:min(c(nrow(out),nf)),]
    return(out)
  })
  
  marker.lst = marker.lst[paste0("CLS",c(5,8,7,3,4,9,2))]
  names(marker.lst) = gsub("CLS","C",names(marker.lst))
  
  spec.sigs = lapply(marker.lst,function(x) unique(x$hgnc_symbol))
  
  print(sapply(spec.sigs,length))
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
    
    # load ESTIMATE results
    estf = "Data/bulk_data/METABRIC.ESTIMATE.txt"
    estimate.res = read.delim(file = estf,sep = "\t",header = TRUE,stringsAsFactors = FALSE,skip = 2)
    emat = as.matrix(estimate.res[,-c(1:2)]);rownames(emat) = estimate.res[[1]]
    emat = t(emat)
    rownames(emat) = gsub("\\.","-",rownames(emat))
    emat = emat[,c("StromalScore","ImmuneScore","ESTIMATEScore")]
    est.cov = data.frame(sid= rownames(emat),as.data.frame(emat))
    metabric.cif$ESTIMATE.score = emat[match(metabric.cif[[1]],rownames(emat)),3]
    
    ### run Shapiro-Wilk Test and fit biomodal normal
    x = metabric.cif$ESTIMATE.score;names(x) = metabric.cif[[1]]
    
    # customized: since TCGA shows single peak, fit single normal distribution
    library(MASS)
    set.seed(99)  
    tiff(file = "METABRIC_ESTIMATE.Mix_Model.tiff",res = 600,width = 7000,height = 3000)
    mixres = model_ESTIMATE_score(x = x,mu.o = c(5000,11000))
    dev.off()
    
    #sample.list = mixres$group.labels
    #sample.list = sample.list[c(3,2,1)]
    #sapply(sample.list,function(x,y) sum(x %in% y),y = colnames(metabric.mat))
    #vec = rep(NA,nrow(metabric.cif));
    #for (i in 1:length(sample.list)) vec[colnames(metabric.mat) %in% sample.list[[i]]] = names(sample.list)[i]
    #metabric.cif$ESTIMATE.group = vec
    
    #metabric.cif = subset(metabric.cif,ER.Expr == "-" & PR.Expr == "-" & Her2.Expr == "-")
    #metabric.cif = subset(metabric.cif,ESTIMATE.group %in% c("TE","TP"))
    #metabric.cif = subset(metabric.cif,ESTIMATE.score < (mixres$MixEM.result$mu[2] - mixres$MixEM.result$sigma[2]))
    #metabric.cif = subset(metabric.cif,ESTIMATE.score < (mixres$MixEM.result$mu[2]))
    metabric.cif = subset(metabric.cif,ESTIMATE.score < median(metabric.cif$ESTIMATE.score,na.rm = T))
    
    common.sid = intersect(metabric.cif[[1]],colnames(metabric.mat))
    metabric.mat = metabric.mat[,common.sid]
    metabric.cif = metabric.cif[match(common.sid,metabric.cif[[1]]),]
    
  }
  
  ############ Run GSVA and plot heatmap
  if (TRUE)
  {
    library(GSVA)
    library(fpc)
    
    ## run GSVA
    gsva.metabric = gsva(expr = metabric.mat,gset.idx.list= spec.sigs,method = "gsva")
    
    d = as.dist(sqrt(2*(1-cor(gsva.metabric))))
    #d = dist(t(gsva.metabric))
    hout = hclust(d,"complete")
    
    # evaluate clusters
    kvec = 2:10;
    eval.mat = matrix(0,nrow = length(kvec),ncol = 4)
    for (i in 1:length(kvec))
    {
      cls.i = cutree(hout,k = kvec[i])
      stat.out = cluster.stats(d = d,clustering = cls.i)
      eval.mat[i,] = do.call('c',stat.out[c("avg.silwidth","pearsongamma","dunn","dunn2")])
      
      if (i == 1) colnames(eval.mat) = c("avg.silwidth","pearsongamma","dunn","dunn2")
    }
    
    tiff(file = "METABRIC_GSVA.k_evaluation.tiff",res = 600,width = 3500,height = 3100)
    plot.new()
    par(mfrow = c(2,1),mar = c(3,5,0.5,0.5)) 
    plot(kvec,eval.mat[,1],xlab = "K",ylab = "Avg. Silhouette Width");abline(v = 5,col = "red")
    plot(kvec,eval.mat[,4],xlab = "K",ylab = "Adjusted Dunn's Index");abline(v = 5,col = "red")
    dev.off()
    
    cls = cutree(hout,k = 5) 
    #cls = cutree(hout,k = kvec[5]) 
    table(cls)
    metabric.cif$GSVA.cluster = factor(cls[metabric.cif[[1]]])
    est.cov$GSVA.cluster = cls[est.cov[[1]]]
    
    #### make heatmap
    library(ComplexHeatmap)
    library(circlize)
    col.df = metabric.cif[,c("lymph_nodes_positive","tnbc_chemo_disease_response",
                             "NOT_IN_OSLOVAL_P53_mutation_status")]
    colnames(col.df)[1:3] =c("Lymph Node","Disease Survival","P53 Mutation")
    
    col.df[[2]] = rep(NA,nrow(col.df))
    col.df[[2]][metabric.cif$DiseaseSpecificSurvival.status == 1 & metabric.cif$DiseaseSpecificSurvival.time <= 365*5] = "Poor"
    col.df[[2]][metabric.cif$DiseaseSpecificSurvival.status == 0 & metabric.cif$DiseaseSpecificSurvival.time > 365*5] = "Good"
    
    #### get top annotation bar
    hc = HeatmapAnnotation(df= col.df,col = list("Lymph Node" = colorRamp2(c(0,4,22),c("chartreuse","salmon","red")),
                                                 "Disease Survival" = c("Good" = "green","Poor" = "red"),
                                                 "P53 Mutation" = c("MUT" = "cyan","WT" = "aquamarine3")))
    
    #### get bottom annotation bar
    hest = HeatmapAnnotation(ESTIMATE = est.cov$ESTIMATEScore[match(colnames(gsva.metabric),est.cov[[1]])],
                             col = list(ESTIMATE = colorRamp2(c(0,5000,10000,20000),c("white","yellow","pink","red"))),
                             show_legend = TRUE)
    
    ht = Heatmap(matrix = gsva.metabric,name = "METABRIC\nGSVA",
                 #col = colorRamp2(c(-3, 0, 3), c("green", "white", "red")),
                 col = colorRamp2(c(-1,0,1), c("green", "white", "red")),
                 show_column_names = FALSE,
                 #cluster_columns = hout,
                 column_split = factor(cls[colnames(gsva.metabric)]),
                 cluster_rows = FALSE,top_annotation = hc,bottom_annotation = hest)
    
    tiff("METABRIC_GSVA.heatmap.tiff",res = 600,width = 6000,height = 3000)
    draw(ht,heatmap_legend_side = "left",annotation_legend_side = "bottom")
    dev.off()
  }
  
  ##### show KM curve 
  if (TRUE)
  {
    
    out = run_kmplot.wt.group(sdf = metabric.cif,group.col = "GSVA.cluster",
                              time.col = "DiseaseSpecificSurvival.time",
                              event.col = "DiseaseSpecificSurvival.status")
    
    # get clusterwise survival stats
    cls.survdiff = get_clusterwise_dev(out$survdiff.obj)
    names(cls.survdiff) = gsub("^(.*)=","",names(cls.survdiff))
    
    # merge with estimate score
    estimate.med = sapply(split(est.cov$ESTIMATEScore,factor(est.cov$GSVA.cluster)),median)
    cell.cluster.score = apply(gsva.metabric,1,function(x) sapply(split(x,factor(cls[colnames(gsva.metabric)])),median))
    
    # surv.table
    surv.tbl = table(cls,col.df$`Disease Survival`,useNA = "ifany")
    colnames(surv.tbl) = c("Good","Poor","NA")
    
    plot.data = data.frame(cluster.name = names(cls.survdiff),event.dev = cls.survdiff,cell.score = cell.cluster.score[,"C9"],
                           ESTIMATE.score = estimate.med,
                           as.data.frame.matrix(surv.tbl[names(cls.survdiff),]),cluster.size = rowSums(surv.tbl[names(cls.survdiff),]))
    colnames(plot.data) = gsub("\\.$","",colnames(plot.data))
    
    dev.obj = ggplot(data = plot.data) + 
      geom_bar(aes(x = cluster.name,y = event.dev,fill = event.dev),stat = "identity",colour = "black") + 
      scale_x_discrete(limits = names(cls.survdiff)[order(cls.survdiff)]) + 
      scale_fill_gradient2(low = "chartreuse",high = "red",mid = "white",midpoint = 0) +
      guides(fill = "none") + 
      labs(title = "Deviation from \nExpected #. Events (\u03C3)",x = "GSVA Cluster",y = "\u03C3") + theme_bw() + 
      theme(axis.title.y = element_text(size = 13),axis.text = element_text(size = 11),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 15,hjust = 0.5),
            legend.position = "bottom",legend.direction = "horizontal",
            legend.title = element_text(size = 17),legend.text = element_text(angle = 45,vjust = 1,hjust = 1,size = 15))
    
    
    
    out$kmplot = out$kmplot + 
      #scale_colour_manual(values = c("1" = "yellow","2" = "red","3" = "chartreuse")) + 
      guides(colour = guide_legend(title = "GSVA\nCluster")) + 
      labs(x = "Time (Days)",y = "Disease-specific survival") + 
      theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18),
            legend.text = element_text(size = 14),legend.title = element_text(size = 16),
            plot.title = element_text(hjust = 0.5,size = 15),
            legend.position = "bottom",legend.direction = "horizontal")
    
    surv.pobj = out$kmplot #+ inset_element(dev.obj, left = 0.67, bottom = 0.5, right = 0.99, top = 0.98)
    
    tiff("METABRIC_GSVA.KM_plot.tiff",res = 600,width = 5000,height = 3700)
    print(surv.pobj)
    dev.off()
  }
}

################################## Analyze TCGA data
if (TRUE)
{
  library(GSVA)
  library(ComplexHeatmap)
  library(circlize)
  library(fpc)
  library(MASS)
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
      cif.df$ESTIMATE.score = est.cov$ESTIMATEScore[match(cif.df[[1]],est.cov[[1]])]
      
      ### run Shapiro-Wilk Test and fit biomodal normal:show that it follow normality
      x = cif.df$ESTIMATE.score;names(x) = cif.df[[1]]
      tiff(file = "TCGA_ESTIMATE.Mix_Model.tiff",res = 600,width = 7000,height = 3000)
      plot.new()
      par(mfrow = c(1,2))
      sw.out = shapiro.test(x)
      qqnorm(x, pch = 1, frame = FALSE)
      qqline(x, col = "steelblue", lwd = 2)
      text(x = 1,y = 5000,label = paste0("Shapiro-Wilk p=",format(signif(sw.out$p.value,3),scientific = T)),cex = 1.5)
      
      fit = fitdistr(x = x,"normal")
      para <- fit$estimate
      
      hist(x, breaks = 10,prob = TRUE,ylim = c(0,0.0002))
      curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
      dev.off()
      
      #cif.df = cif.df[cif.df$ESTIMATE.score < para[1] - para[2],]
      cif.df = cif.df[cif.df$ESTIMATE.score < mean(cif.df$ESTIMATE.score),]
      
      data.mat = data.mat[,cif.df[[1]]]
    }
    
    ## draw heatmap
    if (TRUE)
    {
      ## run GSVA
      gsva.tcga = gsva(expr = data.mat,gset.idx.list= spec.sigs,method = "gsva")
      #gsva.tcga = gsva(expr = data.mat,gset.idx.list= marker.lst,method = "zscore")
      #d = dist(t(gsva.tcga))
      d = as.dist(sqrt(2*(1-cor(gsva.tcga))))
      
      hout = hclust(d,"complete")
      
      # evaluate clusters
      kvec = 2:10;
      eval.mat = matrix(0,nrow = length(kvec),ncol = 4)
      for (i in 1:length(kvec))
      {
        cls.i = cutree(hout,k = kvec[i])
        stat.out = cluster.stats(d = d,clustering = cls.i)
        eval.mat[i,] = do.call('c',stat.out[c("avg.silwidth","pearsongamma","dunn","dunn2")])
        
        if (i == 1) colnames(eval.mat) = c("avg.silwidth","pearsongamma","dunn","dunn2")
      }
      
      tiff(file = "TCGA_GSVA.k_evaluation.tiff",res = 600,width = 3500,height = 3100)
      plot.new()
      par(mfrow = c(2,1),mar = c(3,5,0.5,0.5)) 
      plot(kvec,eval.mat[,1],xlab = "K",ylab = "Avg. Silhouette Width");abline(v = 7,col = "red")
      plot(kvec,eval.mat[,4],xlab = "K",ylab = "Adjusted Dunn's Index");abline(v = 7,col = "red")
      dev.off()
      
      cls = cutree(hout,k = 7) 
      
      table(cls)
      cif.df$'GSVA.cluster' = factor(cls)
      
      ## check if ESTIMATE score and GSVA clusters are associated
      est.cov$GSVA.cluster = cls[match(est.cov$sid,names(cls))]
      
      ## make heatmap
      col.df = cif.df[,c("tumor_stage","tnbc_chemo_overall_response")]
      colnames(col.df) =c("Stage","Overall Survival")
      
      col.df[[2]] = rep(NA,nrow(col.df))
      col.df[[2]][cif.df$survival.overall.event == 1 & cif.df$survival.overall.followup <= 365*5] = "Poor"
      col.df[[2]][cif.df$survival.overall.event == 0 & cif.df$survival.overall.followup > 365*5] = "Good"
      table(col.df[[2]])
      
      col.df$'GSVA cluster' = factor(cls)
      
      # make clinical annotation 
      stage.col = c("not reported" = "grey","stage i" = "cadetblue1","stage ia" = "cyan","stage ii" = "darksalmon",
                    "stage iia" = "darkorange1","stage iib" = "coral2",
                    "stage iiia" = "deeppink","stage iiib" = "brown1","stage iiic" = "brown",
                    "stage iv" = "red",
                    "stage x" = "cornsilk3")
      
      
      hc = HeatmapAnnotation(df = col.df,col = list(Stage = stage.col,'Overall Survival' = c("Good" = "green","Poor" = "red")),
                             #'GSVA cluster' = c("1" = "red","3" = "green","2" = "palegreen1","4" = "pink")),
                             show_legend = c(FALSE,TRUE,TRUE))
      hest = HeatmapAnnotation(ESTIMATE = est.cov$ESTIMATEScore[match(colnames(gsva.tcga),est.cov[[1]])],
                               col = list(ESTIMATE = colorRamp2(c(0,5000,10000,20000),c("white","yellow","pink","red"))),
                               show_legend = FALSE)
      #est.anno = as.matrix(est.cov[,-c(1,4)])
      #colnames(est.anno) = c("Stromal","Immune")
      
      #hest = HeatmapAnnotation(ESTIMATE = anno_barplot(est.anno,gp = gpar(fill = c("deepskyblue","chocolate1"),col = c("deepskyblue","chocolate1"))),
      #                        show_legend = FALSE)
      
      ht = Heatmap(matrix = gsva.tcga,name = "TCGA\nGSVA",
                   #col = colorRamp2(c(-3, 0, 3), c("green", "white", "red")),
                   col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")),
                   column_split = factor(cls[colnames(gsva.tcga)]),
                   show_column_names = FALSE,
                   cluster_rows = FALSE,top_annotation = hc,
                   bottom_annotation = hest)
      
      lgd.est = Legend(col_fun = colorRamp2(c(0,5000,10000,20000),c("white","yellow","pink","red")),
                       title = 'ESTIMATE\nScore',direction = "horizontal") 
      lgd.stage = Legend(labels = names(stage.col),legend_gp = gpar(fill = stage.col),title = "Stage",ncol = 3)
      
      tiff("TCGA_GSVA.heatmap.tiff",res = 600,width = 6000,height = 3000)
      draw(ht,annotation_legend_list = list(lgd.est,lgd.stage),heatmap_legend_side = "left",annotation_legend_side = "bottom")
      dev.off()
    }
    
    ## draw KM curve
    if (TRUE)
    {
      out = run_kmplot.wt.group(sdf = cif.df,group.col = "GSVA.cluster",
                                time.col = time.col,
                                event.col = event.col)
      # get clusterwise survival stats
      cls.survdiff = get_clusterwise_dev(out$survdiff.obj)
      names(cls.survdiff) = gsub("^(.*)=","",names(cls.survdiff))
      
      plot.data = data.frame(cluster.name = names(cls.survdiff),event.dev = cls.survdiff,
                             cluster.size = table(cif.df$GSVA.cluster)[names(cls.survdiff)])
      colnames(plot.data) = gsub("\\.$","",colnames(plot.data))
      
      dev.obj = ggplot(data = plot.data) + 
        geom_bar(aes(x = cluster.name,y = event.dev,fill = event.dev),stat = "identity",colour = "black") + 
        scale_x_discrete(limits = names(cls.survdiff)[order(cls.survdiff)]) + 
        scale_fill_gradient2(low = "chartreuse",high = "red",mid = "white",midpoint = 0) +
        guides(fill = "none") + 
        labs(title = "Deviation from \nExpected #. Events (\u03C3)",x = "GSVA Cluster",y = "\u03C3") + theme_bw() + 
        theme(axis.title.y = element_text(size = 13),axis.text = element_text(size = 11),
              axis.title.x = element_blank(),
              plot.title = element_text(size = 15,hjust = 0.5),
              legend.position = "bottom",legend.direction = "horizontal",
              legend.title = element_text(size = 17),legend.text = element_text(angle = 45,vjust = 1,hjust = 1,size = 15))
      
      out$kmplot = out$kmplot + 
        #scale_colour_manual(values = c("1" = "yellow","2" = "red","3" = "chartreuse")) + 
        guides(colour = guide_legend(title = "GSVA\nCluster")) + 
        labs(x = "Time (Days)",y = "Overall survival") + 
        theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18),
              legend.text = element_text(size = 14),legend.title = element_text(size = 16),
              plot.title = element_text(hjust = 0.5,size = 15),
              legend.position = "bottom",legend.direction = "horizontal")
      
      surv.pobj = out$kmplot #+ inset_element(dev.obj, left = 0.6, bottom = 0.03, right = 0.98, top = 0.55)
      
      tiff("TCGA_GSVA.KM_plot.tiff",res = 600,width = 5000,height = 3700)
      print(surv.pobj)
      dev.off()
    }
    
  }
}