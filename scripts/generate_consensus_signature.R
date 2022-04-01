rm(list = ls())

# signature fold cutoffs
pcut = 0.05;
fc.cut = 2;

################################# Derive TNBC chemoresistance consensus signature 
##### read signatures from tables
if (TRUE)
{
  library(readxl)
  #### in vitro signature
  invitro = read_excel("Data/brca_signature_tables/Signature_Information_Table.xlsx",skip = 2,sheet = 7)
  invitro$DEG.ID[invitro$DEG.ID == "down-regulated"] = "DN"
  invitro$DEG.ID[invitro$DEG.ID == "up-regulated"] = "UP"
  invitro$comparison[invitro$comparison == "Paclitaxel resistant vs parental"] = "PTXR_vs_PTXP"
  invitro$group[invitro$group == "BT20, SUM149, MDA-MB-231, MDA-MB-436 and MDA-MB-468"] = "TNBC.Cells"
  invitro$full.name = paste(invitro$group,"_",invitro$comparison,sep = "")
  
  invitro.sig = subset(invitro,abs(logFC) > log2(fc.cut) & adj.P.Val < pcut)
  
  #### 
  ## single cell signature: chemo-response
  chem.resp = read_excel("Data/brca_signature_tables/Signature_Information_Table.xlsx",sheet = 4,skip = 2)
  
  # summarize stats to merge cross patient data
  genes = unique(chem.resp$Geneid)
  idx = sapply(split(1:nrow(chem.resp),factor(chem.resp$Geneid)),function(x,y) x[which.max(y$PValue[x])],y = chem.resp)
  chem.resp.merged = chem.resp[idx,]
  chem.resp.merged$full.name = "R.A_vs_B"
  chem.resp.merged$DEG.ID = rep(NA,nrow(chem.resp.merged))
  chem.resp.merged$DEG.ID[chem.resp.merged$logFC > 0] = "UP"
  chem.resp.merged$DEG.ID[chem.resp.merged$logFC < 0] = "DN"
  
  chem.resp.merged.sig = subset(chem.resp.merged,abs(logFC) > log2(fc.cut) & FDR < pcut)
  
  ## single cell signature: chemo-resistance
  chem.resist = read_excel("Data/brca_signature_tables/Signature_Information_Table.xlsx",sheet = 2,skip = 2)
  chem.resist$full.name = paste(chem.resist$group,".R_vs_S",sep = "")
  chem.resist$DEG.ID = rep(NA,nrow(chem.resist))
  chem.resist$DEG.ID[chem.resist$logFC > 0] = "UP"
  chem.resist$DEG.ID[chem.resist$logFC < 0] = "DN"
  chem.resist.sig = subset(chem.resist,abs(logFC) > log2(fc.cut) & FDR < pcut)
  
  ## cell cluster specific signature
  cluster.resist = read.delim(file = "Data/brca_signature_tables/Resistant_vs_Sensitive.limma.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  cluster.resist$DEG.ID = rep(NA,nrow(cluster.resist))
  cluster.resist$DEG.ID[cluster.resist$logFC > 0] = "UP"
  cluster.resist$DEG.ID[cluster.resist$logFC < 0] = "DN"
  cluster.resist$full.name = paste(cluster.resist$cluster,".R_vs_S",sep = "")
  cluster.resist.sig = subset(cluster.resist,adj.P.Val < pcut & abs(logFC) > log2(fc.cut))
  
  ###### Now, grab signatures per table
  input.tables = list(inVitro = invitro[,-1],sc.chemo.response = chem.resp.merged.sig,sc.chemo.resistant = chem.resist.sig,sc.cluster.resistant = cluster.resist.sig[,-c(1,2)])
}

##### derive consensus signature
if (TRUE)
{
  # extract signatures
  sig.lst = lapply(input.tables,function(x) lapply(split(x[[1]],factor(paste(x$full.name,x$DEG.ID,sep = "_"))),function(x) unique(gsub("\\|(.*)$","",x))))
  
  # now, get DEG concensus matrix
  up.lst = lapply(sig.lst,function(x) x[grep("_UP$",names(x))]);names(up.lst) = NULL
  up.lst = do.call('c',up.lst)
  
  dn.lst = lapply(sig.lst,function(x) x[grep("_DN$",names(x))]);names(dn.lst) = NULL
  dn.lst = do.call('c',dn.lst)
  
  deg.lst = list(UP = up.lst,DN = dn.lst)
  
  deg.matrix = lapply(deg.lst,function(x) {
    bg = setdiff(Reduce("union",x),NA)
    mat = do.call('cbind',lapply(x,function(y,z) z %in% y,z = bg))
    rownames(mat) = bg;
    colnames(mat) = names(x)
    return(mat)
  })
  
  # get counts
  conc.counts = lapply(deg.matrix,function(x) rowSums(x,na.rm = TRUE))
  
  count.vec = sort(Reduce("union",conc.counts))
  
  conc.signature = lapply(count.vec,function(x,y,z) list(UP = names(y)[y >= x],DN = names(z)[z >= x]),y = conc.counts$UP,z = conc.counts$DN)
  conc.overlap = lapply(conc.signature,function(x) Reduce("intersect",x))
  
  conc.summary = do.call('rbind',lapply(conc.signature,function(x) sapply(x,length)[c("UP","DN")]))
  conc.summary = data.frame(cutoff = count.vec,as.data.frame(conc.summary),intersect = sapply(conc.overlap,length))
  
  output = list(summary = conc.summary,signature = conc.signature)
  
  plot.new()
  with(output$summary, plot(cutoff, UP, col="black", 
                            ylab = "size of union set",xlab = "Concensus cutoff",cex.axis = 1.2,cex.lab = 1.5))
  abline(v = 4,col = "red")
  
  # clean up up-regulated signatures: remove down-regulated signature hits
  output$signature.no.int = lapply(output$signature,function(x) lapply(x,function(p,q) setdiff(p,q),q = Reduce("intersect",x)))
  
  cat("List concensus signatures:\n")
  conc.sigs = output$signature.no.int[[4]]$UP
  print(conc.sigs)
  
  sink("Consensus_signature.txt")
  cat(paste(conc.sigs,collapse = "\n"))
  cat("\n")
  sink()
}

##### derive summarized fold change
if (TRUE)
{
  ### add fold change info for final signature
  require(reshape2)
  conc.sig= output$signature.no.int[[4]]
  
  sig.res = vector("list",length(input.tables));names(sig.res) = names(input.tables)
  for (i in 1:length(input.tables))
  {
    tbl = as.data.frame(input.tables[[i]]);colnames(tbl) = gsub(" ","\\.",colnames(tbl))
    mat.res = vector("list",length(conc.sig));names(mat.res) = names(conc.sig)
    for (j in 1:length(conc.sig))
    {
      sig = conc.sig[[j]]
      ttbl = subset(tbl,DEG.ID == names(conc.sig)[j])
      ttbl = ttbl[ttbl[[1]] %in% sig,]
      
      if (nrow(ttbl) > 0)
      {
        mat = cbind(acast(data = ttbl,formula = as.formula(paste(colnames(ttbl)[1]," ~ full.name",sep = "")),value.var = "logFC",fun.aggregate = function(x) median(x,na.rm = TRUE)))
        mat = cbind(mat[match(sig,rownames(mat)),])
        rownames(mat) = sig
        mat.res[[j]] = mat;
        
      }
      rm(ttbl,sig,mat)
    }
    sig.res[[i]] = mat.res
    rm(mat.res)
  }
  
  fc.results = list(UP = do.call('cbind',lapply(sig.res,function(x) x$UP)),DN = do.call('cbind',lapply(sig.res,function(x) x$DN)))
  fc.summary = list(UP = data.frame(gene = rownames(fc.results$UP),logFC.summarized = apply(fc.results$UP,1,function(x,y) median(x[which(x >= y)],na.rm = TRUE),y = log2(fc.cut)),DEG.ID = rep("UP",nrow(fc.results$UP))),
                    DN = data.frame(gene = rownames(fc.results$DN),logFC.summarized = apply(fc.results$DN,1,function(x,y) median(x[which(x <= y)],na.rm = TRUE),y = log2(fc.cut)),DEG.ID = rep("DN",nrow(fc.results$DN))))
  fc.summary = do.call('rbind.data.frame',fc.summary)
  fc.summary = subset(fc.summary,!is.na(logFC.summarized))
  rownames(fc.summary) = NULL
  
  
  ### output finalized table
  write.table(fc.summary,file = "consensus_signature.summarized_lfc.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
}

##### run MSigDB enrichment analysis
if (TRUE)
{
  ### run MSigDB
  cat("- MSigDB enrichments...\n")
  source("scripts/R_functions/enrichment_functions.R")
  
  # set up 
  gmt.folder = "Data/MSigDB_FrequentSets"
  min.size = 5;
  degsig = split(fc.summary$gene,factor(fc.summary$DEG.ID))
  
  library(SingleCellExperiment)
  sce = readRDS("Data/SCE_object.final.RDS")
  bg = unique(rowData(sce)$hgnc_symbol)
  
  # get respective gmt files
  gmt.files <- list.files(path = gmt.folder,full.names = TRUE,pattern = "\\.gmt$")
  names(gmt.files) <- gsub("^(.*)/|\\.gmt$","",gmt.files)
  gmt.sets <- lapply(gmt.files,function(x) lapply(read.geneSet(x),function(y) y[-1]))
  
  # run MSigDB enrichments
  res.df <- data.frame()
  for (i in 1:length(gmt.sets))
  {
    df <- perform.AllPairs.FET(geneSets1 = degsig,geneSets2 = gmt.sets[[i]],background = bg,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
    df$database = rep(names(gmt.sets)[i],nrow(df))
    df$corrected.FET.pvalue <- p.adjust(df$FET_pvalue,"BH")
    res.df <- rbind.data.frame(res.df,df);rm(df)
  }
  res.df$corrected.FET.pvalue = p.adjust(res.df$FET_pvalue,"BH")
  res.df = res.df[order(res.df$FET_pvalue),]
  
  # write results
  write.table(res.df,file = "Concensus.MSigDB_FET.txt",sep = "\t",row.names = FALSE,col.names = TRUE,
              quote = FALSE)
  
  
}


################################# Look for adequate chemo-resistant TNBC cell lines using GR50 data
if (TRUE)
{
  ##### test how concensus is activated in CCLE data
  if (TRUE)
  {
    require(GSVA)
    
    cell.list = c("MDA-MB-231","HCC1806","MDA-MB-468","HCC1143","MDA-MB-453","MDA-MB-436")
    
    # load CCLE expression for TNBC cells
    CCLE.file <- "Data/CCLE/ccle_tnbc_cpm_tmm_log2.txt"
    CCLE.mat <- as.matrix(read.delim(file = CCLE.file,sep = "\t",header = TRUE,stringsAsFactors = FALSE))
    colnames(CCLE.mat) = gsub("BT","BT-",gsub("HS","Hs ",gsub("MDAMB","MDA-MB-",colnames(CCLE.mat))));
    
    # get CCLE GSVA on concensus signatures
    conc.sig = split(fc.summary[[1]],factor(fc.summary$DEG.ID))["UP"]
    names(conc.sig) = paste("CONCENSUS",names(conc.sig),sep = ".")
    conc.gsva = gsva(expr = CCLE.mat,gset.idx.list = conc.sig,method = "zscore")
    
    require(reshape2)
    conc.gsva.df = melt(conc.gsva);colnames(conc.gsva.df) = c("signature","cell.id","score")
  }
  
  ##### test which cell's chemo-resistance is likely driven by activation of consensus signature
  if (TRUE)
  {
    # get GR50
    GR50.data = read.delim(file = "Data/GR50/DS2_datafile.tsv",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
    CCLE.cells = colnames(CCLE.mat)
    GR50.TNBC = subset(GR50.data,Cell.Name %in% CCLE.cells & Small.Molecule.HMS.LINCS.ID.1 %in% c("Taxol","Doxorubicin"))
    
    ### plot out with ggplot2
    # merge CCLE gsva and GR50 data
    plot.data = cbind.data.frame(GR50.TNBC,data.frame(CONCENSUS.UP = conc.gsva[,match(GR50.TNBC$Cell.Name,colnames(conc.gsva))]))
    plot.data$is.available = plot.data$Cell.Name %in% cell.list
    
    # plot
    require(ggplot2)
    require(ggrepel)
    pobj = ggplot(data = plot.data,aes(x = CONCENSUS.UP,y = GRmax)) + 
      geom_point() + facet_grid(. ~ Small.Molecule.HMS.LINCS.ID.1) + 
      geom_hline(yintercept = median(GR50.TNBC$GRmax),colour = "red") + 
      geom_vline(xintercept = c(-3,3),colour = "red") + 
      geom_text_repel(data = subset(plot.data,GRmax > median(GR50.TNBC$GRmax) & is.available & abs(CONCENSUS.UP) > 3),
                      aes(label = Cell.Name),size = 5) +
      scale_color_gradient2(low = "black",mid = "white",high = "red",midpoint = median(GR50.TNBC$GRmax)) + 
      labs(x = "GSVA z-score:Consensus UP",y = "GRmax") +
      theme_bw() + theme(axis.title = element_text(size = 19),axis.text = element_text(size = 15),strip.text = element_text(size = 17))
    
    print(pobj)
  }

}

################################# Plot out EMUDRA results running FDA drugs repurposed by looking at MDA-MB-231 transcriptome perturbations
if (TRUE)
{
  emudra.df = read.delim(file = "Data/EMUDRA/EMUDRA_FDA_repurposed.txt", sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  approved = readLines("Data/EMUDRA/Approved_drugs_list.txt")[-1]
  drug.list = lapply(strsplit(emudra.df[[1]],";|:|,"),function(x) gsub(" ","",x))
  
  emudra.df$is.approved = sapply(drug.list,function(x,y) any(x %in% y),y = approved)
  emudra.df = emudra.df[order(emudra.df$Zscore.MDAMB231),]
  
  emudra.df$short.name = gsub(";(.*)$|:(.*)$|,(.*)$","",emudra.df[[1]])
  
  #### plot out emudra results
  library(ggplot2)
  library(reshape)
  emu.df = subset(emudra.df,is.approved)
  emu.df = emu.df[order(emu.df$NormalizedScore),]
  
  pdata = as.matrix(emu.df[,grep("Zscore",colnames(emu.df))]);
  colnames(pdata) = gsub("Zscore\\.","",colnames(pdata));rownames(pdata) = emu.df$short.name
  pdat = melt(pdata);
  colnames(pdat) = c("drug","cell.line","Zscore")
  
  pobj = ggplot(data = subset(pdat,cell.line == "MDAMB231"),aes(x = drug,y = Zscore,fill = Zscore)) +
    scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) +  
    geom_bar(stat = "identity") + 
    scale_x_discrete(limits = emu.df$short.name[order(emu.df$Zscore.MDAMB231)]) + 
    labs(x = "FDA-approved drugs",y = "EMUDRA Z-score") + 
    geom_hline(yintercept = -3,colour = "red",linetype = "longdash") + 
    guides(fill = FALSE) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45,vjust=1,hjust = 1,size= 15),axis.text.y = element_text(size = 15),
          axis.title = element_text(size = 18),legend.title = element_text(size = 16),legend.text = element_text(size = 14),
          legend.position = "bottom")
  
  print(pobj)
}
