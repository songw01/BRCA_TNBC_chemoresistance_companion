load.data <- function (filename, gsub.from = "\\.", gsub.to = "-") 
{
    df <- read.delim(file = filename, sep = "\t", header = T)
    out <- as.matrix(df[, 2:ncol(df)])
    rownames(out) <- as.character(df[[1]])
    if (!is.null(gsub.from)) 
        colnames(out) <- gsub(gsub.from, gsub.to, colnames(out))
    return(out)
}
############################# FET pipeline
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

make.Pairwise.Tables <- function(geneSets1,geneSets2,background)
{
 mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
 mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

 d <- t(mem1) %*% mem2;
 b <- abs(t(mem1) %*% (mem2-1))
 c <- abs(t(mem1-1) %*% (mem2))
 a <- t(mem1-1) %*% (mem2-1);

 ij <- do.call(rbind,lapply(1:length(geneSets1),function(i,j) cbind(rep(i,length(j)),j),j = 1:length(geneSets2)))

 pairwise.tables <- lapply(1:nrow(ij),function(i,ij,a,b,c,d) as.table(matrix(c(a[ij[i,1],ij[i,2]],b[ij[i,1],ij[i,2]],c[ij[i,1],ij[i,2]],d[ij[i,1],ij[i,2]]),nrow = 2)),ij = ij,a = a,b = b,c = c,d = d)
 names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
 return(pairwise.tables)
}

do.FisherExactTest <- function(table.count,N_bg = NULL)
{
 if (is.null(N_bg)) N_bg = sum(rowSums(table.count))
 
 out <- fisher.test(x = table.count,or = 1,alternative = "greater")
 odds.ratio <- out$estimate
 p.value <- out$p.value;
 geneSet1.count <- rowSums(table.count)[2]
 geneSet2.count <- colSums(table.count)[2]
 expected.count <- geneSet1.count/N_bg * geneSet2.count
 overlap.count <- table.count[2,2];
 fold.change <- overlap.count/expected.count
 
 out <- c(N_bg,geneSet1.count,geneSet2.count,expected.count,overlap.count,fold.change,odds.ratio,p.value)
 names(out) <- c("Background","set1_size","set2_size","expected.overlap","actual.overlap","enrichment.foldchange","odds.ratio","FET_pvalue")
 return(out)
}

require(parallel)
require(foreach)
require(iterators)
perform.AllPairs.FET <- function(geneSets1,geneSets2,background,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
{
 pairwise.tables <- make.Pairwise.Tables(geneSets1,geneSets2,background)
 
 if (do.multicore)
 {
  #registerDoMC(n.cores)
  set.parallel.backend(n.cores)
  fact <- do.call(c,lapply(1:n.cores,function(i,dn) rep(i,dn),dn = ceiling(length(pairwise.tables)/n.cores)))
  fact <- factor(fact[1:length(pairwise.tables)])
  split.tables <- lapply(split(1:length(pairwise.tables),fact),function(i,obj) obj[i],obj = pairwise.tables);rm(fact)
  output <- foreach (tbl = split.tables,.combine = 'c') %dopar% {
            out <- lapply(tbl,do.FisherExactTest) 
			return(out)
  }
 
 }else{
  output <- lapply(pairwise.tables,do.FisherExactTest)
 }
 # collect outputs
 output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
 int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
 int.element <- do.call(c,int.element)
 if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
 
 return(output)
}

make.paired.Tables <- function(geneSets1,geneSets2,ij,background)
{
 
 #mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
 #mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

 pairwise.tables <- mapply(FUN = function(s1,s2,z) {
                                         v1 <- rep(0,length(z));v1[which(z %in% s1)] <- 1; 
										 v2 <- rep(0,length(z));v2[which(z %in% s2)] <- 1;
										 # n11,n12,n21,n22
										 as.table(matrix(c(sum(abs(v1 - 1) * abs(v2 - 1)),sum(abs(v2 - 1) * v1),sum(abs(v1 - 1) * v2),sum(v1 * v2)),nrow = 2))
                                },s1 = geneSets1[ij[,1]],s2 = geneSets2[ij[,2]],MoreArgs = list(z = background),SIMPLIFY = FALSE)
 names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
 return(pairwise.tables)
}


perform.ijPairs.FET <- function(geneSets1,geneSets2,ij,background,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
{
 require(doParallel)
 pairwise.tables <- make.paired.Tables(geneSets1,geneSets2,ij,background)
 
 if (do.multicore & getDoParWorkers() == 1)
 {
  export.func <- c("do.FisherExactTest","make.Pairwise.Tables","make.paired.Tables")
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  cat(paste("registered:",getDoParWorkers()," cores\n",sep = ""))
  fact <- do.call(c,lapply(1:n.cores,function(i,dn) rep(i,dn),dn = ceiling(length(pairwise.tables)/n.cores)))
  fact <- factor(fact[1:length(pairwise.tables)])
  split.tables <- lapply(split(1:length(pairwise.tables),fact),function(i,obj) obj[i],obj = pairwise.tables);rm(fact)
  cat("commence parallelized enrichment analyses\n")
  output <- foreach (tbl = split.tables,.combine = 'c',.export = export.func) %dopar% {
            out <- lapply(tbl,do.FisherExactTest) 
			return(out)
  }
  
  stopCluster(cl)
 
 }else{
  output <- lapply(pairwise.tables,do.FisherExactTest)
 }
 # collect outputs
 output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
 int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
 int.element <- do.call(c,int.element)
 if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
 
 output <- output[order(output$FET_pvalue),]
 return(output)
}

##### master functions

##### master functions
load.modules <- function(module.file,file.type)
{
 if (file.type == "MEGENA" | file.type == "GMT")
 {
  output <- read.geneSet(module.file)
  if (file.type == "GMT") output <- lapply(output,function(x) x[-1])
 }
 
 if (file.type == "WGCNA" | file.type == "TABLE")
 {
  module.table <- read.delim(file = module.file,sep = "\t",header = T)
  #output <- lapply(split(as.character(module.table$Gene.symbol),factor(module.table$module)),unique)
  output <- lapply(split(as.character(module.table[[1]]),factor(module.table[[ncol(module.table)]])),unique)
 }
 
 output <- lapply(output,function(x) x[!is.na(x) & x != "NA"])
 
 cat("module sizes...\n")
 print(sapply(output,length))
 return(output)
}

perform.FET.output <- function(geneSets1,geneSets2,background,outputfname,adjust.FET.pvalue = T,pvalue.cutoff = 0.05,do.multicore = F,n.cores = NULL)
{
 if (is.null(names(geneSets1))) 
 {
  cat("#### No names were tagged to geneSets1 #####\nProviding names to geneSets1...\n")
  names(geneSets1) <- paste("A",1:length(geneSets1),sep = "")
 }
 
 if (is.null(names(geneSets2))) 
 {
  cat("#### No names were tagged to geneSets2 #####\nProviding names to geneSets1...\n")
  names(geneSets2) <- paste("B",1:length(geneSets2),sep = "")
 }
 cat("###### Performing Fisher Exact Test for over-representation\n")
 cat(paste("- Number of geneSets1:",length(geneSets1),"\n","- Number of geneSets2:",length(geneSets2),"\n",sep = ""))
 
 FET.table <- perform.AllPairs.FET(geneSets1 = geneSets1,geneSets2 = geneSets2,background = background,adjust.FET.pvalue = adjust.FET.pvalue,do.multicore = do.multicore,n.cores = n.cores)
 
 # output gene sets
 cat("###### Outputting files...\n- Output gene sets\n")
 geneSet.files <- paste(outputfname,"_geneSets",c(1,2),".txt",sep = "")
 cat(paste(geneSet.files,"\n",sep = ""))
 output.status <- output.geneSet.file(geneSet = geneSets1,outputfname = geneSet.files[1])
 output.status <- output.geneSet.file(geneSet = geneSets2,outputfname = geneSet.files[2])
 cat("- Output FET output table\n")
 write.table(FET.table,file = paste(outputfname,"_FET-Table.txt",sep = ""),sep = "\t",row.names = F,col.names = T,quote = F)
 cat("- Output significant results\n")
 sig.table <- FET.table[FET.table$corrected.FET.pvalue < pvalue.cutoff,];
 write.table(sig.table,file = paste(outputfname,"_SignificantFET-Table.txt",sep = ""),sep = "\t",row.names = F,col.names = T,quote = F)
 
 return(0) 
}

