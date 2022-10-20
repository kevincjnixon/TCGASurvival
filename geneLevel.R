##### Set up and combine the raw data #####
tpm<-read.delim("rawData/tcga_Kallisto_tpm.gz")
#tpm from Xena browser are log2(tpm+0.001) - reverse that
rownames(tpm)<-tpm$sample
tpm<-tpm[,-1]
tpm<-(2^tpm)-0.001

tx2gene<-read.delim("rawData/probeMap_gencode.v23.annotation.transcript.probemap")

geneLevel<-function(x, tx2gene){
  # idList<-split(tx2gene$id, tx2gene$gene)
  # splitx<-function(ids,mat){
  #   return(x[which(rownames(mat) %in% ids),] )
  # }
  #
  # no_cores<-parallel::detectCores(logical=T)
  # c1<-parallel::makeCluster(no_cores/4)
  # doParallel::registerDoParallel(c1)
  # message("Exporting objects to cluster...")
  # parallel::clusterExport(c1, list('splitx','x'), envir=environment())
  # message("Making tmpList...")
  # tmpList<-c(parallel::parLapply(c1, idList, fun=splitx, mat=x))
  # message("combinging abundances to gene level...")
  # sums<-c(parallel::parLapply(c1, tmpList, fun=colSums))
  # parallel::stopCluster(c1)
  #
  # return(sums)

  genes<-unique(tx2gene$gene)
  y<-c()
  pb<-txtProgressBar(min=0, max=length(genes), style=3)
  for(i in 1:length(genes)){
    ids<-tx2gene$id[which(tx2gene$gene %in% genes[i])]
    tmp<-x[which(rownames(x) %in% ids),]
    y<-rbind(y, colSums(tmp))
    rownames(y)[i]<-genes[i]
    setTxtProgressBar(pb, i)
  }
  return(y)
}

tmp_gene<-geneLevel(tpm, tx2gene)
saveRDS(tmp_gene, "rawData/TPM_geneLevel.RDS")

