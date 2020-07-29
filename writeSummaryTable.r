library(RSQLite)
library(reshape2)
library("vegan")
library(tidyverse)
library(igraph)
source("networkCalculations.r")
args <- commandArgs(trailingOnly = TRUE)
print(args)
wd <- args[1]
prefix<-args[2]
num<-as.integer(args[3])
resultsFolder <-args[4]
scNum<-as.integer(args[5])
r<-as.integer(args[6])

pairwiseDiv<-function(mat){
  newmat<-tcrossprod(mat>0)
  newmat<-newmat/rowSums(mat>0)
  return(1-newmat)  
}

fetchdb<-function(dbname,query,numQuery = 20000000) {
  r<-dbSendQuery(conn=dbname, query)
  er<-dbFetch(r,numQuery)
  while(!dbHasCompleted(r)){
    er <- rbind(er, dbFetch(r, numQuery))
    print(nrow(er))
  }
  dbClearResult(r)
  return(er)
}

buildNetFromAdj<-function(mat, cutoff = NULL) {
  g<-graph.adjacency(mat,weighted=TRUE,mode="directed",diag=F)
  if(is.null(cutoff)){
  cutoff<-quantile(E(g)$weight,0.99)
  }
  g<-delete.edges(g, which(E(g)$weight < cutoff))
  return(g)
}

calStat<-function(prefix = "varMigTest", r, num , s,
                  samplingPeriod = 30, n_host = 10000, cutoffTime = 3600){
	if (s == 0 || s == 4){
	sampleSqlFile <- paste(wd, prefix,"_",num, "_s", s, "_sd.sqlite",sep="")	
	}else{
	sampleSqlFile <- paste(wd, prefix,"_",num, "_s", s,"_r",r, "_sd.sqlite",sep="")		
	}
  #only calculate if the file exists
  
  if (file.exists(sampleSqlFile)){
  db<-dbConnect(SQLite(),dbname = sampleSqlFile)
  
  selMode<-"S"
  if (s>3){
  	selMode<-"G"
  }
  
  if ((s+1)%%4 == 1) {
  	IRS = 0
  }else if ((s+1)%%4 == 2) {
  	IRS = 2
  }else if ((s+1)%%4 == 3) {
  	IRS = 5
  }else{
  	IRS = 10
  }
  
  sc<-"select id, source from sampled_genes"
  geneInfo<-fetchdb(db, sc)
  colnames(geneInfo)[1]<-"gene_id"
  
  sc<-"select id, ind, gene_id from sampled_strains"
  strainInfo<-fetchdb(db, sc)
  colnames(strainInfo)[1]<-"strain_id"
  
  sc<-"select gene_id, locus, allele from sampled_alleles"
  alleleInfo<-fetchdb(db, sc)
  geneAllele<-dcast(alleleInfo, gene_id~locus, value.var="allele")
  colnames(geneAllele)[2:ncol(geneAllele)]<-paste("l",1:(ncol(geneAllele)-1),sep="")
  #only select infections that are expressing gene ids
  
  
  #get all summary info
  sc<-paste("select * from summary  WHERE time >",cutoffTime, sep="")                                                    
  summaryInfo<-fetchdb(db, sc)
  summaryInfo<-summaryInfo %>% mutate(num = num, selMode = selMode, IRS = IRS, r = r)

  summaryTable<-summaryInfo%>%mutate(EIR = n_infected_bites/n_host/samplingPeriod*360,    
                                     Prevalence = n_infected/n_host, 
                                     MOI = n_infections/n_infected)

  sc<-paste("select * from pool_size WHERE time >",cutoffTime, sep="")   
  summaryPoolSize<-fetchdb(db, sc)

  summaryTable<-summaryTable%>%left_join(summaryPoolSize, by = "time")

  sc<-paste("select time, locus, n_circulating_alleles from summary_alleles WHERE time >",cutoffTime, sep="")                                                    
  summaryAlleles<-fetchdb(db, sc)
  
  summaryTable<-summaryTable%>%left_join(summaryAlleles%>%filter(locus==0)%>%select(time, n_circulating_alleles), by = "time")
  summaryTable<-summaryTable%>%left_join(summaryAlleles%>%filter(locus==1)%>%select(time, n_circulating_alleles), by = "time")
	
  #get strain info
  sc<-paste("select time, strain_id from sampled_infections where gene_id > -1 and time >"  ,cutoffTime,sep="")
  sampledInf<-fetchdb(db,sc)
  infStrain<-sampledInf %>% mutate(uniqStrain = 1:nrow(sampledInf)) 
  infStrain<-left_join(infStrain, strainInfo, by="strain_id")
  infStrain<-left_join(infStrain, geneInfo, by="gene_id")
  infStrain<-left_join(infStrain, geneAllele, by="gene_id")
  infStrain<-infStrain %>% mutate(num = num, selMode = selMode, IRS = IRS, r = r)
  
  
  #calculate Niche Divergence                              
  nicheDiv <- data.frame(time = seq(min(infStrain$time),max(infStrain$time),samplingPeriod), nicheDiv = 0,
  alleleShannon=0, alleleSimpson = 0, geneShannon=0, geneSimpson= 0,
  f_01_averageLocalClusteringCoeff=0,f_02_averageLocalClusteringCoeffWeighted=0,
f_03_globalClusteringCoeff=0,f_04_gdensity=0,
f_05_proportionSingletons=0,f_06_proportionEndpoints=0,
f_07_meanDegree=0,f_08_assortativityDegree=0,f_09_meanStrength=0,f_10_straightness=0,
f_11_entropyDegreeDistribution=0,f_12_ratioComponents=0,f_13_giantConnectedRatio=0,
f_14_evennessComponentSize=0,f_15_CentralPointDominance=0,f_16_meanEccentricity=0,
f_17_gdiameter=0,f_18_meanDiameterComponents=0,f_19_globalEfficiency=0,
f_20_averageClosenessCentrality=0,f_21_mot1=0,f_22_mot2=0,f_23_mot3=0,f_24_mot4=0,
f_25_mot5=0,f_26_mot6=0,f_27_mot7=0,f_28_mot8=0,f_29_mot9=0,f_30_mot10=0,f_31_mot11=0,
f_32_mot12=0,f_33_reciprocity=0,f_34_inoutCorrelation=0,f_35_communitySizeEvenness=0,
f_36_numberCommonCommunities=0,f_37_CommunityRatio=0,f_38_modularity=0,f_39_meanFST=0,
f_40_maxFST=0,f_41_minFST=0,f_42_giniFST=0
  )

  #head(nicheDiv)
  for (i in 1:length(unique(infStrain$time))) {
    ## first get all the strain info per layer
    tempInf <- infStrain %>% filter(time == unique(infStrain$time)[i])
    paral1<-acast(tempInf, uniqStrain~l1, length, value.var = "uniqStrain")	
    paral2<-acast(tempInf, uniqStrain~l2, length, value.var = "uniqStrain")	
	paraCb<-cbind(paral1, paral2)
    outMat<-pairwiseDiv(as.matrix(paraCb))
    nicheDiv[i,2] <- mean(c(outMat[lower.tri(outMat)],outMat[upper.tri(outMat)]))
    #diversity measures
    nicheDiv[i,3]<-tempInf%>%count(l1)%>%select(n)%>%vegan::diversity()
    nicheDiv[i,4]<-tempInf%>%count(l1)%>%select(n)%>%vegan::diversity("simpson")
    nicheDiv[i,5]<-tempInf%>%count(gene_id)%>%select(n)%>%vegan::diversity()
    nicheDiv[i,6]<-tempInf%>%count(gene_id)%>%select(n)%>%vegan::diversity("simpson")
    #print(nicheDiv[i,3])
    #print(nicheDiv[i,4])

    #other network property calculations
    
    ng<-buildNetFromAdj((1-outMat))
    ntp<-calculateFeatures(ng, cutoff=0.95)
    nicheDiv[i,7:(length(ntp)+6)]<-ntp
   }

  summaryTable<-summaryTable%>%left_join(nicheDiv, by = "time")

  #write.table(infStrain, file = paste(resultsFolder,prefix, "_", num,"_sampleInfectionTable.txt",sep=""),
  #quote=F,sep="\t",row.names=F,col.names=T,append = T)

  write.table(summaryTable, file = paste(resultsFolder,prefix, "_", num,"_summaryTable.txt",sep=""),
  quote=F,sep="\t",row.names=F,col.names=F,append = T)

  dbDisconnect(db) 
                                     
  
  
  
  
  }
}



calStat(prefix = prefix, r=r, num = num , s = scNum)
