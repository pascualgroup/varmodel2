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
	
  write.table(summaryTable, file = paste(resultsFolder,prefix, "_", num,"_summaryTableFast.txt",sep=""),
  quote=F,sep="\t",row.names=F,col.names=F,append = T)

  dbDisconnect(db) 
                                     
  
  
  
  
  }
}




calStat(prefix = prefix, r=r, num = num , s = scNum)
