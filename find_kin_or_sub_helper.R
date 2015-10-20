find_kins <- function(query){
  nodes <- read.csv("/Users/aolow/Dropbox/MS1_Rcode/NEW_NODES3.csv", as.is=T)
  kinases <- c()
  substrates <- c()
  neither <-c()
  both <- c()
  for(i in 1: length(query)){
    if(query[i] %in% nodes$Kinase & query[i] %in% nodes$Substrate)
    {both <- c(both, query[i])
     print(paste0(query[i], " is both a kinase and a substrate."))}
    else if(query[i] %in% nodes$Kinase)
    {kinases <- c(kinases, query[i])
     print(paste0(query[i], " is a kinase."))}
    else if(query[i] %in% nodes$Substrate)
    {substrates <- c(substrates, query[i])
     print(paste0(query[i], " is a substrate."))}
    else{
      neither <- c(neither, query[i])
      print(paste0(query[i], " wasn't found."))}
  }
}

check <- read.table("/Users/aolow/Desktop/query.txt", sep="\t")
check <- check$V1