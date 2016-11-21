args=commandArgs()

if (length(args)<6) {
  print("Usage: Rscript NodeRemove.R <fileIdentifer> <fileDirectory> [Number of Random files]")
  stop("Insufficient input parameters. Requires both file prefix and directory where ranking files are stored.\n", call.=FALSE)
} 

prefix = args[6]
fileDir = args[7]
location =paste(fileDir,prefix,sep="/")

if(length(args) > 6){
  randomCount = as.numeric(args[8])
}

evcRank = read.delim(paste(sep='',location,"_evcRankRemove.txt"), header=TRUE)
degRank = read.delim(paste(sep='',location,"_degreeRankRemove.txt"), header=TRUE)
avgMIRank = read.delim(paste(sep='',location,"_avgMIRankRemove.txt"), header=TRUE)
btwnRank = read.delim(paste(sep='',location,"_btwnsRankRemove.txt"), header=TRUE)
clsnRank = read.delim(paste(sep='',location,"_clsnsRankRemove.txt"), header=TRUE)
infRank = read.delim(paste(sep='',location,"_infRankRemove.txt"), header=TRUE)
randomSuffix = "_randomRemove"

firstPlot = TRUE
pdfPlot = paste(prefix,"_AvgShortestPathPlot.pdf",sep='')
pdf(pdfPlot,width=12,height=8)

for (i in 0:(randomCount-1)){
  randFile = paste(sep='',location,randomSuffix,"_",i,".txt")
  #print(randomRank)
  randomRank = read.delim(randFile, header=TRUE)
  if(firstPlot){
    #print('Here')
    plot(randomRank$AvgShortestPath, col='green', xlab = "No. of Nodes Removed", ylab='Avg Shortest Path (weighted)', ylim=c(1.9,2.8), main = "Breaking the viral mutation network", type='n')
    firstPlot = FALSE
  }
  points(randomRank$AvgShortestPath, col='green')
}

points(clsnRank$AvgShortestPath, col='darkorchid3')
points(degRank$AvgShortestPath, col='hotpink2')
points(evcRank$AvgShortestPath, col='black')
points(avgMIRank$AvgShortestPath, col='orange')
points(btwnRank$AvgShortestPath, col='cyan3')
points(infRank$AvgShortestPath, col='red')

dev.off()

firstPlot = TRUE
pdfPlot = paste(prefix,"_DiameterPlot.pdf",sep='')
pdf(pdfPlot,width=12,height=8)

for (i in 0:(randomCount-1)){
  randFile = paste(sep='',location,randomSuffix,"_",i,".txt")
  #print(randomRank)
  randomRank = read.delim(randFile, header=TRUE)
  if(firstPlot){
    #print('Here')
    plot(randomRank$Diameter, col='green', xlab = "No. of Nodes Removed", ylab='Network Diameter (weighted)', ylim=c(0,20), main = "Breaking the viral mutation network", type='n')
    firstPlot = FALSE
  }
  points(randomRank$Diameter, col='green')
}

points(clsnRank$Diameter, col='darkorchid3')
points(degRank$Diameter, col='hotpink2')
points(evcRank$Diameter, col='black')
points(avgMIRank$Diameter, col='orange')
points(btwnRank$Diameter, col='cyan3')
points(infRank$Diameter, col='red')

dev.off()
