args=commandArgs()

if (length(args)<6) {
  print("Usage: Rscript NodeRemove.R <fileIdentifer> <fileDirectory> [Number of Random files] [Plot1 Position] [Plot2 Position]")
  
  print(paste("Plot1 and Plot2 Position can be 'topleft', 'topright','bottomleft' or 'bottomright' as an example.",
  "If none are entered, defaults are 'bottomright' for Plot1 and 'topleft' for Plot2.",
  "Check R manual 'Legend' function for full description."))
  
  stop("Insufficient input parameters. Requires both file prefix and directory where ranking files are stored.\n", call.=FALSE)
} 

prefix = args[6]
fileDir = args[7]
location =paste(fileDir,prefix,sep="/")

if(length(args) > 7){
  randomCount = as.numeric(args[8])
}

if(length(args) > 8){
  plot1Pos = args[9]
}else{
  plot1Pos = 'bottomright'
}

if(length(args) > 9){
  plot2Pos = args[10]
}else{
  plot2Pos = 'topleft'
}

evcRank = read.delim(paste(sep='',location,"_evcRankRemove.txt"), header=TRUE)
degRank = read.delim(paste(sep='',location,"_degreeRankRemove.txt"), header=TRUE)
avgMIRank = read.delim(paste(sep='',location,"_avgMIRankRemove.txt"), header=TRUE)
btwnRank = read.delim(paste(sep='',location,"_btwnsRankRemove.txt"), header=TRUE)
clsnRank = read.delim(paste(sep='',location,"_clsnsRankRemove.txt"), header=TRUE)
infRank = read.delim(paste(sep='',location,"_infRankRemove.txt"), header=TRUE)
chatuRank = read.delim(paste(sep='',location,"_chaturankRemove.txt"), header=TRUE)
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
    plot(randomRank$AvgShortestPath, col='green', xlab = "No. of Nodes Removed", ylab='Avg Shortest Path',ylim=c(1.9,2.8), type='n')
    firstPlot = FALSE
  }
  points(randomRank$AvgShortestPath, col='green')
}

points(clsnRank$AvgShortestPath, col='darkorchid3')
points(degRank$AvgShortestPath, col='hotpink2')
points(evcRank$AvgShortestPath, col='black')
points(avgMIRank$AvgShortestPath, col='orange')
points(btwnRank$AvgShortestPath, col='cyan3')
points(chatuRank$AvgShortestPath, col='blue')
points(infRank$AvgShortestPath, col='red')


legend(x=plot1Pos, bty = "n", # places a legend at the appropriate place 
       c('Closeness Centrality','Degree','EVC','NAMI', 'Betweenness Centrality', 'INFR','CR', 'Random'), # puts text in the legend
       lty=c(1), # gives the legend appropriate symbols (lines)
       lwd=c(5),col=c('darkorchid3','hotpink2','black','orange','cyan3','red','blue','green')) # gives the legend lines the correct color and width

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
    plot(randomRank$Diameter, col='green', xlab = "No. of Nodes Removed", ylab='Network Diameter', ylim=c(0,25), type='n')
    firstPlot = FALSE
  }
  points(randomRank$Diameter, col='green')
}

points(clsnRank$Diameter, col='darkorchid3')
points(degRank$Diameter, col='hotpink2')
points(evcRank$Diameter, col='black')
points(avgMIRank$Diameter, col='orange')
points(btwnRank$Diameter, col='cyan3')
points(chatuRank$Diameter, col='blue')
points(infRank$Diameter, col='red')

legend(x=plot2Pos, bty = "n",# places a legend at the appropriate place 
       c('Closeness Centrality','Degree','EVC','NAMI', 'Betweenness Centrality', 'INFR','CR','Random'), # puts text in the legend
       lty=c(1), # gives the legend appropriate symbols (lines)
       lwd=c(5),col=c('darkorchid3','hotpink2','black','orange','cyan3','red','blue','green')) # gives the legend lines the correct color and width

dev.off()
