args=commandArgs()

if (length(args)<6) {
  print("Usage: Rscript NodeRemove_noCRv2.R <fileIdentifer> <fileDirectory> [Number of Random files] [Plot1 Position] [Plot2 Position]")
  
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
  plot2Pos = 'topright'
}

evcRank = read.delim(paste(sep='',location,"_evcRankRemove.txt"), header=TRUE)
degRank = read.delim(paste(sep='',location,"_degreeRankRemove.txt"), header=TRUE)
avgMIRank = read.delim(paste(sep='',location,"_avgMIRankRemove.txt"), header=TRUE)
btwnRank = read.delim(paste(sep='',location,"_btwnsRankRemove.txt"), header=TRUE)
clsnRank = read.delim(paste(sep='',location,"_clsnsRankRemove.txt"), header=TRUE)
infRank = read.delim(paste(sep='',location,"_infRankRemove.txt"), header=TRUE)
# chatuRank = read.delim(paste(sep='',location,"_chaturankRemove.txt"), header=TRUE)
randomSuffix = "_randomRemove"

firstPlot = TRUE
pdfPlot = paste(prefix,"_AvgShortestPathPlot.pdf",sep='')
pdf(pdfPlot,width=12,height=8)

for (i in 0:(randomCount-1)){
  randFile = paste(sep='',location,randomSuffix,"_",i,".txt")
  #print(randFile)
  randomRank = read.delim(randFile, header=TRUE)
  if(firstPlot){
    #print('Here')
    plot(randomRank$AvgShortestPath, col='green', xlab = "No. of Nodes Removed", ylab='Avg Shortest Path',ylim=c(1.9,3.5), type='n')
    firstPlot = FALSE
  }
  points(randomRank$AvgShortestPath, col='green',,pch=16,cex=2)
}

points(clsnRank$AvgShortestPath, col='darkorchid3',pch=16,cex=2)
points(degRank$AvgShortestPath, col='hotpink2',pch=16,cex=2)
points(evcRank$AvgShortestPath, col='black',pch=16,cex=2)
points(avgMIRank$AvgShortestPath, col='orange',pch=16,cex=2)
points(btwnRank$AvgShortestPath, col='cyan3',pch=16,cex=2)
# points(chatuRank$AvgShortestPath, col='blue',pch=16,cex=2)
points(infRank$AvgShortestPath, col='red',pch=16,cex=2)


legend(x=plot1Pos, bty = "n", # places a legend at the appropriate place 
       # c('Closeness Centrality','Degree','EVC','WMI', 'Betweenness Centrality', 'INFR','Random'), # puts text in the legend
       c('Closeness Centrality','Degree','EVC','WMI', 'Betweenness Centrality', 'CI','Random'), # changed INFR to Collective Influence
       # c('Closeness Centrality','Degree','EVC','WMI', 'Betweenness Centrality', 'INFR','CR', 'Random'), # original code
       lty=c(1), # gives the legend appropriate symbols (lines)
       lwd=c(5),col=c('darkorchid3','hotpink2','black','orange','cyan3','red','green')) # gives the legend lines the correct color and width
# lwd=c(5),col=c('darkorchid3','hotpink2','black','orange','cyan3','red','blue','green') # original code

dev.off()

firstPlot = TRUE
pdfPlot = paste(prefix,"_LC_Size.pdf",sep='')
pdf(pdfPlot,width=12,height=8)

for (i in 0:(randomCount-1)){
  randFile = paste(sep='',location,randomSuffix,"_",i,".txt")
  #print(randomRank)
  randomRank = read.delim(randFile, header=TRUE)
  if(firstPlot){
    #print('Here')
    #plot(randomRank$GiantComponent, col='green', xlab = "No. of Nodes Removed", ylab='Largest Component Size', ylim=c(0,1800), type='n') #CoSN1a
    plot(randomRank$GiantComponent, col='green', xlab = "No. of Nodes Removed", ylab='Largest Component Size', ylim=c(0,1400), type='n') #AcSN1a
    #plot(randomRank$GiantComponent, col='green', xlab = "No. of Nodes Removed", ylab='Largest Component Size', ylim=c(0,1500), type='n') #ChSN1a
    #plot(randomRank$GiantComponent, col='green', xlab = "No. of Nodes Removed", ylab='Largest Component Size', ylim=c(0,500), type='n') #ChatuHCV1a
    firstPlot = FALSE
  }
  points(randomRank$GiantComponent, col='green',pch=16,cex=2)
}

points(clsnRank$GiantComponent, col='darkorchid3',pch=16,cex=2)
points(degRank$GiantComponent, col='hotpink2',pch=16,cex=2)
points(evcRank$GiantComponent, col='black',pch=16,cex=2)
points(avgMIRank$GiantComponent, col='orange',pch=16,cex=2)
points(btwnRank$GiantComponent, col='cyan3',pch=16,cex=2)
# points(chatuRank$GiantComponent, col='blue',pch=16,cex=2)
points(infRank$GiantComponent, col='red',pch=16,cex=2)

legend(x=plot2Pos, bty = "n",# places a legend at the appropriate place 
       # c('Closeness Centrality','Degree','EVC','WMI', 'Betweenness Centrality', 'INFR','Random'), # puts text in the legend
       c('Closeness Centrality','Degree','EVC','WMI', 'Betweenness Centrality', 'CI','Random'), # changed INFR to Collective Influence
       # c('Closeness Centrality','Degree','EVC','WMI', 'Betweenness Centrality', 'INFR','CR','Random'), # original code
       lty=c(1), # gives the legend appropriate symbols (lines)
       lwd=c(5),col=c('darkorchid3','hotpink2','black','orange','cyan3','red','green')) # gives the legend lines the correct color and width
# lwd=c(5),col=c('darkorchid3','hotpink2','black','orange','cyan3','red','blue','green')) # original code

dev.off()
