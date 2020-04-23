# Author: PDawgz (~w~)v
# These function will have to be edited accordingly to a user's
# dataset, since weight type and number of weights will differ.

# Import this module inside python console for more customisable
# analysis.

import re, os, argparse, sys
import igraph as ig
import csv, random
from datetime import datetime as dt
import numpy as np
import copy as cp

######################### parse edge data as network #########################
# Parses the file which describes the edges into igraph.
# Typically edges have information about the 2 nodes it is connected to.
# This function does not parse in information about edge direction.
# Refer to iGraph manual for python if directed edges are required.

def parseData(fileName, delim):
    skipHeader = False
    pos_nodeID_dict = {}
    nodeID_pos_dict = {}
    pos_edgeID_dict = {}
    g = ig.Graph()    
    with open(fileName, 'rb') as csvfile:        
        fileHandle = csv.reader(csvfile, delimiter= delim)
        edgeID = 0    
        nodeID = 0        
        #Retrieve node (pos1,pos2) data and
        #edge data(pvalue, mi_weight, one_minus_mi)
        #Customise this part if needed.         
        for row in fileHandle:
            if(not skipHeader):
                skipHeader = True
                continue
            pos1 = row[0]
            pos2 = row[1]
            mi_weight = row[2]
            pValue = row[3]
            one_minus_mi = row[4]            
            # Dealing with 2 nodes, igraph nodeID is sequential
            # starts from 0 to n.
            if(pos1 not in pos_nodeID_dict.keys()):
                pos_nodeID_dict[pos1] = nodeID                
                nodeID_pos_dict[nodeID] = pos1                
                g.add_vertices(1)
                g.vs[nodeID]['Position'] = pos1
                nodeID += 1
            if(pos2 not in pos_nodeID_dict.keys()):
                pos_nodeID_dict[pos2] = nodeID
                nodeID_pos_dict[nodeID] = pos2                
                g.add_vertices(1)
                g.vs[nodeID]['Position'] = pos2
                nodeID += 1            
            #Dealing with the edge between the 2 nodes            
            if((pos1,pos2) not in pos_edgeID_dict.keys()):
                pos_edgeID_dict[(pos1,pos2)] = edgeID                
            #Stored the id of nodes of positions and
            #stored the id of edges between positions.
            #Time to add this edge into the graph.            
            g.add_edges([(pos_nodeID_dict[pos1],pos_nodeID_dict[pos2])])            
            #Add some attributes to the edge
            #Customise this part if needed.
            g.es[edgeID]['MI'] = float(mi_weight)
            g.es[edgeID]['PValue'] = float(pValue)
            g.es[edgeID]['1-MI'] = float(one_minus_mi)
            edgeID += 1
    return g



######################### calculates avg MI per node #########################

def node_avgMI(nodeID, network):
    
    miSum = (sum(network.es.select(_source=nodeID)['MI']) + 
             sum(network.es.select(_target=nodeID)['MI'])) 
    
    avgMI = miSum / network.vs[nodeID].degree()
    return avgMI

######################### calculates avg MI for all nodes in graph #########################

def allNode_avgMI(network):
    outputList = []
    nodeIDlist = [v.index for v in network.vs]
    
    for i in nodeIDlist:
        outputList.append(node_avgMI(i, network))
    return outputList

######################### output network statistics #########################
# EVC -> EigenVectorCentrality
# Columns are:
#   Col 01: Position
#   Col 02: NodeAvgMI
#   Col 03: Degree
#   Col 04: EVC Score (MI weighted)
#   Col 05: InfRank
#   Col 06: Betweenness (1-MI weighted)
#   Col 07: Closeness (1-MI weighted)
#   Col 08: Cluster Number
#   Col 09: ChatuRank

def nodeStats2(nw, nodeStatDir,deRefINFRank, delim):
    with open(nodeStatDir,'wb') as wHandle:
        evc_score = nw.eigenvector_centrality(directed=False, weights=nw.es['MI'])
        between3 = nw.betweenness(weights=nw.es['1-MI'],directed=False)
        closeness3 = nw.closeness(weights=nw.es['1-MI'])
        nwClusters = nw.community_leading_eigenvector(weights=nw.es['MI'])
        numClusters = len(nwClusters)
        nodeRank = chatuRankAll(nw, 'MI', 50)
        infDict = infRankMatch(deRefINFRank, nw, delim)
                
        wHandle.write('Position\tNodeAvgMI\tDegree\tEVC-Score(MI)\t')
        wHandle.write('InfRank\tBetweenness(1-MI)\t')
        wHandle.write('Closeness(1-MI)\tCluster\tRank\n')        
        for index in range(0,len(nw.vs)):

            currPos = nw.vs[index]['Position']
            
            for clusterID in range(0, numClusters):
                subG = nwClusters.subgraph(clusterID)
                vertSeqObj = subG.vs.select(Position_eq = currPos)
                if(len(vertSeqObj) == 1):
                    wHandle.write('{}\t{}\t'.format(currPos, node_avgMI(index,nw)))
                    wHandle.write('{}\t{}\t'.format(nw.vs[index].degree(),evc_score[index]))
                    if(infDict[int(currPos)] == -1):
                        wHandle.write('NA\t')
                    else:
                        wHandle.write('{}\t'.format(infDict[int(currPos)]))
                    wHandle.write('{}\t'.format(between3[index]))
                    wHandle.write('{}\t'.format(closeness3[index]))
                    wHandle.write('{}\t{}\n'.format(clusterID,nodeRank[currPos]))
                else:
                    # currPos is not in this cluster, check next cluster 
                    pass



######################### matching influential rank to nodes #########################
# some nodes are not ranked, therefore 'NA' symbol replaced blank ranking
# Uses the infRank output from tool from the Morone paper
# note: writes to file. So everything becomes a str type

def infRankMatch(deRefINFRank, network, delim):
    infRankDict = {}
    with open(deRefINFRank, 'rb') as fileHandle:
        readFile = csv.reader(fileHandle, delimiter= delim)
        for line in readFile:        
            infRankDict[int(line[0])] = int(line[1])
    for node in [network.vs[v.index]['Position'] for v in network.vs]:
        key = node
        if(int(key) not in infRankDict.keys()):
            infRankDict[int(key)] = -1
    return infRankDict #dictionary with int(aapos) -> infRank
      
######################### remove NA and retrive aa pos using infRank #########################
# Using the AA Position <-> InfRank dictionary, we transform this to nodeID<->InfRank 
# Dictionary, so we can remove nodes from network using nodeID, which is sorted
# by infRank(Used for noderemove test).

def transformInfRank(network, infRankDict):
    NAcount = 0
    nodeID2InfRank = {}
    for key in infRankDict.keys(): #links nodeID -> inf rank
        #print key
        nodeID = network.vs.select(Position_eq = str(key))[0].index
        nodeID2InfRank[nodeID] = infRankDict[key]
    #now we have a list of infRank, ordered by nodeID    
    infRankList = [nodeID2InfRank[v] for v in sorted(nodeID2InfRank.keys())]
    #list of nodeID sorted by infRank
    sortedIndex = np.array(infRankList).argsort() 
    tempSorted = cp.copy(list(sortedIndex))
    #getting rid of those nodeIDs without a rank -> anything -1
    for index in sortedIndex:
        if (infRankList[index] == -1):
            tempSorted.remove(index)
    aapos = [network.vs[v]['Position'] for v in tempSorted]
    
    print "AA position list length = {}".format(len(aapos))
    return aapos    
    
######################### Update Immune targeted sites #########################

def handleMutSites(fileA, fileB, outFile):
    aHash = {}
    with open(fileA, 'rb') as csvfile:
        readfile = csv.reader(csvfile)
        for line in readfile:
            aHash[int(line[0])] = line[1:3]
    
    with open(fileB, 'rb') as csvfile:
        readfile = csv.reader(csvfile)
        with open(outFile, 'wb') as writeHandle:
            for line in readfile:
                pos = line[0]
                if(int(pos) in aHash.keys()):
                    writeHandle.write('{}\t{}\t{}\n'.format(int(pos), aHash[int(pos)][0], aHash[int(pos)][1]))
                else:
                    writeHandle.write('{}\t0\t0\n'.format(int(pos)))

######################### randomise network #########################
# This function is to take a real network and randomly redraw edges.
# The aim is to break up the "real" partners and see if the nodes
# still connect the same way (e.g. form similar clusters). Returns 
# a simplified network with one node attribute and weight attribute
# specified by user.

def randomiseEdge(network, weightAttr, nodeAttr):    
    weightDist = network.es[weightAttr]
    numEdge = len(weightDist)
    nodeList = network.vs[nodeAttr] #contains attribute info of nodes
    numNode = len(nodeList)
    numPotentialEdge = (numNode*(numNode-1))/2
    edgeChance = float(numEdge)/numPotentialEdge    
    randG = network.copy()
    randG.delete_edges(range(0,numEdge))    
    #Adding edge to graph randG randomly
    for i in range(0,numNode):
        j = i+1
        while(j < numNode):
            random.seed(dt.now().time())
            if(random.uniform(0,1) <= edgeChance):
                #add edge between i and j
                randG.add_edges([(i,j)])
                edgeID = randG.get_eid(i,j)
                # randomly select a weight and add it
                # to the edge
                randWeight = weightDist[random.randint(0,numEdge-1)]
                randG.es[edgeID][weightAttr] = randWeight
            j+=1
    return randG

######################### removes node to check stats #########################
# removes node using a list sorted by user
# idList = [v.index for v in nw.vs.select()]

def nodeRemove(network, rankAttr, outFile, weightAttr, isInfRank = False):

    nw = network.copy()    
    aaPosList = []
        
    if(isInfRank):
        aaPosList = transformInfRank(nw, rankAttr)
    else:
        sortedIndex = np.array(rankAttr).argsort()[::-1]
        for i in sortedIndex:
            aaPosList.append(nw.vs[i]['Position'])
    print "Starting to remove nodes and write output to:"
    print "{}".format(outFile)
    
    with open(outFile, 'wb') as writeFile:
        # writeFile.write('Pos\tAvgPathLen\tAvgShortestPath\tDiameter\tRemainingNodes\n')
        writeFile.write('Pos\tAvgPathLen\tAvgShortestPath\tDiameter\tGiantComponent\tRemainingNodes\n')
        printCounter = 0    
        for pos in aaPosList:
            nodeID = [v.index for v in nw.vs.select(Position_eq=pos)][0]
            nw.delete_vertices(nodeID)
            if(printCounter%10 == 0 and printCounter > 0):
                print "Removed {} out of {} nodes.".format(printCounter,len(aaPosList))
            if(len(nw.vs) >= 1):
                writeFile.write('{}\t{}\t'.format(pos,nw.average_path_length()))
                writeFile.write('{}\t{}\t'.format(avgShortestPath(nw,weightAttr),nw.diameter(weights=nw.es[weightAttr])))
                writeFile.write('{}\t'.format(len(nw.components().giant().vs)))
                writeFile.write('{}\n'.format(len(nw.vs)))
            else:
                break
            printCounter += 1
    print "Done."

######################### removes node to check stats 2 #########################
# removes nodes randomly.

def randNodeRemove(network, outFile, weightAttr):
    
    nw = network.copy()
    aaPosList = []
        
    for i in [v.index for v in nw.vs]:
        aaPosList.append(nw.vs[i]['Position'])
    
    random.seed(dt.now().time())
    random.shuffle(aaPosList)
    
    print "Starting to remove nodes and write output to:"
    print "{}".format(outFile)
    
    with open(outFile, 'wb') as writeFile:
        writeFile.write('Pos\tAvgPathLen\tAvgShortestPath\tDiameter\tGiantComponent\tRemainingNodes\n')
        printCounter = 0    
        for pos in aaPosList:
            nodeID = [v.index for v in nw.vs.select(Position_eq=pos)][0]
            nw.delete_vertices(nodeID)
            if(printCounter%10 == 0 and printCounter > 0):
                print "Removed {} out of {} nodes.".format(printCounter,len(aaPosList))
            if(len(nw.vs) >= 1):
                writeFile.write('{}\t{}\t'.format(pos,nw.average_path_length()))
                writeFile.write('{}\t{}\t'.format(avgShortestPath(nw,weightAttr),nw.diameter(weights=nw.es[weightAttr])))
                writeFile.write('{}\t'.format(len(nw.components().giant().vs)))
                writeFile.write('{}\n'.format(len(nw.vs)))
            else:
                break
            printCounter += 1
    print "Done."


######################### avg shortest path #########################
# Uses get_shortest_paths function to try and retrieve avg shortest
# path of graph using user-defined attribute as weights
# Used in nodeRemove, returns a numpy float type value

def avgShortestPath(network,attrName,replaceValue=False):
    weightMatrix = []
    if(replaceValue):
        for i in range(0,len(network.vs)):    
            sPathList = network.shortest_paths(weights=network.es[attrName], source = i)
            weightMatrix.append([replaceValue if x==float('inf') else x for x in sPathList])
    else:
        weightMatrix = network.shortest_paths(weights=network.es[attrName])
    avgShortPath = np.sum(weightMatrix)/(len(weightMatrix)*(len(weightMatrix)-1))
    return avgShortPath

# lighter version
def avgShortestPathV2(nw,wAttr):
    n = nw.shortest_paths(weights = wAttr)
    n2 = [[x for x in i if x != float('inf')] for i in n]
    return np.mean(sum(n2,[]))

######################### generate X rounds of random remove #########################

def testRandom(network, outputName,weight, rounds):
    for i in range(0,rounds):
        output = outputName+"_"+str(i)+".txt"
        randNodeRemove(network, output, weight)

######################### Plot distribution using random networks #########################

def randomDist(network, weight, rounds):    
    shortestPathList = []
    for i in range(0,rounds):
        if(i%10 == 0 and i > 0):
            print "{} rounds completed.".format(i)
        randG = randomiseEdge(network, weight, 'Position')
        shortestPathList.append(avgShortestPath(randG, weight))        
    return shortestPathList

######################### Edge file edits for Gephi #########################
# Using the current gephi edge file, outputs a new edited file to
# showcase how the network looks like with edge removed until network broke

def removeEdgeFromFile(rankRemove, edgeFile, outputFile, delim):
    skipHeader = False
    toRemove = []
    toWrite = []
    firstInf = True    
    # Collect rankRemove data
    with open(rankRemove, 'rb') as fileHandle:
        rankFile = csv.reader(fileHandle, delimiter = delim)
        for line in rankFile:
            if(line[2] == 'inf'):
                if(firstInf):
                    firstInf = False
                else:
                    break
            toRemove.append(line[0]) #list of strings
    # Check the edge to see if it connects to any of the toRemove nodes
    with open(edgeFile, 'rb') as fileHandle:
        edges = csv.reader(fileHandle, delimiter = delim)
        for line in edges:            
            if(not skipHeader):
                skipHeader = True                
                toWrite.append(line)
                continue
            pos1 = line[0]
            pos2 = line[1]
            if(pos1 in toRemove or pos2 in toRemove):
                continue
            elif(pos1 in toRemove and pos2  in toRemove):
                continue
            else:
                toWrite.append(line)
    with open(outputFile, 'wb') as outHandle:
        for line in toWrite:
            outHandle.write("{}\t{}\t{}\t".format(line[0],line[1],line[2]))
            outHandle.write("{}\t{}\t{}\t".format(line[3],line[4],line[5]))
            outHandle.write("{}\t{}\t{}\n".format(line[6],line[7],line[8]))

######################### Node file edits for Gephi #########################
# Using the current gephi node file, outputs a new edited file to
# showcase how the network looks like with node removed until network broke

def removeNodeFromFile(rankRemove, nodeFile, outputFile, delim):    
    toRemove = []
    toWrite = []    
    skipHeader = False
    firstInf = True       
    # Collect rankRemove data
    with open(rankRemove, 'rb') as fileHandle:
        rankFile = csv.reader(fileHandle, delimiter = delim)
        for line in rankFile:
            if(line[2] == 'inf'):
                if(firstInf):
                    firstInf = False
                else:
                    break
            toRemove.append(line[0])
            
    #check which nodes to remove
    with open(nodeFile, 'rb') as fileHandle:
        nodes = csv.reader(fileHandle)
        for line in nodes:
            if(not skipHeader):
                skipHeader = True                
                toWrite.append(line)
                continue
            pos = line[0]
            if(pos in toRemove):
                continue
            else:
                toWrite.append(line)
    with open(outputFile, 'wb') as outHandle:
        for line in toWrite:
            outHandle.write("{},{},{},".format(line[0],line[1],line[2]))
            outHandle.write("{},{},{},".format(line[3],line[4],line[5]))            
            outHandle.write("{},{},{},".format(line[6],line[7],line[8]))
            outHandle.write("{},{},{},".format(line[9],line[10],line[11]))
            outHandle.write("{},{},{},".format(line[12],line[13],line[14]))
            outHandle.write("{},{}\n".format(line[15],line[16]))

######################### Clustering stats #########################            
# Uses newman leading eigenvector to calculate clusters
# Also gives clustering coefficient

def clusterStats(network, weight, outputName):
    nw = network.copy()
    nwClusters = nw.community_leading_eigenvector(weights = nw.es[weight])
    with open (outputName, 'wb') as writeFile:    
        writeFile.write("Network modularity score = {}\n".format(nwClusters.q))
        for i in range(0,len(nwClusters)):
            subG = nwClusters.subgraph(i)    
            coeffList = subG.transitivity_local_undirected(weights=subG.es[weight])    
            subG_coeff = np.average(coeffList)
            writeFile.write("Subgraph {}\tClusteringCoeff = {}\t".format(i, subG_coeff))
            writeFile.write("Subgraph size = {}\n".format(len(subG.vs)))


######################### Scores for nodes and node ranking method #########################            
# Our in-house idea of ranking nodes in the network based on neighbours using 
# the combination method of what has been implemented. First function sets a penalty score
# based on neighbour and is used in the second function to provide ranks. 
# Use in nodeStat function. Note, even though functions deals with to/from, network can
# still be undirected (for this case it is undirected).

def setNodePenalty(percentile,rankList,m=100):
    benchMark = np.percentile(rankList,percentile)
    lower = min(rankList) - benchMark
    higher = max(rankList) - benchMark
    transformed = []
    
    for i in range(0,len(rankList)):
        s = rankList[i] - benchMark
        transformed.append(m*((s-lower)/(higher - lower)))
    return transformed

#--------------------------------------------------------------------------#

def chatuRankAll(network, weight, percentile, listOutput=False):
    nw = network.copy()
    avgMIRank = allNode_avgMI(nw)
    penaltyList = setNodePenalty(percentile,avgMIRank)

    if(not listOutput):
        rankDict = {}
        for nodeIndex in range(0,len(avgMIRank)):
            aapos = nw.vs[nodeIndex]['Position']
            penaltySum = 0
            sourceList = nw.es.select(_source=nodeIndex) #all the paths coming out of aapos
            targetList = nw.es.select(_target=nodeIndex) #all the paths going into aapos
            for neighbourIndex in [x.target for x in sourceList]:
                penaltySum += penaltyList[neighbourIndex]
            for neighbourIndex in [x.source for x in targetList]:
                penaltySum += penaltyList[neighbourIndex]
            rankDict[aapos] = penaltySum/nw.vs[nodeIndex].degree()
        return rankDict
    elif(listOutput):
        rankList = []
        for nodeIndex in range(0,len(avgMIRank)):
            #aapos = nw.vs[nodeIndex]['Position']
            penaltySum = 0
            sourceList = nw.es.select(_source=nodeIndex) #all the paths coming out of aapos
            targetList = nw.es.select(_target=nodeIndex) #all the paths going into aapos
            for neighbourIndex in [x.target for x in sourceList]:
                penaltySum += penaltyList[neighbourIndex]
            for neighbourIndex in [x.source for x in targetList]:
                penaltySum += penaltyList[neighbourIndex]
            rankList.append(penaltySum/nw.vs[nodeIndex].degree())
        return rankList


######################### Check directly connected neighbours #########################
# Helper functions to aid in retrieve stats more easily
# because iGraph nodeID is a pain to deal with. Note, even though
# functions deals with to/from, network can still be undirected
# (for this case it is undirected).
# returns a list of neighbours given the position.

def checkNeighbours(network, pos):
    nodeID = network.vs.select(Position_eq = str(pos))[0].index

    #list of nodes that nodeID points to
    from_nodeID = [e.target for e in network.es.select(_source=nodeID)]

    #list of nodes that points to nodeID
    to_nodeID = [e.source for e in network.es.select(_target=nodeID)]

    #transform back to positions and put into dictionary
    from_nodeID2 = [network.vs[v]['Position'] for v in from_nodeID]
    to_nodeID2 = [network.vs[v]['Position'] for v in to_nodeID]
    return from_nodeID2 + to_nodeID2



######################### Track avg short path length #########################
# track avg path length while adding nodes
# also tracks number of edges added while adding nodes
# The latter is interesting because once a new node is added,
# all other edges connecting to nodes in the graph will be sequetially added
# between an addition of a new node, edge count will rise without node count increasing.


# update coulnd't get this to work 100%. Some flaws. probalby won't use.
def trackAvgShortPathLen(fileName, delim):
    skipHeader = False
    pos_nodeID_dict = {}
    #nodeID_pos_dict = {}
    pos_edgeID_dict = {}
    g = ig.Graph()

    #avg shortest path list. index of list is number of nodes in the graph.
    #The avg shortest path len is prefilled with two zeroes, because
    #no nodes and one node in graph can't calculate avg path length.
    avgSPList = [0,0]
    listCount = len(avgSPList)-1 #

    edgeCountList = [0,0]

    with open(fileName, 'rb') as csvfile:
        fileHandle = csv.reader(csvfile, delimiter= delim)
        edgeID = 0
        nodeID = 0

        #Retrieve node (pos1,pos2) data and
        #edge data(pvalue, mi_weight, one_minus_mi)
        #Customise this part if needed.
        for row in fileHandle:
            if(not skipHeader):
                skipHeader = True
                continue
            pos1 = row[0]
            pos2 = row[1]
            mi_weight = row[2]
            pValue = row[3]
            one_minus_mi = row[4]

            # Dealing with 2 nodes, igraph nodeID is sequential
            # starts from 0 to n. Makes sure that the nodes added
            # to graph have not been added before
            if(pos1 not in pos_nodeID_dict.keys()):
                pos_nodeID_dict[pos1] = nodeID
                #nodeID_pos_dict[nodeID] = pos1
                g.add_vertices(1)
                g.vs[nodeID]['Position'] = pos1
                nodeID += 1
            if(pos2 not in pos_nodeID_dict.keys()):
                pos_nodeID_dict[pos2] = nodeID
                #nodeID_pos_dict[nodeID] = pos2
                g.add_vertices(1)
                g.vs[nodeID]['Position'] = pos2
                nodeID += 1

            #Dealing with the edge between the 2 nodes
            if((pos1,pos2) not in pos_edgeID_dict.keys()):
                pos_edgeID_dict[(pos1,pos2)] = edgeID

            #Stored the id of nodes of positions and
            #stored the id of edges between positions.
            #Time to add this edge into the graph.
            g.add_edges([(pos_nodeID_dict[pos1],pos_nodeID_dict[pos2])])

            #Add some attributes to the edge
            #Customise this part if needed.
            #g.es[edgeID]['MI'] = float(mi_weight)
            #g.es[edgeID]['PValue'] = float(pValue)
            g.es[edgeID]['1-MI'] = float(one_minus_mi)
            edgeID += 1

            
            #check if there's a new node in the network
            if(len(pos_nodeID_dict) > listCount):
                edgeCountList.append(len(pos_edgeID_dict))
                avgSPList.append(avgShortestPath(g,"1-MI"))
                listCount = len(pos_nodeID_dict)


    return (avgSPList, edgeCountList)


######################### print avg short path length stats #########################

def writeNetworkBuildStat(avgSPList, edgeCountList, filePath):

    assert(len(avgSPList) == len(edgeCountList))
    with open (filePath, 'wb') as writeFile:
        writeFile.write("NumNodes,AvgShortPath,NumEdges\n")
        for i in range(0,len(avgSPList)):
            writeFile.write(str(i)+","+str(avgSPList[i])+","+str(edgeCountList[i])+"\n")





