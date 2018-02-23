import ChaToolNetworkFunctions as cnf
import os, argparse, sys
import igraph as ig

#==========================="Global Variables"===========================#

Delim = '\t'
Rounds = 10
    
#==========================="Main"===========================#

if __name__ == '__main__':
    
    program_function = """
    *****
    
    Tests network robustness by removing nodes sorted by:
    
    -> Eigen Vector Centrality
    -> Node Average Mutual Information (MI)
    -> Degree
    -> Betweenness Centrality
    -> Closeness Centrality
    -> Influential Ranking
    -> Random
    
    Ver. 0.01
    
    Written by PresDawgz (~w~)v
    
    *****
    """
    parser = argparse.ArgumentParser()

    #required arguments
    parser.add_argument('-f', '--file', help='Filename file containing edge information.',  required=True)    
    parser.add_argument('-inf', '--inf_file', help='De-referenced influential node file.', required = True)
    parser.add_argument('-o', '--output', help='File to write output to.', required = True)

    #optional arguments
    parser.add_argument('-d', '--delimiter', help='Delimiter, default "\t".')
    parser.add_argument('-rds', '--rounds',help='Number times a random removal of network is done, default = 10', type = int)

    if len(sys.argv) <= 1:
            print(program_function)
            parser.print_help()
            sys.exit()
    
    args=parser.parse_args()
    
    if(args.delimiter):
        Delim = args.delimiter
    
    if(args.rounds):
        Rounds = args.rounds
        
    network = cnf.parseData(args.file, Delim)

    evcScore = network.eigenvector_centrality(directed=False, weights=network.es['MI'])
    degRank = network.vs.degree()
    avgMIRank = cnf.allNode_avgMI(network)
    btwnRank = network.betweenness(weights=network.es['1-MI'],directed=False)
    clsnsRank = network.closeness(weights=network.es['1-MI'])
    infRank = cnf.infRankMatch(args.inf_file, network, '\t')
    chatuRank = cnf.chatuRankAll(network, 'MI', 50,True)

    print "Stage 1: Removing nodes according to Eigen Vector Centrality."
    cnf.nodeRemove(network, evcScore, args.output+"_evcRankRemove.txt", '1-MI')
    print "Stage 2: Removing nodes according to Degree Rank."
    cnf.nodeRemove(network, degRank, args.output+"_degreeRankRemove.txt", '1-MI')
    print "Stage 3: Removing nodes according to Node Avg MI Rank."
    cnf.nodeRemove(network, avgMIRank, args.output+"_avgMIRankRemove.txt", '1-MI')
    print "Stage 4: Removing nodes according to Betweenness Centrality."
    cnf.nodeRemove(network, btwnRank, args.output+"_btwnsRankRemove.txt", '1-MI')
    print "Stage 5: Removing nodes according to Closeness Centrality."
    cnf.nodeRemove(network, clsnsRank, args.output+"_clsnsRankRemove.txt", '1-MI')
    print "Stage 6: Removing nodes according to Influential Rank."
    cnf.nodeRemove(network, infRank, args.output+"_infRankRemove.txt", '1-MI', isInfRank = True)
    print "Stage 7: Removing nodes according to CR Rank."
    cnf.nodeRemove(network, chatuRank, args.output+"_chaturankRemove.txt", '1-MI')
    print "Stage 8: Removing nodes according to Random Selection, repeat {} times.".format(Rounds)
    cnf.testRandom(network,args.output+"_randomRemove",'1-MI',Rounds) 
    
    startAvgPathLen = network.average_path_length()
    startAvgShortPath = cnf.avgShortestPath(network,'1-MI')
    startDiameter = network.diameter(weights=network.es['1-MI'])
    print "Initial avg path length (non-weighted): {}".format(startAvgPathLen)
    print "Initial avg shortest path (weighted): {}".format(startAvgShortPath)
    print "Initial diameter (weighted): {}".format(startDiameter) 
    
