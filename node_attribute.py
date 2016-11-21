import ChaToolNetworkFunctions as cnf
import os, argparse, sys
import igraph as ig

#==========================="Global Variables"===========================#

Delim = '\t'
    
#==========================="Main"===========================#

if __name__ == '__main__':
    
    program_function = """
    *****
    
    Calculates node attribute statistics.    
                 
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

    if len(sys.argv) <= 1:
            print(program_function)
            parser.print_help()
            sys.exit()
    
    args=parser.parse_args()
    
    if(args.delimiter):
        Delim = args.delimiter
    
    network = cnf.parseData(args.file, Delim)
    
    cnf.nodeStats2(network, args.output+'_nodeAttributes.txt',args.inf_file, Delim)
    cnf.clusterStats(network,'MI',args.output+'_clusterStats.txt')


