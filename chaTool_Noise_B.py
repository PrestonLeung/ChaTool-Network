#!/usr/local/bin/python

import re, os, argparse, sys
import CommonFunctions as CF
import random
import sklearn.metrics as skm
from datetime import datetime as dt
import numpy as np


#==========================="Global Variables"===========================#

NumOfShuffle = 1000
PValue = 0.05
ShufflePortion = 1.0
Window = 100
NoHeader = False

#==========================="parse fasta file"===========================#
# Make fasta record as a list instead of a generator.

def parseFasta(fastaFile):
    
    fastaList = []    
    fastaRecord = CF.getFile(fastaFile)
    
    for line in fastaRecord:        
        fastaList.append(list(line.seq))
    return fastaList

#==========================="shuffle MI"===========================#
# Takes a pair of positions as lists. Then shuffles ShufflePortion% of the
# list n times and recalculates MI score for each shuffle. Compares 
# observed MI score versus shuffled MI score. The total count of shuffled
# MI score which is greater than observed MI score, is less than Pvalue 
# over n times of shuffling, then the correlation is deemed non random.


def shuffle_MI(fastaList, start, end):

    seqLen = len(fastaList[0])
    fastaLen = len(fastaList)
    
    with open(str(start+1)+"_"+str(end)+".txt", 'wb') as writeHandle:
    
        if(not NoHeader):
            writeHandle.write("Pos1\tPos2\tMI\tP-Value\t1-MI\n")
        
        for indexA in range(start,end):
            if(indexA < seqLen):                
                indexB = indexA + 1                
                vectorX = list(np.array(fastaList)[:,indexA])
                
                while(indexB < seqLen):
                    vectorY = list(np.array(fastaList)[:,indexB])                    
                    originalMI = skm.mutual_info_score(vectorX,vectorY)                    
                    countMI = 0
                     
                    for trials in range(0, NumOfShuffle):
                        
                        if(ShufflePortion < 1.0):                        
                            tempVectorX = shuffleVector(vectorX)
                            tempVectorY = shuffleVector(vectorY)
                        else:
                            tempVectorX = list(vectorX)
                            tempVectorY = list(vectorY)                        
                            random.seed(dt.now().time())
                            random.shuffle(tempVectorX)
                            random.shuffle(tempVectorY)                        
                        
                        randomMI = skm.mutual_info_score(tempVectorX, tempVectorY)
                        
                        if(randomMI >= originalMI):
                            countMI += 1
                            
                    if(originalMI > 0 and (float(countMI)/NumOfShuffle) < PValue):
                        writeHandle.write(str(indexA+1)+"\t"+str(indexB+1)+"\t")
                        writeHandle.write(str(originalMI)+"\t"+str(float(countMI)/NumOfShuffle)+"\t")
                        writeHandle.write(str(1-originalMI)+"\n")
                    indexB += 1    

        
#==========================="shuffle vector"===========================#

def shuffleVector(vectorS):
    vlen = len(vectorS)
    portionToShuffle = int(ShufflePortion * vlen)
    tempVec = list(vectorS)
    
    if(portionToShuffle < 1):
        raise CF.InputError("Given portion to shuffle is less than 1.")
    else:
        random.seed(dt.now().time())
        portionIndex = random.sample(range(0, vlen), portionToShuffle)
        shuffleV = list(np.array(vectorS)[portionIndex])
        random.shuffle(shuffleV)
        
        for i in range(0,len(portionIndex)):
            tempVec[portionIndex[i]] = shuffleV[i]
    
    return tempVec
    
#==========================="main function"===========================#

if __name__ == '__main__':
    
    bigHash = {}
    seqLength = None
    start_pos = 1
    program_function = """
    *****
    
        Takes a pair of positions as lists. Then shuffles the list n times
        and recalculates MI score for each shuffle. Compares observed MI
        score versus shuffled MI score. If the total number of count of 
        shuffled MI score which is greater than observed MI score, is less than
        5% [default] over n times of shuffling, then the correlation is deemed 
        non random. This is done across entire sequence.
        
        Note:   Parameter for start is a number that split up the pair-wise 
                comparison so multiple instances of the tool can be run in 
                different sections given a simple series of integers of start.
                Best to be run on qsub. This reduces long waits for results. 
                
                E.g. Start list = [1,2,3,4,.....,30]
                Window = 100 [Default]
                S = (Start*100)-100
                E = S*100
                
                
                At Start = 1:
                ------------------------------------------------------------
                S             E
                
                Section Covered:
                S-----------------------------------------------------------
                |2----------------------------------------------------------
                | 3---------------------------------------------------------
                |             .
                |             . 
                |             .
                <---Window--->E---------------------------------------------
                
                Stop.
                
                At Start = 2:
                ------------------------------------------------------------
                               S             E
                Section Covered at Start = 2:
                               S--------------------------------------------
                               |2-------------------------------------------
                               | 3------------------------------------------
                               |             .
                               |             . 
                               |             .
                               <---Window--->E------------------------------
                Stop.
                
                
                Hence MI has been tested from S to E across the sequence until 
                end of start list.
                 
    Ver. 0.01
    
    Written by PresDawgz (~w~)v
    
    *****
    """
    parser = argparse.ArgumentParser()
    #required arguments
    parser.add_argument('-f', '--file', help='Filename of FASTA file.',  required=True)    
    parser.add_argument('-s', '--start', help='Start position.', required=True, type=int)

    #optional arguments
    parser.add_argument('-sf', '--shuffle', help="Number of iterations to perform. Default = 1000.", type=int)
    parser.add_argument('-p','--pvalue', help="Set p-value. Default 0.05.", type=float)
    parser.add_argument('-sp', '--shuffle_portion', help='Declare how much percent of shuffle to do. Default = 1.0.', type=float)  
    parser.add_argument('-w', '--window', help='Length to do the shuffling. Min window is 1. Default = 100', type=int)
    parser.add_argument('-nh','--no_header',help = 'Exclude header.', action='store_true')
    
    if len(sys.argv) <= 1:
            print(program_function)
            parser.print_help()
            sys.exit()
    
    args=parser.parse_args()
    if(args.shuffle):
        NumOfShuffle = args.shuffle
    
    if(args.pvalue):
        if(args.pvalue > 0.0):
            PValue = args.pvalue
        else:
            print "Given p-value is {} and is invalid. Proceeding with default value".format(args.pvalue)    
        
    if(args.shuffle):
        if(args.shuffle_portion > 0.0):
            ShufflePortion = args.shuffle_portion
        else:
            print "Given shuffle portion is {} and is invalid. Proceeding with default value".format(args.shuffle_portion)    
    if(args.window):
        if(args.window >= 1):
            Window = args.window
        else:
            print "Given window is {} and is invalid. Proceeding with default value".format(args.window)    
    if(args.no_header):
        NoHeader = True
    
    fastaList = parseFasta(args.file)
    
    end=args.start*Window
    start=args.start*Window-Window    
    
    shuffle_MI(fastaList,start, end)
    
    
    
    

    
  
  
