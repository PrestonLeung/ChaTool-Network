# CoMNet (Co-occurring Mutation NETwork)

This python tool takes in aligned fasta sequences of the same length and calculates pairwise mutual information between sites that mutate across the alignment using each position in the alignment as a variable. It performs an exhaustive pairwise mutual information calculation from position 1->N, and then 2->N (where N is the length of the alignment) etc, to a total of N(N-1)/2 calculations. These values are essentially used as edge weights in the network between two positions. The tool also filters out mutual information where connections between two positions are considered random (MIEdgeMaker.py.

Node attributes, such as node centralities within the network are also performed (nodeattribute.py). This tool will provide a table (.csv), of every node in the network and its node centrality calculations.

NodeRemove.py will test the the network by sequentually removing nodes according the node centralities calculated in nodeattribute.py from the most important to the least. NodeRemove.py will measure average shortest path and network diameter per each node removed.

NodeRemovePlot.R is just an extra script I wrote to plot out what happens to average shortest path and network diameter when you remove nodes.
