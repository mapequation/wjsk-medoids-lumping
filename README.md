# Divisive hierarchical state node lumping with k-medoids++ based on weighted Jensen-Shannon divergence

## Author:

Martin Rosvall

For contact information, see http://www.mapequation.org/about


## Getting started:

In a terminal with the GNU Compiler Collection installed,
just run 'make' in the current directory to compile the
code with the included Makefile.


Call: ./dangling-lumping [-s \<seed\>]  [-N \<number of attempts\>] [-k \<number of clusters\>] [-d \<number of clusters in each division (>= 2)\>] [--tune] [--batchoutput] input_state_network.net output_state_network.net  
seed: Any positive integer.  
number of clusters: The preferred number of clusters per physical node. Default is 100.  
nunmber of attempts: The number of attempts to optimize the cluster assignments of each physical node. Default is 1.  
number of clusters in each division (>= 2): The number of clusters the cluster with highest weighted Jensen-Shannon divergence will be divided into. Default is 2.  
--tune: Iteratively tunes medoids in deepest hierarchical level with pivot method. Default is false.  
--batchoutput: Writes the output in batched format if input is in batched format. Default is false.  
number of random state nodes in medoid update: Size of randomly sampled set of state nodes in medoid for finding a new center. 0 for skipping completely. 200 seems to give a good balance between speed and accuracy. Default is full medoid size.   
input_state_network.net: A state network without dangling state nodes, batched or non-batched, looks for repeated instances of *States, *Links, and *Contexts  
output_state_network.net: The lumped state network where all state nodes in each physical node have been lumped into at most k clusters with k-medoids++ based on weighted Jensen-Shannon divergence  