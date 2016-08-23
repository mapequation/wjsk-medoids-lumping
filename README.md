# Hierarchical state node lumping with k-medoids++ based on weighted Jensen-Shannon divergence

## Author:

Martin Rosvall

For contact information, see http://www.mapequation.org/about


## Getting started:

In a terminal with the GNU Compiler Collection installed,
just run 'make' in the current directory to compile the
code with the included Makefile.


Call: ./dangling-lumping [-s \<seed\>] -k \<number of clusters\>  -l \<number of hierarchical levels\> --num-update-states \<number of random state nodes in medoid update\> --batchoutput input_state_network.net output_state_network.net  
seed: Any positive integer.  
number of clusters: The preferred number of clusters per physical node  
number of hierarchical levels: The number of hierarchical rounds to reach the preferred number of clusters per physical node such that the number of medoids will multiply by (number of clusters)^(1/number of hierarchical levels) in each round. Default is 1.  
--batchoutput: Writes the output in batched format if input is in batched format  
number of random state nodes in medoid update: Size of randomly sampled set of state nodes in medoid for finding a new center. 0 for skipping completely. 200 seems to give a good balance between speed and accuracy. Default is full medoid size.   
input_state_network.net: A state network without dangling state nodes, batched or non-batched, looks for repeated instances of *States, *Links, and *Contexts  
output_state_network.net: The lumped state network where all state nodes in each physical node have been lumped into at most k clusters with k-medoids++ based on weighted Jensen-Shannon divergence  