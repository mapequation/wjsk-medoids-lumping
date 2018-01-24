# Divisive hierarchical state node lumping with k-medoids++ based on weighted Jensen-Shannon divergence

## Author:

Martin Rosvall

For contact information, see http://www.mapequation.org/about


## Getting started:

In a terminal with the GNU Compiler Collection installed,
just run 'make' in the current directory to compile the
code with the included Makefile.


Call: ./dangling-lumping [-s \<seed\>]  [-N \<number of attempts\>] [-k \<max number of clusters\>] [-e \<max entropy in lumped state node\>] [--max-order \<max Markov order\>] [-d \<number of clusters in each division (>= 2)\>] [--fast] [--batchoutput] [--state-containers state_assignments.txt] [--context-containers lumped_state_network.net] input_state_network.net output_state_network.net  
seed: Any positive integer.  
max number of clusters: The max number of clusters per physical node if max entropy is not met. Default is 100.
max entropy in lumped state node: The max entropy rate in lumped state node if max number of clusters is not met. Default is no limit.  
max Markov order: First lump based on contexts to max Markov order.  
nunmber of attempts: The number of attempts to optimize the cluster assignments of each physical node. Default is 1.  
number of clusters in each division (>= 2): The number of clusters the cluster with highest weighted Jensen-Shannon divergence will be divided into. Default is 2.  
--fast: Will use center to approximate meodid (no lumping during cluster assignment). Default is false.  
--batchoutput: Writes the output in batched format if input is in batched format. Default is false.  
number of random state nodes in medoid update: Size of randomly sampled set of state nodes in medoid for finding a new center. 0 for skipping completely. 200 seems to give a good balance between speed and accuracy. Default is full medoid size.  
state_assignments.txt: Lump in containers set by state-to-container assignemnts. Default is to lump in physical nodes.  
lumped_state_network.net: Lump in containers set by context assignemnts. Default is to lump in physical nodes.  
input_state_network.net: A state network without dangling state nodes, batched or non-batched, looks for repeated instances of *States, *Links, and *Contexts  
output_state_network.net: The lumped state network where all state nodes in each physical node have been lumped into at most k clusters with k-medoids++ based on weighted Jensen-Shannon divergence  