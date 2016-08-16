# State node lumping with k-medoids++ based on weighted Jensen-Shannon divergence

## Author:

Martin Rosvall

For contact information, see http://www.mapequation.org/about


## Getting started:

In a terminal with the GNU Compiler Collection installed,
just run 'make' in the current directory to compile the
code with the included Makefile.


Call: ./dangling-lumping [-s \<seed\>] -k \<number of clusters\>  --batchoutput --fast input_state_network.net output_state_network.net  
seed: Any positive integer.  
number of clusters: The preferred number of clusters per physical node  
--batchoutput: Writes the output in batched format if input is in batched format
--fast: Skip updating medoids
input_state_network.net: A state network without dangling state nodes, batched or non-batched, looks for repeated instances of *States, *Links, and *Contexts  
output_state_network.net: The lumped state network where all state nodes in each physical node have been lumped into at most k clusters with k-medoids++ based on weighted Jensen-Shannon divergence  