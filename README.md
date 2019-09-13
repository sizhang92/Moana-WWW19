# Multilevel-network-alignment-(Moana)
This algorithm is to align networks not only at the node level, but also at other coarse levels to discover the correspondences between clusters at different resolutions.
## Overview
The package contains the following files:
- Moana.m: the main function of the proposed Moana algorithm for multilevel network alignment
- pMMF.m: the function of parallel multiresolution matrix factorization
- FINAL_core.m: the algorithm for alignment at the coarest level
- cluster_collection.m: the function to collect cluster assignments at different levels
- evaluate_cluster_align.m: the function to evaluate cluster alignment accuracy
- FINAL_P.m: the single node-level alignment counterpart (exact version)
- greedy_match.m: greedy matching algorithm to obtain one-to-one mappings
- toy_demo.m: the demo code for users to play with
- CAGrQc.mat: the dataset of CAGrQc networks
- graclus1.2(linux): Compiled graclus algorithm for Linux users
- graclus1.2(windows): Compiled graclus algorithm for Windows users

## Usage
Please refer to the demo code file toy_demo.m and the descriptions in each file for the detailed information. 
The code can be only used for academic purpose and please kindly cite our published paper.

## Reference
Zhang, Si, Hanghang Tong, Ross Maciejewski, and Tina Eliassi-Rad. "Multilevel Network Alignment." In The World Wide Web Conference, pp. 2344-2354. ACM, 2019.
