### OVERVIEW
This repository contains code to reproduce the experiments in the NeurIPS 2021 paper [1]. If you use this code in your own work, you are required to cite [1].

[1] F. Bernard, D. Cremers, J. Thunberg. Sparse Quadratic Optimisation over the Stiefel Manifold with Application to Permutation Synchronisation. NeurIPS 2021

### RUNNING THE CODE  
 Simply run the file
     run_evaluation_house.m
It should not take more than 10-15 minutes to complete, depending on your hardware 
specifications.

The code was successfully tested on Mac OS (11.2).


### EXTERNAL CODE AND LICENSES  
This package includes (parts of) the following external code:
- MatchEig
     - source: http://www.diegm.uniud.it/fusiello/demo/mvm/
- MatchALS: 
     - source: https://github.com/zju-3dv/multiway
     - notes: by default MatchALS is not evaluated since it is slow. If you'd
       like to also run MatchALS, we recommed to install the MOSEK optimisation 
       toolbox ( https://www.mosek.com ), which provides an efficient implementation
       of the function quadprog(), which MatchALS uses internally in proj2kav().
       If MOSEK is not installed, you will not be able to reproduce the reported 
       runtimes of MatchALS, as MatchALS will be much slower.
- NmfSync: 
     - source: https://github.com/fbernardpi/NmfSync
- Dataset:
     - source: http://www.cs.cmu.edu/afs/cs/project/vision/vasc/idb/www/html/motion/ and http://pages.cs.wisc.edu/~pachauri/perm-sync/



### ISSUES 
If you face proplems on Mac OS with unverified mex files, you can resolve 
these errors by calling the following command from the main folder. 
     sudo xattr -r -d com.apple.quarantine ./external/NmfSync-master/fastAuction_v2.6/
     sudo spctl --add ./external/NmfSync-master/fastAuction_v2.6/auctionAlgorithmSparseMex.mexmaci64

Alternatively, you can recompile the file auctionAlgorithmSparseMex.cpp in 
the folder NmfSync-master/fastAuction_v2.6 on your own using the command
     mex -largeArrayDims auctionAlgorithmSparseMex.cpp -lut
