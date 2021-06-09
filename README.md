# SBM - Stochastic Boltzman Machine Inference
Stochastic Quasi-Newton Optimizer using Boltzman Machine Updates for Protein Models

## Install
To use the SBM on Matlab one must perform Mex operation on the c++ files in the Functions folder.
in matlab, when in the Functions folder run 
```
mex C_MSA_MonteCarlo.cpp
```

In MacOs machines there might be a compatability issue with the compiler (Clang), this will menifest in Matlab crashing when running this.
To overcome this we must force the compilation to use gcc. In the file mex_openmp.sh change the `MATLABROOT` value to whatever the matlabroot is,
you can see that by running `matlabroot` in matlab. Afterwards, in the terminal, with appropriate permissions, run 
```
./mex_openmp.sh C_MSA_MonteCarlo.cpp
``` 

## Alignment
The code includes a line to load a .mat file with an MSA in a form of a numeric matrix called "align" in the data folder.
To create this alignment from a fasta file just run in the matlab command line
```
align=fasta2num(fasta_file);
```
and save the align into a matlab file placed inside the data folder.

Alternatively, just replace line 26 in the RunFile with the line above to circumvent the need for the mat file and use only fasta files instead. 
