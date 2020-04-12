# Multi-objective-inversion-of-CSEM-data

This repository contains the source code for implementing the multi-objective inversion of CSEM data using Non dominated sorting genetic algorithms (NSGA-II).

The results of this code are published in the corresponding article:https://academic.oup.com/gji/article-abstract/220/2/1066/5613371

This code uses the forward modeling subroutines from the Occam1DCSEM as described by Key (2009). The following is the link to the paper: https://library.seg.org/doi/abs/10.1190/1.3058434

## Pre-requisites
This code uses MPI libraries,therefore they should be installed first on your local machine. (I used 'brew install open-mpi' on my mac to install these libraries)

## Code compilation instruction for mac/linux
1. Download the code and open the terminal. cd to the source folder.

2. type $make  on your terminal and it should create the exectuble named "main". Ideally it should be compiled without error. However, if you get the following error:
  "FATAL:/opt/local/bin/../libexec/as/x86_64/as: I don't understand 'm' flag!"

  Copy the following lines on your terminal:
  
  directory_to_remove=/opt/local/bin 
 
  PATH=:$PATH: 
 
  PATH=${PATH//:$directory_to_remove:/:} 
 
  PATH=${PATH#:}; PATH=${PATH%:}
 
  Then type $make and it should be good to go.

3. type $mpirun -np 4 ./main ( You can use more than 4 cores if you are running in on the cluster to speed up the  performance)

### input files to the executable
RUNFILE: this file contains all the input parameters like Transmitter co-ordinates, No of frequencies used, no of layers and their corresponding resistivity values.

