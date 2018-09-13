# Parallelized-Eager-Prim
SURG2018 Project

1. Introduction
As a project for SURG 2018 undergraduate Research, I parallelized Eager Prim by design the algorithm that build multiple partial minimum spanning tree from different start points with different threads, and merge threads when any 2 threads meet. 
Code written in C++, OpenMP is used for parallelization. 

2. Files
The cpp file and the header file are in the "Code" folder, some example test cases are provided in the "test-cases" folder.
Number of cores, droprates, debug sessions are defined as Macros. To make it easy for individual test cases, this version of the code use the naive way to read files. The test-case filename have to be changed manually in the scanf statement. 

3. Results
Basic design and performance analysis are shown in the Poster "JingPosterV6.pdf". 
