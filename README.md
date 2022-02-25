# SIMD-instructions/parallism
>   Parallelism standards for accelerating performance on calculations



## Table of contents
* [General Info](#general-information)
* [Features](#features)
* [Prerequisites](#prerequisites)
* [Setup](#setup)
* [How to run](#how-to-run)
* [Acknowledgements](#acknowledgements)


## General Information
Streaming SIMD Extensions (SSE), MPI and Pthreads 
for parallelization the calculation of a simplified form of OMEGA
statistical, to detect positive selection in DNA sequences.

Repeated for a set of N DNA sites, extract statistical OMEGA outputs

## Features
* Reference
* SSE instruction application
* Parallel standards (pthreads/MPI) + SSE


Benchmarked on Intel(R) Core(TM) i7-1065G7 @ 1.30GHz 1.50 GHz with 8GB DDR3 memory.

## Prerequisites 




## How to run

1. GCC installation
```
$ gcc --version
$ sudo apt install gcc
```
### Reference

1. Compile .c file
```
gcc -o newserial newserial.c
```

### Pthreads

### MPI



# Setup
Script  variables initialized as:
* N = 10000000. 
* Threads = [2 4]
* Processors = [2 4].


## Acknowledgements
* This project was created for the requirements of the lesson Architecture of Parallel and Distributed Computers

