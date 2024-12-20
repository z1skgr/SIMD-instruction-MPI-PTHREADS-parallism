# SIMD-instructions/parallism
>   Parallelism standards for accelerating performance on calculations for detection of positive DNA selection



## Table of contents
* [General Info](#general-information)
* [Features](#features)
* [Prerequisites](#prerequisites)
* [Setup](#setup)
* [How to run](#how-to-run)
* [Acknowledgements](#acknowledgements)


## General Information
__Streaming SIMD Extensions (SSE), MPI and Pthreads__
for parallelization the calculation of a simplified form of OMEGA
statistical, to detect positive selection in DNA sequences.

Repeated for a set of _N DNA_ sites, extract statistical OMEGA outputs

| $$num = \frac{L+R}{\left(\frac{m*(m-1.0)}{2.0}\right) + \left(\frac{n*(n-1.0)}{2.0}\right)}$$ |
|---|
| $$den= \frac{C-L-R}{m*n}$$  |
| $$\omega= \frac{num}{den+0.01}$$  |

The above calculations are repeated iteratively for a set of N DNA sites, 
$N$ is a variable for which $N ≥ 1$ . All input data
is arranged in arrays of length N (user variable).

Interested for 
| $$max_{\omega}$$   | $$min_{\omega}$$  | $$avg_{\omega}$$  |
|---|---|---|



_BONUS:_ Implementation of a different memory layout for better 
performance with _SSE_ commands. 

## Features
* Reference
* SSE instructions 
* Parallel standards (pthreads/MPI) + SSE

The implementation was guided by the reference file
Benchmarked on Intel(R) Core(TM) i7-1065G7 @ 1.30GHz 1.50 GHz with 8GB DDR3 memory.

## Prerequisites 
[SIMD](https://software.intel.com/sites/landingpage/IntrinsicsGuide) <br>
[MPI](http://mpitutorial.com/tutorials/)



## How to run

1. GCC Installation
```
$ gcc --version
$ sudo apt install gcc
$ sudo apt install lam-runtime

```
Supposing you can navigate the proper folder to compile the desired source code from the terminal

2. MPI Installation
Download the source code [MPI](https://www.mpich.org/downloads/) release
```
$ tar -xzf mpich-X-X-X.tar.gz
$ cd mpich-X-X-X

```
where `X-X-X` mpi version

```
$ ./configure
$ make; sudo make install
```

3. Check MPI version

```
$ mpiexec --version
```

4. LAM installation
```
sudo apt-get install lam-runtime
```

### Reference

```
gcc -o reference reference.c
./reference 'N'
```

### SSE
```
gcc -o jam jam.c
./jam 'N'
gcc -o unroll unroll.c
./unroll 'N'
```
Add -msse4.2 library
```
gcc -o SSE SSE.c -msse4.2
./SSE 'N'
```

### Pthreads
```
gcc -pthread SSE_pthreads.c -o SSE_pthreads -msse4.2
./SSE_pthreads 'N' 'N_threads'
```
### MPI
```
mpicc -pthread SSE_Mpthreads.c -o SSE_Mpthreads -msse4.2 -lm
lamboot -v host
mpiexec -n 'N_proc' ./SSE_Mpthreads 'N' 'N_threads'
```

### Bonus 
```
gcc -o bonus bonus.c -msse4.2
./bonus 'N'
```
## Setup
Script  variables initialized as:
* N = 10000000. 
* Threads = [2 4]
* Processors = [2 4].


## Acknowledgements
* This project was created for the requirements of the lesson Architecture of Parallel and Distributed Computers

