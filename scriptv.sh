#!/bin/bash
# Description :	Compile and run all the variations of the algorithm witn the same inputs.

#### Find relative directory
echo $0
full_path=$(realpath $0)
dir_path=$(dirname $full_path)

######

echo ""
echo "-------->The script starts<-----"
echo ""

#Compile the code.
gcc -o $dir_path/reference/reference $dir_path/reference/reference.c
gcc -o $dir_path/SSE/jam $dir_path/SSE/jam.c
gcc -o $dir_path/SSE/unroll $dir_path/SSE/unroll.c
gcc -o $dir_path/SSE/SSE $dir_path/SSE/SSE.c -msse4.2
gcc -pthread $dir_path/pthreads/SSE_pthreads.c -o $dir_path/pthreads/SSE_pthreads -msse4.2
mpicc -pthread $dir_path/MPI/SSE_Mpthreads.c -o $dir_path/MPI/SSE_Mpthreads -msse4.2 -lm
gcc -o $dir_path/bonus/bonus $dir_path/bonus/bonus.c -msse4.2

#Run the code.
echo ""
echo "------------>Reference execution<-------------"
echo ""


$dir_path/reference/reference '10000000'
echo ""
echo "------------>Jam execution<-------------"
echo ""
echo ""

$dir_path/SSE/jam '10000000'
echo ""
echo ""

echo "------------>Unroll execution<-------------"
echo ""

$dir_path/SSE/unroll '10000000'
echo ""
echo ""


echo "------------>SSE executions.<-------------"
echo ""


$dir_path/SSE/SSE '10000000'
echo ""
echo ""
echo "------------>SSE executions + Pthread<-------------"
echo "------------->2 threads<-------------"
echo ""
$dir_path/pthreads/SSE_pthreads '10000000' '2'
echo ""

echo "-------------> 4 threads<-------------"
echo ""
$dir_path/pthreads/SSE_pthreads '10000000' '4'



lamboot 
echo ""
echo "------------->SSE executions + Pthread + MPI<-------------"
echo "------------>   2 threads 2 procs<------------"
echo ""
mpiexec -n 2 $dir_path/MPI/SSE_Mpthreads '10000000' '2'
echo ""
echo ""

echo '------------>  2 threads 4 procs<------------'
echo ""
mpiexec -n 4 $dir_path/MPI/SSE_Mpthreads '10000000' '2'
echo ""
echo ""

echo '------------>  4 threads 2 procs<------------'
echo ""
mpiexec -n 2 $dir_path/MPI/SSE_Mpthreads '10000000' '4'
echo ""
echo ""

echo "------------>4 threads 4 procs<------------"
echo ""
mpiexec -n 4 $dir_path/MPI/SSE_Mpthreads '10000000' '4'
echo ""

echo "-------------> Bonus <-------------"
echo "-----------------------------------"
$dir_path/bonus/bonus '10000000'
echo ""

echo ""
echo "-------->The script ends<-----"
echo ""


lamhalt
#Shutdown the LAM/MPI run-time environment. 
