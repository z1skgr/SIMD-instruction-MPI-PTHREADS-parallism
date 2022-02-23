#!/bin/bash
# Description :	Compile and run all the variations of the algorithm witn the same inputs.

echo ""
echo "-------->The script starts<-----"
echo ""

#Compile the code.
gcc -o reference reference.c
gcc -o jam jam.c
gcc -o unroll unroll.c
gcc -o SSE SSE.c -msse4.2
gcc -pthread SSE_pthreads.c -o SSE_pthreads -msse4.2
mpicc -pthread SSE_Mpthreads.c -o SSE_Mpthreads -msse4.2 -lm
gcc -o bonus bonus.c -msse4.2

echo "------------>Reference<-------------"
#Run the code.
echo ""
echo "------------>Reference execution<-------------"
echo ""


./reference '10000000'
echo ""
echo "------------>Jam execution<-------------"
echo ""
echo ""

./jam '10000000'
echo ""
echo ""

echo "------------>Unroll execution<-------------"
echo ""

./unroll '10000000'
echo ""
echo ""


echo "------------>SSE executions.<-------------"
echo ""


./SSE '10000000'
echo ""
echo ""
echo "------------>SSE executions + Pthread<-------------"
echo "------------->2 threads<-------------"
echo ""
./SSE_pthreads '10000000' '2'
echo ""

echo "-------------> 4 threads<-------------"
echo ""
./SSE_pthreads '10000000' '4'



lamboot -v host
echo ""
echo "------------->SSE executions + Pthread + MPI<-------------"
echo "------------>   2 threads 2 procs<------------"
echo ""
mpiexec -n 2 ./SSE_Mpthreads '10000000' '2'
echo ""
echo ""

echo '------------>  2 threads 4 procs<------------'
echo ""
mpiexec -n 4 ./SSE_Mpthreads '10000000' '2'
echo ""
echo ""

echo '------------>  4 threads 2 procs<------------'
echo ""
mpiexec -n 2 ./SSE_Mpthreads '10000000' '4'
echo ""
echo ""

echo "------------>4 threads 4 procs<------------"
echo ""
mpiexec -n 4 ./SSE_Mpthreads '10000000' '4'
echo ""

echo "-------------> Bonus <-------------"
echo "-----------------------------------"
./bonus '10000000'
echo ""

echo ""
echo "-------->The script ends<-----"
echo ""


lamhalt
#Shutdown the LAM/MPI run-time environment. 