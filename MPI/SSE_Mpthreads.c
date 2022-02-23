/*
            ======================================================
            + Name        : SSE_Mpthreads.c                      +
            + Author      : Christos Z, 			             +
            + Description : SSE extension pthreads 	mpi		     +
            ======================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <float.h>
#include <pthread.h>
#include "xmmintrin.h"
#include <mpi.h>

#define MINSNPS_B 5
#define MAXSNPS_E 20

#define BUSYWAIT 0
#define EXECUTE 1
#define EXIT 127



double gettime(void);
float randpval (void);

//Thead stuff
static pthread_t * workerThread;
static pthread_barrier_t barrier;


#define MAX(a,b) (((a)>(b))?(a):(b))
#define MAX4(a,b,c,d) MAX(0, MAX(a,MAX(b,MAX(c,d))))

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MIN4(a,b,c,d) MIN(0, MIN(a,MIN(b,MIN(c,d))))


double gettime(void)
{
    struct timeval ttime;
    gettimeofday(&ttime , NULL);
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

float randpval (void)
{
    int vr = rand();
    int vm = rand()%vr;
    float r = ((float)vm)/(float)vr;
    assert(r>=0.0f && r<=1.00001f);
    return r;
}

typedef struct
{
    //Thread Data
    int threadID;
    int threadTOTAL;
    int threadBARRIER;
    int threadOPERATION;

    //Where the loop begins and ends
    int threadSTART;
    int threadEND;

    //All the data
    __m128 * mVec;
    __m128 * nVec;
    __m128 * LVec;
    __m128 * RVec;
    __m128 * CVec;
    __m128 * FVec;

    //All the necessary values to get results
    float * maxval;
    float * minval;
    float * average;

    __m128 * avgF;
    __m128 * maxF;
    __m128 * minF;

    float thread_max;
    float thread_min;
    float thread_sum;


}dData_t;


typedef struct
{
    //Thread Data
    int threadID;
    int threadTOTAL;
    int threadBARRIER;
    int threadOPERATION;

    //Where the loop begins and ends
    int threadSTART;
    int threadEND;


    float thread_max;
    float thread_min;
    float thread_sum;

    dData_t dt;


}threadData_t;




void initializeThreadData(threadData_t * threadData, int i, int threads,int n,float * mVec1,float * nVec1,float * LVec1,float * RVec1,float* CVec1,float* FVec1)
{
    threadData->threadID=i;
    threadData->threadTOTAL=threads;
    threadData->threadBARRIER=0;
    threadData->threadOPERATION=BUSYWAIT;

    ///We find the limits of loop

    threadData->threadSTART=(n/threads)*i;
    threadData->threadEND =(n/threads)*(i+1);

    threadData->dt.mVec= (__m128 *) mVec1;
    threadData->dt.nVec= (__m128 *) nVec1;
    threadData->dt.LVec= (__m128 *) LVec1;
    threadData->dt.RVec= (__m128 *) RVec1;
    threadData->dt.CVec= (__m128 *) CVec1;
    threadData->dt.FVec= (__m128 *) FVec1;

    threadData->dt.maxval = (float*)_mm_malloc(sizeof(float),32);
    threadData->dt.minval = (float*)_mm_malloc(sizeof(float),32);
    threadData->dt.average = (float*)_mm_malloc(sizeof(float),32);


    for (int i = 0; i < 4; ++i)
    {
        threadData->dt.maxval[i]=0.0f;
        threadData->dt.minval[i]=FLT_MAX;
        threadData->dt.average[i]=0.0f;
    }

    threadData->dt.avgF= (__m128 *) threadData->dt.average;
    threadData->dt.maxF= (__m128 *) threadData->dt.maxval;
    threadData->dt.minF= (__m128 *) threadData->dt.minval;

    threadData->thread_max=0.0f;
    threadData->thread_min=FLT_MAX;
    threadData->thread_sum=0.0f;
}


void updateThread(threadData_t * threadData)
{
    for (int unsigned i=0;i<(threadData->threadTOTAL);i++){

        for (int l = 0; l < 4; ++l)
        {
            threadData[i].dt.maxval[l]=FLT_MIN;
            threadData[i].dt.minval[l]=FLT_MAX;
            threadData[i].dt.average[l]=0.0f;
        }

        threadData[i].thread_max=FLT_MIN;
        threadData[i].thread_min=FLT_MAX;
        threadData[i].thread_sum=0.0f;


    }
}



void compute(threadData_t * threadData)
{
    float avgF = 0.0f;
    float maxF = FLT_MIN;
    float minF = FLT_MAX;


    int      i= threadData->threadSTART ;
    int      end= threadData->threadEND ;

    __m128 * mVec128=threadData->dt.mVec;
    __m128 * nVec128=threadData->dt.nVec;
    __m128 * LVec128=threadData->dt.LVec;
    __m128 * RVec128=threadData->dt.RVec;
    __m128 * CVec128=threadData->dt.CVec;
    __m128 * FVec128=threadData->dt.FVec;


    //Vectors with constants.;
    __m128 vec1 = _mm_set1_ps(1.0f);
    __m128 vec2 = _mm_set1_ps(2.0f);
    __m128 vec001 = _mm_set_ps1(0.01f);

    //	Vectors for max,min,avg.
    __m128 maxFVec128 = _mm_set_ps1(FLT_MIN);
    __m128 minFVec128 = _mm_set_ps1(FLT_MAX);
    __m128 avgFVec128 = _mm_set_ps1(0);
    //	Vectors for flags.
    __m128 comparemaxFlags;
    __m128 compareminFlags;

    for(i=i/4;i<end/4;i+=1)
    {
        //	Same calculations with the serial execution. Floating point operations are non-associative and their order is respected.
        __m128 num_0 = _mm_add_ps( LVec128[i], RVec128[i]);
        __m128 num_1 = _mm_div_ps( _mm_mul_ps( mVec128[i], _mm_sub_ps( mVec128[i], vec1)), vec2);
        __m128 num_2 = _mm_div_ps( _mm_mul_ps( nVec128[i], _mm_sub_ps( nVec128[i], vec1)), vec2);
        __m128 num = _mm_div_ps( num_0, _mm_add_ps( num_1, num_2));

        __m128 den_0 = _mm_sub_ps( _mm_sub_ps( CVec128[i], LVec128[i]), RVec128[i]);
        __m128 den_1 = _mm_mul_ps( mVec128[i], nVec128[i]);
        __m128 den = _mm_div_ps( den_0, den_1);

        FVec128[i] = _mm_div_ps(num, _mm_add_ps(den, vec001));

        //	Compare verticaly and keep larger floats in a vector.
        comparemaxFlags = _mm_cmpgt_ps( FVec128[i], maxFVec128);
        compareminFlags=  _mm_cmplt_ps( FVec128[i], minFVec128);
        maxFVec128 = _mm_or_ps(_mm_and_ps( comparemaxFlags, FVec128[i]),  _mm_andnot_ps(comparemaxFlags, maxFVec128));
        minFVec128 = _mm_or_ps(_mm_and_ps( compareminFlags, FVec128[i]),  _mm_andnot_ps(compareminFlags, minFVec128));
        avgFVec128 = _mm_add_ps(avgFVec128,FVec128[i]);

    }

    __m128 tempMax = _mm_max_ps(maxFVec128, _mm_shuffle_ps(maxFVec128, maxFVec128, _MM_SHUFFLE(0,0,3,2)));
    maxF = _mm_cvtss_f32(_mm_max_ps(tempMax, _mm_shuffle_ps(tempMax, tempMax, _MM_SHUFFLE(0,0,0,1))));
    __m128 tempMin = _mm_min_ps(minFVec128, _mm_shuffle_ps(minFVec128, minFVec128, _MM_SHUFFLE(0,0,3,2)));
    minF = _mm_cvtss_f32(_mm_min_ps(tempMin, _mm_shuffle_ps(tempMin, tempMin, _MM_SHUFFLE(0,0,0,1))));
    __m128 tempAvg = _mm_add_ps(avgFVec128, _mm_shuffle_ps(avgFVec128, avgFVec128, _MM_SHUFFLE(0,0,3,2)));
    avgF = _mm_cvtss_f32(_mm_add_ps(tempAvg, _mm_shuffle_ps(tempAvg, tempAvg, _MM_SHUFFLE(0,0,0,1))));


    threadData->thread_max = maxF;
    threadData->thread_min = minF;
    threadData->thread_sum = avgF;

}

void syncThreadsBARRIER(threadData_t * threadData)
{
    threadData[0].threadOPERATION=BUSYWAIT;
    pthread_barrier_wait(&barrier);
}

void setThreadOperation(threadData_t * threadData, int operation)
{
    int i, threads=threadData[0].threadTOTAL;

    for(i=0;i<threads;i++)
        threadData[i].threadOPERATION = operation;
}


void startThreadOperations(threadData_t * threadData, int operation)
{
    setThreadOperation(threadData, operation);

    if(operation==EXECUTE)
        compute(&threadData[0]);

    threadData[0].threadBARRIER=1;
    syncThreadsBARRIER(threadData);
}

void terminateWorkerThreads(pthread_t * workerThreadL, threadData_t * threadData)
{
    int i, threads=threadData[0].threadTOTAL;

    for(i=0;i<threads;i++)
        threadData[i].threadOPERATION = EXIT;

    for(i=1;i<threads;i++)
        pthread_join(workerThreadL[i-1],NULL);
}


void * thread (void * x)
{
    threadData_t * currentThread = (threadData_t *) x;

    while(1)
    {
        __sync_synchronize();

        if(currentThread->threadOPERATION==EXIT)
            return NULL;

        if(currentThread->threadOPERATION==EXECUTE)
        {
            compute (currentThread);
            currentThread->threadOPERATION=BUSYWAIT;
            pthread_barrier_wait(&barrier);
        }
    }

    return NULL;
}

int main(int argc, char ** argv)
{
    assert(argc==3);

    float avgF = 0.0f;
    float maxF = 0.0f;
    float minF = FLT_MAX;
    float *max, *min  = NULL;


    int world_size;
    int world_rank;
    MPI_Init(NULL, NULL);
    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


    double timeTotalMainStart = gettime();
    unsigned int N = (unsigned int)atoi(argv[1]);
    int threads = (int)atoi(argv[2]);
    unsigned int iters = 1;
    int k=N%4;




    srand(1);
    float * mVec = (float*)_mm_malloc(sizeof(float)*N,32);
    assert(mVec!=NULL);
    float * nVec = (float*)_mm_malloc(sizeof(float)*N,32);
    assert(nVec!=NULL);
    float * LVec = (float*)_mm_malloc(sizeof(float)*N,32);
    assert(LVec!=NULL);
    float * RVec = (float*)_mm_malloc(sizeof(float)*N,32);
    assert(RVec!=NULL);
    float * CVec = (float*)_mm_malloc(sizeof(float)*N,32);
    assert(CVec!=NULL);
    float * FVec = (float*)_mm_malloc(sizeof(float)*N,32);
    assert(FVec!=NULL);



    //Threads
    ////////////////////////////////////////////////////////////

    //We are going to use Pthread-Barrier as our choise of synchronization
    int s = pthread_barrier_init(&barrier,NULL,(unsigned int)threads);	//We initialize the barrier
    assert(s == 0);


    workerThread=NULL;
    workerThread=(pthread_t *) malloc (sizeof(pthread_t)*((unsigned long)threads-1)); //We write the workerThread as thread-1 cause the one thread is the master

    threadData_t * threadData = (threadData_t *) malloc (sizeof(threadData_t)*((unsigned long)threads));
    assert(threadData!=NULL);

    for (int i = 0; i < threads; i++)
    {
        initializeThreadData(&threadData[i],i,threads,N,mVec,nVec,LVec,RVec,CVec,FVec);
        if(i>0)
            pthread_create(&workerThread[i-1],NULL,thread,(void*)(&threadData[i]));
    }


    ////////////////////////////////////////////////////////////

    //Initialize data
    ////////////////////////////////////////////////////////////
    for(unsigned int i=0;i<N;i++)
    {
        mVec[i] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
        nVec[i] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
        LVec[i] = randpval()*mVec[i];
        RVec[i] = randpval()*nVec[i];
        CVec[i] = randpval()*mVec[i]*nVec[i];
        FVec[i] = 0.0;

        assert(mVec[i]>=MINSNPS_B && mVec[i]<=(MINSNPS_B+MAXSNPS_E));
        assert(nVec[i]>=MINSNPS_B && nVec[i]<=(MINSNPS_B+MAXSNPS_E));
        assert(LVec[i]>=0.0f && LVec[i]<=1.0f*mVec[i]);
        assert(RVec[i]>=0.0f && RVec[i]<=1.0f*nVec[i]);
        assert(CVec[i]>=0.0f && CVec[i]<=1.0f*mVec[i]*nVec[i]);
    }
    ////////////////////////////////////////////////////////////


    double timeOmegaTotalStart = gettime();

    for(unsigned int j=0;j<1;j++)
    {
        //printf("%d \n", threadData->threadTOTAL);
        avgF=0;
        maxF=0;
        minF=0;
        updateThread(threadData);
        startThreadOperations(threadData, EXECUTE);
        for (unsigned int k=0;k<threads;k++) {
            maxF = (&threadData[k])->thread_max > maxF ? (&threadData[k])->thread_max : maxF;
            minF = (&threadData[k])->thread_min < minF ? (&threadData[k])->thread_min : minF;
            avgF += (&threadData[k])->thread_sum;

        }


        if (!world_rank) {

            int start = N - k;
            int end = N;
            for (int i = start; i < end; i++) {
                //use scalar (traditional) way to compute remaining of array ( N%4 iterations )
                float num_0 = LVec[i] + RVec[i];
                float num_1 = mVec[i] * (mVec[i] - 1.0f) / 2.0f;
                float num_2 = nVec[i] * (nVec[i] - 1.0f) / 2.0f;
                float num = num_0 / (num_1 + num_2);

                float den_0 = CVec[i] - LVec[i] - RVec[i];
                float den_1 = mVec[i] * nVec[i];
                float den = den_0 / den_1;

                FVec[i] = num / (den + 0.01f);

                maxF = FVec[i] > maxF ? FVec[i] : maxF;
                minF = FVec[i] < minF ? FVec[i] : minF;
                avgF += FVec[i];

            }
        }
    }
    terminateWorkerThreads(workerThread,threadData);
    pthread_barrier_destroy(&barrier);


    MPI_Gather(&maxF, 1, MPI_FLOAT, max, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&minF, 1, MPI_FLOAT, min, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        max = (float *) malloc(world_size * sizeof(float));
        assert(max != NULL);
        min = (float *) malloc(world_size * sizeof(float));
        assert(min != NULL);
    }







    double timeOmegaTotal = gettime()-timeOmegaTotalStart;
    double timeTotalMainStop = gettime();

    // Sync and finalize the MPI environment. No more MPI calls can be made after this
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    if (world_rank == 0) {
        float maxF = 0.0f;
        float minF = FLT_MAX;

        int k=0;
        while(k<world_size){
            maxF = MAX(max[k], maxF);
            minF = MIN(min[k], minF);
            k++;
        }
        printf("Omega time %fs\nTotal time %fs\nMin %e\nMax %e\nAvg %e\n", timeOmegaTotal/iters, timeTotalMainStop-timeTotalMainStart, (double)minF, (double)maxF,(double)avgF/N);
    }

 //   printf("Omega time %fs \n Total time %fs \n Min %e \n Max %e \n Avg %e\n",
   //        timeOmegaTotal/iters, timeTotalMainStop-timeTotalMainStart, (double)minF, (double)maxF,(double)avgF/N);



    _mm_free(mVec);
    _mm_free(nVec);
    _mm_free(LVec);
    _mm_free(RVec);
    _mm_free(CVec);
    _mm_free(FVec);
    pthread_exit(EXIT_SUCCESS);
}