/*
            ======================================================
            + Name        : bonus.c 		                         +
            + Author      : Christos Z, 			             +
            + Description : SSE extension					     +
            ======================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <float.h>
#include <xmmintrin.h>
#include <pmmintrin.h>

#define MINSNPS_B 5
#define MAXSNPS_E 20

double gettime(void);
float randpval (void);

double gettime(void){
    struct timeval ttime;
    gettimeofday(&ttime , NULL);
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

float randpval (void){
    int vr = rand();
    int vm = rand()%vr;
    float r = ((float)vm)/(float)vr;
    assert(r>=0.0f && r<=1.00001f);
    return r;
}

int main(int argc,char** argv){

    assert(argc==2);
    double timeTotalMainStart = gettime();

    unsigned int N = (unsigned int)atoi(argv[1]);
    unsigned int iters = 10;

    unsigned int leftovers=N%4;

    srand(1);

    float avgF = 0.0f;
    float maxF = FLT_MIN;
    float minF = FLT_MAX;


    float * maxVec = (float*)_mm_malloc(sizeof(float),32);
    assert(maxVec!=NULL);
    float * minVec = (float*)_mm_malloc(sizeof(float),32);
    assert(minVec!=NULL);
    float * avgVec = (float*)_mm_malloc(sizeof(float),32);
    assert(avgVec!=NULL);


    for(int i=0;i<4;i++){
        maxVec[i]=FLT_MIN;
        minVec[i]=FLT_MAX;
        avgVec[i]=avgF;
    }

    // mVec,nVec,LVec,RVec,CVec,FVec
    float * allVec = (float*)_mm_malloc(sizeof(float)*N*6,32);
    assert(allVec!=NULL);



/*	float * mVec = (float*)_mm_malloc(sizeof(float)*N,32);
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
	assert(FVec!=NULL); */

    //Initialize mVec,,nVec,LVec,RVec,CVec,FVec

    double timeOmegaTotalStart = gettime();
	
    for(unsigned int i=0;i<(N*6);i+=24){
        allVec[i] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
        allVec[i+4] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
        allVec[i+8] = randpval()*allVec[i]*allVec[i+4];
        allVec[i+12] = randpval()*allVec[i];
        allVec[i+16] = randpval()*allVec[i+4];
        allVec[i+20] = 0.0;

        allVec[i+1] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
        allVec[i+5] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
        allVec[i+9] = randpval()*allVec[i+1]*allVec[i+5];
        allVec[i+13] = randpval()*allVec[i+1];
        allVec[i+17] = randpval()*allVec[i+5];
        allVec[i+21] = 0.0;


        allVec[i+2] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
        allVec[i+6] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
        allVec[i+10] = randpval()*allVec[i+2]*allVec[i+6];
        allVec[i+14] = randpval()*allVec[i+2];
        allVec[i+18] = randpval()*allVec[i+6];
        allVec[i+22] = 0.0;

        allVec[i+3] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
        allVec[i+7] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
        allVec[i+11] = randpval()*allVec[i+3]*allVec[i+7];
        allVec[i+15] = randpval()*allVec[i+3];
        allVec[i+19] = randpval()*allVec[i+7];
        allVec[i+23] = 0.0;

    }

    //Vectors with constants.;

    //Vector all data
    __m128 * allFVec128 = (__m128 *) allVec;

    //	Vectors for max,min,avg.
	__m128 maxFVec128 = _mm_set_ps1(FLT_MIN);
	__m128 minFVec128 = _mm_set_ps1(FLT_MAX);
	__m128 avgFVec128 = _mm_set_ps1(0);


    __m128 num_0;
    __m128 num_1;
    __m128 num_2;
    __m128 num;
    __m128 den_0;
    __m128 den_1;
    __m128 den;
    
        __m128 vec1 = _mm_set1_ps(1.0f);
    __m128 vec2 = _mm_set1_ps(2.0f);
    __m128 vec3 = _mm_set_ps1(0.01f);

	__m128 comparemaxFlags;
	__m128 compareminFlags;



	for(unsigned int j=0;j<iters;j++){
		avgF = 0.0f;
		maxF = 0.0f;
		minF = FLT_MAX;
		
		for(unsigned int i=0;i<(N*6)/4;i+=6){
			num_0 = _mm_add_ps(allFVec128[i+3], allFVec128[i+4]);
			
			num_1 = _mm_div_ps( _mm_mul_ps( allFVec128[i], _mm_sub_ps( allFVec128[i], vec1)), vec2);
			num_2 = _mm_div_ps( _mm_mul_ps( allFVec128[i+1], _mm_sub_ps( allFVec128[i+1], vec1)), vec2);
			
			num = _mm_div_ps( num_0, _mm_add_ps( num_1, num_2));
			
			den_0 = _mm_sub_ps( _mm_sub_ps( allFVec128[i+2], allFVec128[i+3]), allFVec128[i+4]);
			den_1 = _mm_mul_ps( allFVec128[i], allFVec128[i+1]);
			den = _mm_div_ps( den_0, den_1);
			
			allFVec128[i+5] = _mm_div_ps(num, _mm_add_ps(den, vec3));


			comparemaxFlags = _mm_cmpgt_ps( allFVec128[i+5], maxFVec128);
			compareminFlags=  _mm_cmplt_ps( allFVec128[i+5], minFVec128);
			maxFVec128 = _mm_or_ps(_mm_and_ps( comparemaxFlags, allFVec128[i+5]),  _mm_andnot_ps(comparemaxFlags, maxFVec128));
			minFVec128 = _mm_or_ps(_mm_and_ps( compareminFlags, allFVec128[i+5]),  _mm_andnot_ps(compareminFlags, minFVec128));
			avgFVec128 = _mm_add_ps(avgFVec128,allFVec128[i+5]);
		//	maxFVec128= _mm_max_ps(*maxFVec128,allFVec128[i+5]);
		//	minFVec128= _mm_min_ps(allFVec128[i+5],*minFVec128);
	//		avgFVec128= _mm_add_ps(allFVec128[i+5],*avgFVec128);
			
			
		
		} 
	/*	maxF = maxVec[0];
   		maxF = maxVec[1] > maxF ? maxVec[1] : maxF;
   		maxF = maxVec[2] > maxF ? maxVec[2] : maxF;
   		maxF = maxVec[3] > maxF ? maxVec[3] : maxF;

   		minF = minVec[0];
   		minF = minVec[1] < minF ? minVec[1] : minF;
   		minF = minVec[2] < minF ? minVec[2] : minF;
   		minF = minVec[3] < minF ? minVec[3] : minF;

   		avgF = avgVec[0] + avgVec[1] + avgVec[2] + avgVec[3]; 	*/

			//	Horizontal max.
		__m128 tempMax = _mm_max_ps(maxFVec128, _mm_shuffle_ps(maxFVec128, maxFVec128, _MM_SHUFFLE(0,0,3,2)));
		maxF = _mm_cvtss_f32(_mm_max_ps(tempMax, _mm_shuffle_ps(tempMax, tempMax, _MM_SHUFFLE(0,0,0,1))));
		__m128 tempMin = _mm_min_ps(minFVec128, _mm_shuffle_ps(minFVec128, minFVec128, _MM_SHUFFLE(0,0,3,2)));
		minF = _mm_cvtss_f32(_mm_min_ps(tempMin, _mm_shuffle_ps(tempMin, tempMin, _MM_SHUFFLE(0,0,0,1))));
		__m128 tempAvg = _mm_add_ps(avgFVec128, _mm_shuffle_ps(avgFVec128, avgFVec128, _MM_SHUFFLE(0,0,3,2)));
		avgF = _mm_cvtss_f32(_mm_add_ps(tempAvg, _mm_shuffle_ps(tempAvg, tempAvg, _MM_SHUFFLE(0,0,0,1))));


			  //Leftovers;
	   	for(unsigned int i=N-leftovers;i<N;i+=6)
		{
				float num_0 = allVec[i+3]+allVec[i+4];
				float num_1 = allVec[i]*(allVec[i]-1.0f)/2.0f;
				float num_2 = allVec[i+1]*(allVec[i+1]-1.0f)/2.0f;
				float num = num_0/(num_1+num_2);
	
				float den_0 = allVec[i+2]-allVec[i+3]-allVec[i+4];
				float den_1 = allVec[i]*allVec[i+1];
				float den = den_0/den_1;
	
				allVec[i+5] = num/(den+0.01f);
				
				maxF = allVec[i+5]>maxF?allVec[i+5]:maxF;
				minF = allVec[i+5]<minF?allVec[i+5]:minF;
				avgF += allVec[i+5];
		}	

		
	}








    double timeOmegaTotal = gettime()-timeOmegaTotalStart;
    double timeTotalMainStop = gettime();
    printf("Omega time %fs\nTotal time %fs\nMin %e\nMax %e\nAvg %e\n", timeOmegaTotal/iters, timeTotalMainStop-timeTotalMainStart, (double)minF, (double)maxF,(double)avgF/N);
    _mm_free(allVec);

}

