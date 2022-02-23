/*
            ======================================================
            + Name        : SSE.c 		                         +
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
	float avgF = 0.0f;
	float maxF = FLT_MIN;
	float minF = FLT_MAX;
	unsigned int N = (unsigned int)atoi(argv[1]);
	unsigned int iters = 10;
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

	for(unsigned int i=0;i<N;i++){
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

	//Vector pointers to float arrays.
	__m128* LVec128 = (__m128*) LVec;
	__m128* RVec128 = (__m128*) RVec;
	__m128* mVec128 = (__m128*) mVec;
	__m128* nVec128 = (__m128*) nVec;
	__m128* CVec128 = (__m128*) CVec;
	__m128* FVec128 = (__m128*) FVec;

	double timeOmegaTotalStart = gettime();

	for(unsigned int j=0;j<iters;j++){
		__m128 avgFVec128 = _mm_set_ps1(0);
		unsigned int i;
		for(i=0; i < N/4; i++){
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
			//minFVec128 = _mm_min_ps(minFVec128,FVec128[i]);
			//float result[4];
			//_mm_store_ps(result, FVec128[i]);
			// avg=avg+(result[0]+result[1]+ result[2]+ result[3]);

		}
		//	Horizontal max. Look bottom for a more detailed explanation.
		__m128 tempMax = _mm_max_ps(maxFVec128, _mm_shuffle_ps(maxFVec128, maxFVec128, _MM_SHUFFLE(0,0,3,2)));
		maxF = _mm_cvtss_f32(_mm_max_ps(tempMax, _mm_shuffle_ps(tempMax, tempMax, _MM_SHUFFLE(0,0,0,1))));
		__m128 tempMin = _mm_min_ps(minFVec128, _mm_shuffle_ps(minFVec128, minFVec128, _MM_SHUFFLE(0,0,3,2)));
		minF = _mm_cvtss_f32(_mm_min_ps(tempMin, _mm_shuffle_ps(tempMin, tempMin, _MM_SHUFFLE(0,0,0,1))));
		__m128 tempAvg = _mm_add_ps(avgFVec128, _mm_shuffle_ps(avgFVec128, avgFVec128, _MM_SHUFFLE(0,0,3,2)));
		avgF = _mm_cvtss_f32(_mm_add_ps(tempAvg, _mm_shuffle_ps(tempAvg, tempAvg, _MM_SHUFFLE(0,0,0,1))));


		for(i=i*4;i<N;i++)
		{
			float num_0 = LVec[i]+RVec[i];
			float num_1 = mVec[i]*(mVec[i]-1.0f)/2.0f;
			float num_2 = nVec[i]*(nVec[i]-1.0f)/2.0f;
			float num = num_0/(num_1+num_2);

			float den_0 = CVec[i]-LVec[i]-RVec[i];
			float den_1 = mVec[i]*nVec[i];
			float den = den_0/den_1;
			
			FVec[i] = num/(den+0.01f);
			
			maxF = FVec[i]>maxF?FVec[i]:maxF;
			minF = FVec[i]<minF?FVec[i]:minF;
			avgF += FVec[i];
		}
	}
	double timeOmegaTotal = gettime()-timeOmegaTotalStart;
	double timeTotalMainStop = gettime();
	printf("Omega time %fs\nTotal time %fs\nMin %e\nMax %e\nAvg %e\n", timeOmegaTotal/iters, timeTotalMainStop-timeTotalMainStart, (double)minF, (double)maxF,(double)avgF/N);
	_mm_free(mVec128);
	_mm_free(nVec128);
	_mm_free(LVec128);
	_mm_free(RVec128);
	_mm_free(CVec128);
	_mm_free(FVec128);
}
