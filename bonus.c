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
	__m128 comparemaxFlags[4];
	__m128 compareminFlags[4];

	//Vector pointers to float arrays.
	__m128* LVec128 = (__m128*) LVec;
	__m128* RVec128 = (__m128*) RVec;
	__m128* mVec128 = (__m128*) mVec;
	__m128* nVec128 = (__m128*) nVec;
	__m128* CVec128 = (__m128*) CVec;
	__m128* FVec128 = (__m128*) FVec;

	__m128 num_0[4];
	__m128 num_1[4];
	__m128 num_2[4];
	__m128 num[4];
	__m128 den_0[4];
	__m128 den_1[4];
	__m128 den[4];

	double timeOmegaTotalStart = gettime();

	for(unsigned int j=0;j<iters;j++){
		__m128 avgFVec128 = _mm_set_ps1(0);
		maxF=0;
		minF=0;
		avgF=0;
		unsigned int i;
		for(i=0; i < N/16; i++){
			//	Same calculations with the serial execution. Floating point operations are non-associative and their order is respected.
			num_0[0] = _mm_add_ps( LVec128[i], RVec128[i]);
			num_0[1] = _mm_add_ps( LVec128[(N/16)+i], RVec128[(N/16)+i]);
			num_0[2] = _mm_add_ps( LVec128[(2*N/16)+i], RVec128[(2*N/16)+i]);
			num_0[3] = _mm_add_ps( LVec128[(3*N/16)+i], RVec128[(3*N/16)+i]);

			num_1[0] = _mm_div_ps( _mm_mul_ps( mVec128[i], _mm_sub_ps( mVec128[i], vec1)), vec2);
			num_1[1] = _mm_div_ps( _mm_mul_ps( mVec128[(N/16)+i], _mm_sub_ps( mVec128[(N/16)+i], vec1)), vec2);
			num_1[2] = _mm_div_ps( _mm_mul_ps( mVec128[(2*N/16)+i], _mm_sub_ps( mVec128[(2*N/16)+i], vec1)), vec2);
			num_1[3] = _mm_div_ps( _mm_mul_ps( mVec128[(3*N/16)+i], _mm_sub_ps( mVec128[(3*N/16)+i], vec1)), vec2);

			num_2[0] = _mm_div_ps( _mm_mul_ps( nVec128[i], _mm_sub_ps( nVec128[i], vec1)), vec2);
			num_2[1] = _mm_div_ps( _mm_mul_ps( nVec128[(N/16)+i], _mm_sub_ps( nVec128[(N/16)+i], vec1)), vec2);
			num_2[2] = _mm_div_ps( _mm_mul_ps( nVec128[(2*N/16)+i], _mm_sub_ps( nVec128[(2*N/16)+i], vec1)), vec2);
			num_2[3] = _mm_div_ps( _mm_mul_ps( nVec128[(3*N/16)+i], _mm_sub_ps( nVec128[(3*N/16)+i], vec1)), vec2);


			num[0] = _mm_div_ps( num_0[0], _mm_add_ps( num_1[0], num_2[0]));
			num[1] = _mm_div_ps( num_0[1], _mm_add_ps( num_1[1], num_2[1]));
			num[2] = _mm_div_ps( num_0[2], _mm_add_ps( num_1[2], num_2[2]));
			num[3] = _mm_div_ps( num_0[3], _mm_add_ps( num_1[3], num_2[3]));

			den_0[0] = _mm_sub_ps( _mm_sub_ps( CVec128[i], LVec128[i]), RVec128[i]);
			den_0[1] = _mm_sub_ps( _mm_sub_ps( CVec128[(N/16)+i], LVec128[(N/16)+i]), RVec128[(N/16)+i]);
			den_0[2] = _mm_sub_ps( _mm_sub_ps( CVec128[(2*N/16)+i], LVec128[(2*N/16)+i]), RVec128[(2*N/16)+i]);
			den_0[3] = _mm_sub_ps( _mm_sub_ps( CVec128[(3*N/16)+i], LVec128[(3*N/16)+i]), RVec128[(3*N/16)+i]);


			den_1[0] = _mm_mul_ps( mVec128[i], nVec128[i]);
			den_1[1] = _mm_mul_ps( mVec128[(N/16)+i], nVec128[(N/16)+i]);
			den_1[2] = _mm_mul_ps( mVec128[(2*N/16)+i], nVec128[(2*N/16)+i]);
			den_1[3] = _mm_mul_ps( mVec128[(3*N/16)+i], nVec128[(3*N/16)+i]);

			den[0] = _mm_div_ps( den_0[0], den_1[0]);
			den[1] = _mm_div_ps( den_0[1], den_1[1]);
			den[2] = _mm_div_ps( den_0[2], den_1[2]);
			den[3] = _mm_div_ps( den_0[3], den_1[3]);

			FVec128[i] = _mm_div_ps(num[0], _mm_add_ps(den[0], vec001));
			FVec128[(N/16)+i] = _mm_div_ps(num[1], _mm_add_ps(den[1], vec001));
			FVec128[(2*N/16)+i] = _mm_div_ps(num[2], _mm_add_ps(den[2], vec001));
			FVec128[(3*N/16)+i] = _mm_div_ps(num[3], _mm_add_ps(den[3], vec001));

			comparemaxFlags[0] = _mm_cmpgt_ps( FVec128[i], maxFVec128);
			comparemaxFlags[1] = _mm_cmpgt_ps( FVec128[(N/16)+i], maxFVec128);
			comparemaxFlags[2] = _mm_cmpgt_ps( FVec128[(2*N/16)+i], maxFVec128);
			comparemaxFlags[3] = _mm_cmpgt_ps( FVec128[(3*N/16)+i], maxFVec128);

			compareminFlags[0] =  _mm_cmplt_ps( FVec128[i], minFVec128);
			compareminFlags[1] =  _mm_cmplt_ps( FVec128[(N/16)+i], minFVec128);
			compareminFlags[2]=  _mm_cmplt_ps( FVec128[(2*N/16)+i], minFVec128);
			compareminFlags[3]=  _mm_cmplt_ps( FVec128[(3*N/16)+i], minFVec128);

			maxFVec128 = _mm_or_ps(_mm_and_ps( comparemaxFlags[0], FVec128[i]),  _mm_andnot_ps(comparemaxFlags[0], maxFVec128));
			maxFVec128 = _mm_or_ps(_mm_and_ps( comparemaxFlags[1], FVec128[(N/16)+i]),  _mm_andnot_ps(comparemaxFlags[1], maxFVec128));
			maxFVec128 = _mm_or_ps(_mm_and_ps( comparemaxFlags[2], FVec128[(2*N/16)+i]),  _mm_andnot_ps(comparemaxFlags[2], maxFVec128));
			maxFVec128 = _mm_or_ps(_mm_and_ps( comparemaxFlags[3], FVec128[(3*N/16)+i]),  _mm_andnot_ps(comparemaxFlags[3], maxFVec128));

			minFVec128 = _mm_or_ps(_mm_and_ps( compareminFlags[0], FVec128[i]),  _mm_andnot_ps(compareminFlags[0], minFVec128));
			minFVec128 = _mm_or_ps(_mm_and_ps( compareminFlags[1], FVec128[(N/16)+i]),  _mm_andnot_ps(compareminFlags[1], minFVec128));
			minFVec128 = _mm_or_ps(_mm_and_ps( compareminFlags[2], FVec128[(2*N/16)+i]),  _mm_andnot_ps(compareminFlags[2], minFVec128));
			minFVec128 = _mm_or_ps(_mm_and_ps( compareminFlags[3], FVec128[(3*N/16)+i]),  _mm_andnot_ps(compareminFlags[3], minFVec128));

			avgFVec128 = _mm_add_ps(avgFVec128,FVec128[i]);			
			avgFVec128 = _mm_add_ps(avgFVec128,FVec128[(N/16)+i]);
			avgFVec128 = _mm_add_ps(avgFVec128,FVec128[(2*N/16)+i]);
			avgFVec128 = _mm_add_ps(avgFVec128,FVec128[(3*N/16)+i]);

		}
		//	Horizontal max. Look bottom for a more detailed explanation.
		__m128 tempMax = _mm_max_ps(maxFVec128, _mm_shuffle_ps(maxFVec128, maxFVec128, _MM_SHUFFLE(0,0,3,2)));
		maxF = _mm_cvtss_f32(_mm_max_ps(tempMax, _mm_shuffle_ps(tempMax, tempMax, _MM_SHUFFLE(0,0,0,1))));
		__m128 tempMin = _mm_min_ps(minFVec128, _mm_shuffle_ps(minFVec128, minFVec128, _MM_SHUFFLE(0,0,3,2)));
		minF = _mm_cvtss_f32(_mm_min_ps(tempMin, _mm_shuffle_ps(tempMin, tempMin, _MM_SHUFFLE(0,0,0,1))));
		__m128 tempAvg = _mm_add_ps(avgFVec128, _mm_shuffle_ps(avgFVec128, avgFVec128, _MM_SHUFFLE(0,0,3,2)));
		avgF = _mm_cvtss_f32(_mm_add_ps(tempAvg, _mm_shuffle_ps(tempAvg, tempAvg, _MM_SHUFFLE(0,0,0,1))));


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
