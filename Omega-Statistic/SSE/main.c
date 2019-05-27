#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <float.h>
//ADDED INCLUDE FOR SSE
#include <xmmintrin.h>
#define MINSNPS_B 5
#define MAXSNPS_E 20

double gettime(void);
float randpval (void);

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

float maxCalc(float x, float y){
    if(x > y){
        return x;
    }
    return y;
}

float minCalc(float x, float y){
    if(x < y){
        return x;
    }
    return y;
}

int main(int argc, char ** argv)
{
    assert(argc==2);

    double timeTotalMainStart = gettime();
    float avgF = 0.0f;
    float maxF = 0.0f;
    float minF = FLT_MAX;
    unsigned int N = (unsigned int)atoi(argv[1]);
    unsigned int iters = 10;

    srand(1);

    float * mVec = (float*)_mm_malloc(sizeof(float)*N,16);
    assert(mVec!=NULL);

    float * nVec = (float*)_mm_malloc(sizeof(float)*N,16);
    assert(nVec!=NULL);

    float * LVec = (float*)_mm_malloc(sizeof(float)*N,16);
    assert(LVec!=NULL);

    float * RVec = (float*)_mm_malloc(sizeof(float)*N,16);
    assert(RVec!=NULL);

    float * CVec = (float*)_mm_malloc(sizeof(float)*N,16);
    assert(CVec!=NULL);

    float * FVec = (float*)_mm_malloc(sizeof(float)*N,16);
    assert(FVec!=NULL);


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
        //Modified initial assignment as advised
        assert(LVec[i]>=0.0f && LVec[i]<=1.0f*mVec[i]);
        assert(RVec[i]>=0.0f && RVec[i]<=1.0f*nVec[i]);
        assert(CVec[i]>=0.0f && CVec[i]<=1.0f*mVec[i]*nVec[i]);
    }

    /*
     * SSE CHANGES
     */

    //Implementation with pointers.

    __m128 *mVec_ptr = (__m128 *) mVec;
    __m128 *nVec_ptr = (__m128 *) nVec;
    __m128 *LVec_ptr = (__m128 *) LVec;
    __m128 *RVec_ptr = (__m128 *) RVec;
    __m128 *CVec_ptr = (__m128 *) CVec;
    __m128 *FVec_ptr = (__m128 *) FVec;
    __m128 maxF_vec = _mm_setzero_ps();
    __m128 avgF_vec = _mm_setzero_ps();
    __m128 minF_vec = _mm_set_ps1(FLT_MAX); //TODO CHANGE THAT
    //__m128 *maxv__m128 = (__m128 *) maxv;



    __m128 temp_num_0;
    __m128 temp_num_1;
    __m128 temp_num_2;
    __m128 temp_num;
    __m128 temp_den_0;
    __m128 temp_den_1;
    __m128 temp_den;

    __m128 vec_num;
    __m128 vec_num_0;
    __m128 vec_num_1;
    __m128 vec_num_2;
    __m128 vec_den;
    __m128 vec_den_0;
    __m128 vec_den_1;

    /*
     * We prefer to declare the follow const values as __m128 and
     * not as row data as sse handles better them that way
     */
    const __m128 temp_one = _mm_set1_ps(1.0f);
    const __m128 __temp_one = _mm_set1_ps(0.01f);
    const __m128 temp_two = _mm_set1_ps(2.0f);

    /*
     * List of instructions used as found from Intrinsics Guide. All the instructions
     * we use are based on float as these are demanded
     *
     * -addps - Adds 4 single-precision (32bit) floating-point values to 4 other single-precision floating-point values
     * -subps - Subtracts 4 single-precision floating-point values from 4 other single-precision floating-point values
     * -mulps - Multiplies 4 single-precision floating-point values with 4 other single-precision values.
     * -divps - Divides 4 single-precision floating-point values by 4 other single-precision floating-point values.
     * -maxps - Returns maximum of 2 values in each of 4 single-precision values.
     * -minps - Returns minimum of 2 values in each of 4 single-precision values.
     *
     */
    double timeOmegaTotalStart = gettime();
    for(unsigned int j=0;j<iters;j++)
    {
        avgF = 0.0f;
        maxF = 0.0f;
        minF = FLT_MAX;
        for(unsigned int i=0; i<N/4 ;i++)
        {
            //Replace each statement step by step using the instructions as mentioned before

            //float num_0 = LVec[i] + RVec[i];
            temp_num_0 = _mm_add_ps(LVec_ptr[i],RVec_ptr[i]);

            //float num_1 = mVec[i]*(mVec[i]-1.0f)/2.0f;
            temp_num_1 = _mm_sub_ps(mVec_ptr[i], temp_one);
            temp_num_1 = _mm_mul_ps(mVec_ptr[i], temp_num_1);
            temp_num_1 = _mm_div_ps(temp_num_1, temp_two);


            //float num_2 = nVec[i]*(nVec[i]-1.0f)/2.0f;
            temp_num_2 = _mm_sub_ps(nVec_ptr[i], temp_one);
            temp_num_2 = _mm_mul_ps(nVec_ptr[i], temp_num_2);
            temp_num_2 = _mm_div_ps(temp_num_2, temp_two);

            //float num = num_0/(num_1+num_2);
            temp_num = _mm_add_ps(temp_num_1,temp_num_2);
            temp_num = _mm_div_ps(temp_num_0, temp_num);

            //float den_0 = CVec[i]-LVec[i]-RVec[i];
            temp_den_0 = _mm_sub_ps(CVec_ptr[i], LVec_ptr[i]);
            temp_den_0 = _mm_sub_ps(temp_den_0, RVec_ptr[i]);

            //float den_1 = mVec[i]*nVec[i];
            temp_den_1 = _mm_mul_ps(mVec_ptr[i], nVec_ptr[i]);

            //float den = den_0/den_1;
            temp_den = _mm_div_ps(temp_den_0, temp_den_1);

            //FVec[i] = num/(den+0.01f);
            FVec_ptr[i] = _mm_add_ps(temp_den, __temp_one);
            FVec_ptr[i] = _mm_div_ps(temp_num, FVec_ptr[i]);

            //maxF = FVec[i]>maxF?FVec[i]:maxF;
            maxF_vec = _mm_max_ps(FVec_ptr[i], maxF_vec);

            //minF = FVec[i]<minF?FVec[i]:minF;
            minF_vec = _mm_min_ps(FVec_ptr[i], minF_vec);

            //avgF += FVec[i];
            avgF_vec = _mm_add_ps(avgF_vec, FVec_ptr[i]);

        }

        //_mm_store_ps(max, _maxF);
        //_mm_store_ps(max, maxFV);
        maxF = maxF_vec[0];
        maxF = maxCalc(maxF_vec[1],maxF);
        maxF = maxCalc(maxF_vec[2],maxF);
        maxF = maxCalc(maxF_vec[3],maxF);

        minF = minF_vec[0];
        minF = minCalc(minF_vec[1],minF);
        minF = minCalc(minF_vec[2],minF);
        minF = minCalc(minF_vec[3],minF);

        avgF = avgF_vec[0]+avgF_vec[1]+avgF_vec[2]+avgF_vec[3];

        for (int j = (N - N % 4); j < N; j++) {
            // Scalar computation of whatever its left
            float num_0 = LVec[j] + RVec[j];
            float num_1 = mVec[j] * (mVec[j] - 1.0f) / 2.0f;
            float num_2 = nVec[j] * (nVec[j] - 1.0f) / 2.0f;
            float num = num_0 / (num_1 + num_2);
            float den_0 = CVec[j] - LVec[j] - RVec[j];
            float den_1 = mVec[j] * nVec[j];
            float den = den_0 / den_1;

            FVec[j] = num / (den + 0.01f);
            maxF = FVec[j] > maxF ? FVec[j] : maxF;
            minF = FVec[j]<minF?FVec[j]:minF;
            avgF += FVec[j];
        }

    }
    double timeOmegaTotal = gettime()-timeOmegaTotalStart;
    double timeTotalMainStop = gettime();

    printf("Omega time %fs - Total time %fs - Min %e - Max %e - Avg %e\n",
           timeOmegaTotal/iters, timeTotalMainStop-timeTotalMainStart, (double)minF, (double)maxF,
           (double)avgF/N);

    free(mVec);
    free(nVec);
    free(LVec);
    free(RVec);
    free(CVec);
    free(FVec);
}

//
//#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>
//#include <sys/time.h>
//#include <assert.h>
//#include <float.h>
////ADDED INCLUDE FOR SSE
//#include <xmmintrin.h>
//#define MINSNPS_B 5
//#define MAXSNPS_E 20
//
//double gettime(void);
//float randpval (void);
//
//double gettime(void)
//{
//    struct timeval ttime;
//    gettimeofday(&ttime , NULL);
//    return ttime.tv_sec + ttime.tv_usec * 0.000001;
//}
//
//float randpval (void)
//{
//    int vr = rand();
//    int vm = rand()%vr;
//    float r = ((float)vm)/(float)vr;
//    assert(r>=0.0f && r<=1.00001f);
//    return r;
//}
//
//float maxCalc(float x, float y){
//    if(x > y){
//        return x;
//    }
//    return y;
//}
//
//float minCalc(float x, float y){
//    if(x < y){
//        return x;
//    }
//    return y;
//}
//
//int main(int argc, char ** argv)
//{
//    assert(argc==2);
//
//    double timeTotalMainStart = gettime();
//    float avgF = 0.0f;
//    float maxF = 0.0f;
//    float minF = FLT_MAX;
//    unsigned int N = (unsigned int)atoi(argv[1]);
//    unsigned int iters = 10;
//
//    srand(1);
//
//    float * mVec = (float*)_mm_malloc(sizeof(float)*N,16);
//    assert(mVec!=NULL);
//
//    float * nVec = (float*)_mm_malloc(sizeof(float)*N,16);
//    assert(nVec!=NULL);
//
//    float * LVec = (float*)_mm_malloc(sizeof(float)*N,16);
//    assert(LVec!=NULL);
//
//    float * RVec = (float*)_mm_malloc(sizeof(float)*N,16);
//    assert(RVec!=NULL);
//
//    float * CVec = (float*)_mm_malloc(sizeof(float)*N,16);
//    assert(CVec!=NULL);
//
//    float * FVec = (float*)_mm_malloc(sizeof(float)*N,16);
//    assert(FVec!=NULL);
//
//
//    for(unsigned int i=0;i<N;i++)
//    {
//        mVec[i] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
//        nVec[i] = (float)(MINSNPS_B+rand()%MAXSNPS_E);
//        LVec[i] = randpval()*mVec[i];
//        RVec[i] = randpval()*nVec[i];
//        CVec[i] = randpval()*mVec[i]*nVec[i];
//        FVec[i] = 0.0;
//
//        assert(mVec[i]>=MINSNPS_B && mVec[i]<=(MINSNPS_B+MAXSNPS_E));
//        assert(nVec[i]>=MINSNPS_B && nVec[i]<=(MINSNPS_B+MAXSNPS_E));
//        //Modified initial assignment as advised
//        assert(LVec[i]>=0.0f && LVec[i]<=1.0f*mVec[i]);
//        assert(RVec[i]>=0.0f && RVec[i]<=1.0f*nVec[i]);
//        assert(CVec[i]>=0.0f && CVec[i]<=1.0f*mVec[i]*nVec[i]);
//    }
//
//    /*
//     * SSE CHANGES
//     */
//
//    //Implementation with pointers.
//
//    __m128 *mVec_ptr = (__m128 *) mVec;
//    __m128 *nVec_ptr = (__m128 *) nVec;
//    __m128 *LVec_ptr = (__m128 *) LVec;
//    __m128 *RVec_ptr = (__m128 *) RVec;
//    __m128 *CVec_ptr = (__m128 *) CVec;
//    __m128 *FVec_ptr = (__m128 *) FVec;
//    __m128 maxF_vec = _mm_setzero_ps();
//    __m128 avgF_vec = _mm_setzero_ps();
//    __m128 avgF_vec_1 = _mm_setzero_ps();
//    __m128 minF_vec = _mm_set_ps1(FLT_MAX); //TODO CHANGE THAT
//    //__m128 *maxv__m128 = (__m128 *) maxv;
//
//
//
//    __m128 temp_num_0;
//    __m128 temp_num_0_1;
//    __m128 _num_0;
//    __m128 temp_num_1;
//    __m128 temp_num_1_1;
//    __m128 _num_1;
//    __m128 temp_num_2;
//    __m128 temp_num_2_1;
//    __m128 _num_2;
//    __m128 temp_num;
//    __m128 temp_num_;
//    __m128 _num;
//    __m128 temp_den_0;
//    __m128 temp_den_0_1;
//    __m128 _den_0;
//    __m128 temp_den_1;
//    __m128 temp_den_1_1;
//    __m128 _den_1;
//    __m128 temp_den;
//    __m128 temp_den_;
//    __m128 _den;
//
//    __m128 vec_num;
//    __m128 vec_num_0;
//    __m128 vec_num_1;
//    __m128 vec_num_2;
//    __m128 vec_den;
//    __m128 vec_den_0;
//    __m128 vec_den_1;
//
//    /*
//     * We prefer to declare the follow const values as __m128 and
//     * not as row data as sse handles better them that way
//     */
//    const __m128 temp_one = _mm_set1_ps(1.0f);
//    const __m128 __temp_one = _mm_set1_ps(0.01f);
//    const __m128 temp_two = _mm_set1_ps(2.0f);
//
//    /*
//     * List of instructions used as found from Intrinsics Guide. All the instructions
//     * we use are based on float as these are demanded
//     *
//     * -addps - Adds 4 single-precision (32bit) floating-point values to 4 other single-precision floating-point values
//     * -subps - Subtracts 4 single-precision floating-point values from 4 other single-precision floating-point values
//     * -mulps - Multiplies 4 single-precision floating-point values with 4 other single-precision values.
//     * -divps - Divides 4 single-precision floating-point values by 4 other single-precision floating-point values.
//     * -maxps - Returns maximum of 2 values in each of 4 single-precision values.
//     * -minps - Returns minimum of 2 values in each of 4 single-precision values.
//     *
//     */
//    double timeOmegaTotalStart = gettime();
//    for(unsigned int j=0;j<iters;j++)
//    {
//        avgF = 0.0f;
//        maxF = 0.0f;
//        minF = FLT_MAX;
//        for(unsigned int i=0; i<N/8 ;i++)
//        {
//            //Replace each statement step by step using the instructions as mentioned before
//
//            //float num_0 = LVec[i] + RVec[i];
//            temp_num_0 = _mm_add_ps(LVec_ptr[i],RVec_ptr[i]);
//            temp_num_0_1 = _mm_add_ps(LVec_ptr[i+4],RVec_ptr[i+4]);
//
//
//            //float num_1 = mVec[i]*(mVec[i]-1.0f)/2.0f;
//            temp_num_1 = _mm_sub_ps(mVec_ptr[i], temp_one);
//            temp_num_1 = _mm_mul_ps(mVec_ptr[i], temp_num_1);
//            temp_num_1 = _mm_div_ps(temp_num_1, temp_two);
//
//            temp_num_1_1 = _mm_sub_ps(mVec_ptr[i+4], temp_one);
//            temp_num_1_1 = _mm_mul_ps(mVec_ptr[i+4], temp_num_1_1);
//            temp_num_1_1 = _mm_div_ps(temp_num_1_1, temp_two);
//
//
//            //float num_2 = nVec[i]*(nVec[i]-1.0f)/2.0f;
//            temp_num_2 = _mm_sub_ps(nVec_ptr[i], temp_one);
//            temp_num_2 = _mm_mul_ps(nVec_ptr[i], temp_num_2);
//            temp_num_2 = _mm_div_ps(temp_num_2, temp_two);
//
//            temp_num_2_1 = _mm_sub_ps(nVec_ptr[i+4], temp_one);
//            temp_num_2_1 = _mm_mul_ps(nVec_ptr[i+4], temp_num_2_1);
//            temp_num_2_1 = _mm_div_ps(temp_num_2_1, temp_two);
//
//            //float num = num_0/(num_1+num_2);
//            temp_num = _mm_add_ps(temp_num_1,temp_num_2);
//            temp_num = _mm_div_ps(temp_num_0, temp_num);
//
//            temp_num_ = _mm_add_ps(temp_num_1_1,temp_num_2_1);
//            temp_num_ = _mm_div_ps(temp_num_0_1, temp_num_);
//
//            //float den_0 = CVec[i]-LVec[i]-RVec[i];
//            temp_den_0 = _mm_sub_ps(CVec_ptr[i], LVec_ptr[i]);
//            temp_den_0 = _mm_sub_ps(temp_den_0, RVec_ptr[i]);
//
//            temp_den_0_1 = _mm_sub_ps(CVec_ptr[i+4], LVec_ptr[i+4]);
//            temp_den_0_1 = _mm_sub_ps(temp_den_0_1, RVec_ptr[i+4]);
//
//            //float den_1 = mVec[i]*nVec[i];
//            temp_den_1 = _mm_mul_ps(mVec_ptr[i], nVec_ptr[i]);
//
//            temp_den_1_1 = _mm_mul_ps(mVec_ptr[i+4], nVec_ptr[i+4]);
//
//            //float den = den_0/den_1;
//            temp_den = _mm_div_ps(temp_den_0, temp_den_1);
//
//            temp_den_ = _mm_div_ps(temp_den_0_1, temp_den_1_1);
//
//            //FVec[i] = num/(den+0.01f);
//            FVec_ptr[i] = _mm_add_ps(temp_den, __temp_one);
//            FVec_ptr[i] = _mm_div_ps(temp_num, FVec_ptr[i]);
//
//            FVec_ptr[i+4] = _mm_add_ps(temp_den_, __temp_one);
//            FVec_ptr[i+4] = _mm_div_ps(temp_num_, FVec_ptr[i+4]);
//
//            //maxF = FVec[i]>maxF?FVec[i]:maxF;
//            maxF_vec = _mm_max_ps(FVec_ptr[i], maxF_vec);
//
//            //minF = FVec[i]<minF?FVec[i]:minF;
//            minF_vec = _mm_min_ps(FVec_ptr[i], minF_vec);
//
//            //avgF += FVec[i];
//            avgF_vec = _mm_add_ps(avgF_vec, FVec_ptr[i]);
//            avgF_vec_1 = _mm_add_ps(avgF_vec, FVec_ptr[i+4]);
//
//        }
//
//        //_mm_store_ps(max, _maxF);
//        //_mm_store_ps(max, maxFV);
////        maxF = maxF_vec[0];
////        maxF = maxCalc(maxF_vec[1],maxF);
////        maxF = maxCalc(maxF_vec[2],maxF);
////        maxF = maxCalc(maxF_vec[3],maxF);
////
////        minF = minF_vec[0];
////        minF = minCalc(minF_vec[1],minF);
////        minF = minCalc(minF_vec[2],minF);
////        minF = minCalc(minF_vec[3],minF);
//
//        avgF = avgF_vec[0]+avgF_vec[1]+avgF_vec[2]+avgF_vec[3]+avgF_vec_1[0]+avgF_vec_1[1]+avgF_vec_1[2]+avgF_vec_1[3];
//
//        for (int j = (N - N % 8); j < N; j++) {
//            // Scalar computation of whatever its left
//            float num_0 = LVec[j] + RVec[j];
//            float num_1 = mVec[j] * (mVec[j] - 1.0f) / 2.0f;
//            float num_2 = nVec[j] * (nVec[j] - 1.0f) / 2.0f;
//            float num = num_0 / (num_1 + num_2);
//            float den_0 = CVec[j] - LVec[j] - RVec[j];
//            float den_1 = mVec[j] * nVec[j];
//            float den = den_0 / den_1;
//
//            FVec[j] = num / (den + 0.01f);
//            maxF = FVec[j] > maxF ? FVec[j] : maxF;
//            minF = FVec[j]<minF?FVec[j]:minF;
//            avgF += FVec[j];
//        }
//
//    }
//    double timeOmegaTotal = gettime()-timeOmegaTotalStart;
//    double timeTotalMainStop = gettime();
//
//    printf("Omega time %fs - Total time %fs - Min %e - Max %e - Avg %e\n",
//           timeOmegaTotal/iters, timeTotalMainStop-timeTotalMainStart, (double)minF, (double)maxF,
//           (double)avgF/N);
//
//    free(mVec);
//    free(nVec);
//    free(LVec);
//    free(RVec);
//    free(CVec);
//    free(FVec);
//}