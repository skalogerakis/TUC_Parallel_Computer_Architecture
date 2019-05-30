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

    /*
     * Three malloc used in this memory layout. We merged arrays [L,R] , [m,n] , [F,C] from the initial code
     */
    float * LRVec = (float*)_mm_malloc(sizeof(float)*2*N,16);
    assert(LRVec!=NULL);

    float * mnVec = (float*)_mm_malloc(sizeof(float)*2*N,16);
    assert(mnVec!=NULL);

    float * FCVec = (float*)_mm_malloc(sizeof(float)*2*N,16);
    assert(FCVec!=NULL);

    /*
     * Initialize all variables. Initialization changed for different layout implementation
     */
    unsigned int i=0;
    while(i < 2 * N){

        /*
         * When we reach each time the fourth element jump 8 positions so that you do not overwrite anything else
         */
        if((i+4)%8 == 0){
            i = i+4;
            continue;
        }else{
            mnVec[i] = (float)(MINSNPS_B+rand()%MAXSNPS_E); //mVec
            mnVec[i+4] = (float)(MINSNPS_B+rand()%MAXSNPS_E);   //nVec

            LRVec[i] = randpval()*mnVec[i];  //LVec
            LRVec[i+4] = randpval()*mnVec[i+4]; //RVec

            FCVec[i+4] = randpval()*mnVec[i]*mnVec[i+4];  //CVec
            FCVec[i] = 0.0; //Fvec

            assert(mnVec[i]>=MINSNPS_B &&mnVec[i]<=(MINSNPS_B+MAXSNPS_E));
            assert(mnVec[i+4]>=MINSNPS_B && mnVec[i+4]<=(MINSNPS_B+MAXSNPS_E));
            //Modified initial assignment as advised
            assert(LRVec[i]>=0.0f && LRVec[i]<=1.0f*mnVec[i]);
            assert(LRVec[i+4]>=0.0f && LRVec[i+4]<=1.0f*mnVec[i+4]);
            assert(FCVec[i+4]>=0.0f && FCVec[i+4]<=1.0f*mnVec[i]*mnVec[i+4]);
        }
        i++;
    }

    /*
     * SSE CHANGES
     */

    //Implementation with pointers.
    __m128 *mnVec_ptr = (__m128 *) mnVec;
    __m128 *LRVec_ptr = (__m128 *) LRVec;
    __m128 *FCVec_ptr = (__m128 *) FCVec;

    __m128 maxF_vec = _mm_setzero_ps();
    __m128 avgF_vec = _mm_setzero_ps();
    __m128 minF_vec = _mm_set_ps1(FLT_MAX); //TODO CHANGE THAT

    __m128 temp_num_0;
    __m128 temp_num_1;
    __m128 temp_num_2;
    __m128 temp_num;
    __m128 temp_den_0;
    __m128 temp_den_1;
    __m128 temp_den;
    __m128 temp_fvec;

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
        avgF_vec = _mm_setzero_ps();
        for(unsigned int i=0; i< 2*N/4 ;i+=2)
        {
            //Replace each statement step by step using the instructions as mentioned before

            /*
             * Unrolling and jamming not shown step by step by the idea is the following:
             * SSE works with vectors of 4 elements. So every operation needs four elements.
             * ex. when we add LVec_ptr[i] and RVec_ptr[i] sse adds the elements
             * LVec_ptr[i+0]+RVec_ptr[i+0], LVec_ptr[i+1]+RVec_ptr[i+1], LVec_ptr[i+2]+RVec_ptr[i+2],
             * LVec_ptr[i+3]+RVec_ptr[i+3]. That's why we use division by 4(every i corresponds to 4
             * elements). We see first the referrence code in comment and sse implementation follows
             */

            //float num_0 = LVec[i] + RVec[i];
            temp_num_0 = _mm_add_ps(LRVec_ptr[i],LRVec_ptr[i+1]);

            //float num_1 = mVec[i]*(mVec[i]-1.0f)/2.0f;
            temp_num_1 = _mm_sub_ps(mnVec_ptr[i], temp_one);
            temp_num_1 = _mm_mul_ps(mnVec_ptr[i], temp_num_1);
            temp_num_1 = _mm_div_ps(temp_num_1, temp_two);


            //float num_2 = nVec[i]*(nVec[i]-1.0f)/2.0f;
            temp_num_2 = _mm_sub_ps(mnVec_ptr[i+1], temp_one);
            temp_num_2 = _mm_mul_ps(mnVec_ptr[i+1], temp_num_2);
            temp_num_2 = _mm_div_ps(temp_num_2, temp_two);


            //float num = num_0/(num_1+num_2);
            temp_num = _mm_add_ps(temp_num_1,temp_num_2);
            temp_num = _mm_div_ps(temp_num_0, temp_num);

            //float den_0 = CVec[i]-LVec[i]-RVec[i];
            temp_den_0 = _mm_sub_ps(FCVec_ptr[i+1], LRVec_ptr[i]);
            temp_den_0 = _mm_sub_ps(temp_den_0, LRVec_ptr[i+1]);


            //float den_1 = mVec[i]*nVec[i];
            temp_den_1 = _mm_mul_ps(mnVec_ptr[i], mnVec_ptr[i+1]);

            //float den = den_0/den_1;
            temp_den = _mm_div_ps(temp_den_0, temp_den_1);

            //FVec[i] = num/(den+0.01f);
            temp_fvec = _mm_add_ps(temp_den, __temp_one);
            FCVec_ptr[i] = _mm_div_ps(temp_num, temp_fvec);

            //maxF = FVec[i]>maxF?FVec[i]:maxF;
            maxF_vec = _mm_max_ps(FCVec_ptr[i], maxF_vec);

            //minF = FVec[i]<minF?FVec[i]:minF;
            minF_vec = _mm_min_ps(FCVec_ptr[i], minF_vec);

            //avgF += FVec[i];
            avgF_vec = _mm_add_ps(FCVec_ptr[i], avgF_vec );
            //printf("HU %e, %e, %e, %e\n",(double)avgF_vec[0], (double)avgF_vec[1], (double)avgF_vec[2], (double)avgF_vec[3]);
        }

        /*
         * Prefer not to use for loop but unrolled version, so access all the values hardcoded
         */
        maxF = maxF_vec[0];
        maxF = maxCalc(maxF_vec[1],maxF);
        maxF = maxCalc(maxF_vec[2],maxF);
        maxF = maxCalc(maxF_vec[3],maxF);

        minF = minF_vec[0];
        minF = minCalc(minF_vec[1],minF);
        minF = minCalc(minF_vec[2],minF);
        minF = minCalc(minF_vec[3],minF);

        avgF =  avgF_vec[0] + avgF_vec[1] + avgF_vec[2] + avgF_vec[3];

        /*
         * Use a simple scalar computation to find whatever is left to compute.
         * NOTE: In our we use N = 10000000 so it is divisable by 4 and nothing
         * is left to compute. This code is added for generalization purposes
         * so that we have a correct result for a different N
         */
        //TODO ADD THAT
        for (int j = (2 * N - 2 * N % 8); j < 2* N; j+=2) {
            float num_0 = LRVec[j] + LRVec[j+1];
            float num_1 = mnVec[j] * (mnVec[j] - 1.0f) / 2.0f;
            float num_2 = mnVec[j+1] * (mnVec[j+1] - 1.0f) / 2.0f;
            float num = num_0 / (num_1 + num_2);
            float den_0 = FCVec[j+1] - LRVec[i] - LRVec[j+1];
            float den_1 = mnVec[j] * mnVec[j+1];
            float den = den_0 / den_1;

            FCVec[j] = num / (den + 0.01f);
            maxF = FCVec[j] > maxF ? FCVec[j] : maxF;
            minF = FCVec[j]<minF?FCVec[j]:minF;
            avgF += FCVec[j];
        }


    }
    double timeOmegaTotal = gettime()-timeOmegaTotalStart;
    double timeTotalMainStop = gettime();

    printf("Omega time %fs - Total time %fs - Min %e - Max %e - Avg %e\n",
           timeOmegaTotal/iters, timeTotalMainStop-timeTotalMainStart, (double)minF, (double)maxF,
           (double)avgF/N);

    /*
     * Free everything
     */
    _mm_free(mnVec);
    _mm_free(LRVec);
    _mm_free(FCVec);
}