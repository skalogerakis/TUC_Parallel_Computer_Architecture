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


#include <pthread.h>
#include <unistd.h>

#define EXIT 127
#define BUSYWAIT 0
#define COMPUTE_ANTIDIAGONAL 1
#define COMPUTE_ANTIDIAGONAL_P2_A 2
#define COMPUTE_ANTIDIAGONAL_P2_B 3
#define COMPUTE_ANTIDIAGONAL_P3 4


typedef struct
{
    float * mVec;
    float * nVec;
    float * LVec;
    float * RVec;
    float * CVec;
    float * FVec;
    unsigned int N;
    float avgF;
    float maxF;
    float minF;
    //TODO MAYBE ADD __m128

} arrayCalc;


typedef struct
{
    int threadID;
    int threadTOTAL;
    int threadBARRIER;
    int threadOPERATION;

    arrayCalc arrayCalc;

} threadData_t;

static pthread_t * workerThread;


static pthread_barrier_t barrier;

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
/*
 * Finds the max value between two values. No need for ready function
 */
float maxCalc(float x, float y){
    if(x > y){
        return x;
    }
    return y;
}

/*
 * Finds the min value between two values. No need for ready function
 */
float minCalc(float x, float y){
    if(x < y){
        return x;
    }
    return y;
}



void setThreadArguments_COMPUTE_ANTIDIAGONAL_operation (threadData_t * threadData, int tid, int diagSize, int diagDoffset, int diagIndex, int jstart0)
{
//    threadData[tid].antidiagonalData.diagSize=diagSize;
//    threadData[tid].antidiagonalData.diagDoffset=diagDoffset;
//    threadData[tid].antidiagonalData.diagIndex=diagIndex;
//    threadData[tid].antidiagonalData.jstart0=jstart0;
}

void updateThreadArguments_COMPUTE_ANTIDIAGONAL_operation (threadData_t * threadData, int diagSize, int diagDoffset, int diagIndex, int jstart0)
{
    int threadIndex = 0;
    for(threadIndex=0;threadIndex<threadData->threadTOTAL;threadIndex++)
        setThreadArguments_COMPUTE_ANTIDIAGONAL_operation (threadData, threadIndex, diagSize, diagDoffset, diagIndex, jstart0);

}

void initializeThreadData(threadData_t * cur, int i, int threads, float *mVec, float *nVec, float *LVec, float *RVec, float *CVec, float  *FVec,
                          unsigned int N, float avgF, float maxF, float minF)
{
    printf("INIT %d\n",i);
    cur->threadID=i;
    cur->threadTOTAL=threads;
    cur->threadBARRIER=0;
    cur->threadOPERATION=BUSYWAIT;
    cur->arrayCalc.mVec = mVec;
    cur->arrayCalc.nVec = nVec;
    cur->arrayCalc.LVec = LVec;
    cur->arrayCalc.RVec = RVec;
    cur->arrayCalc.CVec = CVec;
    cur->arrayCalc.FVec = FVec;
    cur->arrayCalc.N = N;
    cur->arrayCalc.avgF = avgF;
    cur->arrayCalc.maxF = maxF;
    cur->arrayCalc.minF = minF;
}


void * thread (void * x)
{
    threadData_t * currentThread = (threadData_t *) x;

    int tid = currentThread->threadID;

    int threads = currentThread->threadTOTAL;
    printf("HI THREAD %d\n",tid);
    while(1)
    {
        __sync_synchronize();

        if(currentThread->threadOPERATION==EXIT)
            return NULL;

        //TODO CHANGE ALL THAT
        if(currentThread->threadOPERATION==COMPUTE_ANTIDIAGONAL)
        {
            //computeAntidiagonal (currentThread);

            currentThread->threadOPERATION=BUSYWAIT;

            pthread_barrier_wait(&barrier);

        }
    }

    return NULL;
}

static inline void setThreadOperation(threadData_t * threadData, int operation)
{
    int i, threads=threadData[0].threadTOTAL;

    for(i=0;i<threads;i++)
        threadData[i].threadOPERATION = operation;
}

static inline void execFunctionMaster(threadData_t * threadData, int operation)
{
//    if(operation==COMPUTE_ANTIDIAGONAL)
//        computeAntidiagonal(&threadData[0]);
//
//    if(operation==COMPUTE_ANTIDIAGONAL_P2_A)
//        computeAntidiagonal_P2_A(&threadData[0]);
//
//    if(operation==COMPUTE_ANTIDIAGONAL_P2_B)
//        computeAntidiagonal_P2_B(&threadData[0]);
//
//    if(operation==COMPUTE_ANTIDIAGONAL_P3)
//        computeAntidiagonal_P3(&threadData[0]);

}

static inline void syncThreadsBARRIER(threadData_t * threadData)
{
    int i, threads = threadData[0].threadTOTAL, barrierS=0;

    threadData[0].threadOPERATION=BUSYWAIT;

    pthread_barrier_wait(&barrier);

}

void startThreadOperations(threadData_t * threadData, int operation)
{
    setThreadOperation(threadData, operation);

    execFunctionMaster(threadData,operation);

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


int main(int argc, char ** argv)
{
    assert(argc==2);

    double timeTotalMainStart = gettime();
    float avgF = 0.0f;
    float maxF = 0.0f;
    float minF = FLT_MAX;
    unsigned int N = (unsigned int)atoi(argv[1]);
    unsigned int iters = 10;


    //THREAD PART STARTS HERE
    int threads = 4;




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

    __m128 temp_num_0 = _mm_setzero_ps();
    __m128 temp_num_1 = _mm_setzero_ps();
    __m128 temp_num_2 = _mm_setzero_ps();
    __m128 temp_num = _mm_setzero_ps();
    __m128 temp_den_0 = _mm_setzero_ps();
    __m128 temp_den_1 = _mm_setzero_ps();
    __m128 temp_den = _mm_setzero_ps();


    int s = pthread_barrier_init(&barrier, NULL, (unsigned int)threads);
    assert(s == 0);


    workerThread = NULL;
    workerThread = (pthread_t *) malloc (sizeof(pthread_t)*((unsigned long)(threads-1)));

    threadData_t * threadData = (threadData_t *) malloc (sizeof(threadData_t)*((unsigned long)threads));
    assert(threadData!=NULL);

    for(int i=0;i<threads;i++)
        initializeThreadData(&threadData[i],i,threads, mVec, nVec, LVec, RVec, CVec, FVec, N, avgF, maxF, minF );

    for(int i=1;i<threads;i++)
        pthread_create (&workerThread[i-1], NULL, thread, (void *) (&threadData[i]));

    /*
     * Initialize all variables
     */
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
        //Re-initialize in each iteration.
        avgF = 0.0f;
        maxF = 0.0f;
        minF = FLT_MAX;
        avgF_vec = _mm_setzero_ps();
        maxF_vec = _mm_setzero_ps();
        minF_vec = _mm_set_ps1(FLT_MAX);

        //Initialize for all threads.
        for(int i=0;i<threads;i++)
            initializeThreadData(&threadData[i],i,threads, mVec, nVec, LVec, RVec, CVec, FVec, N, avgF, maxF, minF );

        for(unsigned int i=0; i<N/4 ;i++)
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
            avgF_vec = _mm_add_ps(FVec_ptr[i], avgF_vec );
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
        for (int j = (N - N % 4); j < N; j++) {
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

    //TODO ADDED
    terminateWorkerThreads(workerThread,threadData);
    pthread_barrier_destroy(&barrier);
    if(threadData!=NULL)
        free(threadData);

    threadData = NULL;

    _mm_free(mVec);
    _mm_free(nVec);
    _mm_free(LVec);
    _mm_free(RVec);
    _mm_free(CVec);
    _mm_free(FVec);
}