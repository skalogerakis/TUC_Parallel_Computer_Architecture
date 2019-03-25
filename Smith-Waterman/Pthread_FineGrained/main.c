#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <memory.h>
#include <stdbool.h>
#include <stdint.h>
#include <pthread.h>


u_int64_t _cellVal = 0;
u_int64_t _traceSteps = 0;
double _totalTime = 0;
double timeTotal=0;
double _totalCellTime = 0;
double _totalTraceTime = 0;
double _CUPSTotal = 0;
double _CUPSCell = 0;

static char const UPCODE = 'U';
static char const LEFTCODE = 'L';
static char const DIAGCODE = 'D';
static char const ZEROCODE = 'Z';


int _MAX_SIMILARITY;
int SerialFlag = 0;
int balanceCount =0;
int barrierCounter = 0;

int Q_len = 0, D_len = 0 ;

int MATCH;
int MISMATCH;
int GAP;
int THREADS;

bool nameBool, inputBool, matchBool, misBool, gapBool,threadBool = false;

int pairs= -1;
int qMin = -1;
int qMax = -1;
int dSize = -1;

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutexattr_t mutexattr;
//pthread_mutex_t *swmutex;
pthread_barrier_t barrier;

//long long int antiDiagLen = 0;
long long int antiDiagNum = 1 ;

struct thread_data{
    long long int xCalc;
    long long int yCalc;
    int thread_id;
    long long int antiLen;
};

typedef struct thread_data _td;


FILE *finFile;

int **ScoreTable = NULL;
char * Q = NULL;
char * D = NULL;


/*
 * Add more input variables if you wish and update commandChecker function. ADDED EXTRA FUNCTIONALITY FOR THREADS
 */
const char* inVariable[] = {"-name", "-input", "-match", "-mismatch", "-gap" , "-threads"};

/*
 * This function is responsible to print all the demanded elements at the end of execution
 */
void terminalPrinter(){
    printf("\n\n");
    printf("NUMBER OF Q-D PAIRS: %d\n", pairs);
    printf("TOTAL NUMBER OF CELLS UPDATED: %lu\n",_cellVal);
    printf("TRACEBACK STEPS: %lu\n",_traceSteps);
    printf("TOTAL TIME : %lf\n", _totalTime);
    printf("TOTAL CELL TIME : %lf\n", _totalCellTime);
    printf("TOTAL TRACEBACK TIME : %lf\n", _totalTraceTime);
    printf("CUPS TOTAL TIME: %.3lf\n", _cellVal/_totalTime);
    printf("CUPS CELL TIME: %.3lf\n", _cellVal/_totalCellTime);

}

/*
 * This function calculates the current time.Used to calculate elapsed time
 */
double gettime(void)
{
    struct timeval ttime;
    gettimeofday(&ttime, NULL);
    return ttime.tv_sec+ttime.tv_usec * 0.000001;
}

/*
 * This function is responsible to check parameters from command line.
 * In case, the demanded criteria are not met then the program terminates.
 * NOTE: This function takes into consideration only the flags and takes as
 * their value next argv value. It does not check that this value is correct
 */
void commandChecker(int argc, char * argv[]){
    //printf("Argument count=%d\n",argc);

    /*
     * If you wish to add new case for command line do it here.
     */
    for (uint8_t i = 0; i < argc; i++) {
        if(strcmp(argv[i],inVariable[0]) == 0){
            nameBool = true;
        }
        if(strcmp(argv[i],inVariable[1]) == 0){
            inputBool = true;
        }
        if(strcmp(argv[i],inVariable[2]) == 0) {
            MATCH = atoi(argv[++i]);
            matchBool = true;
        }
        if(strcmp(argv[i],inVariable[3]) == 0) {
            MISMATCH = atoi(argv[++i]);
            misBool = true;
        }
        if(strcmp(argv[i],inVariable[4]) == 0) {
            GAP = atoi(argv[++i]);
            gapBool = true;
        }
        //THREAD CHECKER
        if(strcmp(argv[i],inVariable[5]) == 0) {
            THREADS = atoi(argv[++i]);
            threadBool = true;
        }
    }

    if( !nameBool | !inputBool | !matchBool | !misBool | !gapBool | !threadBool){
        printf("NOT ALL DEMANDED INPUT VARIABLES WERE GIVEN IN COMMAND LINE. EXITING.......");
        exit(-10);
    }else{
        printf("COMMAND LINE VARIABLES AS DEMANDED. CONTINUE\n");
    }
}

/*
 * This functions is used to check that all the criteria demanded from requirements
 * are met. In case we find an error terminate the program.
 */
void ErrorCode(int checker){
    if(checker < 0){
        printf("Error code . The testing file does not have the requested format. Terminating\n");
        exit(-9);
    }
}

/*
 * This function is responsible to reverse a string. Used at the end of backtracking algorithm
 * as the final string to be printed needs to be reversed
 */
char* reverseArr(char *str, size_t len) {
    size_t i = 0;
    while (len > i) {
        char tmp = str[--len];
        str[len] = str[i];
        str[i++] = tmp;
    }
    return str;
}

/*
 * This function is responsible to update all header file value.
 * In case the format is not correct then we exit the program
 */
void fileHeaderValues(FILE *fp){
    fscanf(fp, "Pairs: %d\n", &pairs);
    ErrorCode(pairs);
    //printf("\t\tNUMBER OF Q-D PAIRS: %d\n", pairs);

    fscanf(fp, "Q_Sz_Min: %d\n", &qMin);
    ErrorCode(qMin);
    //printf("Show pairs %d\n", qMin);

    fscanf(fp, "Q_Sz_Max: %d\n", &qMax);
    ErrorCode(qMax);
    //printf("Show pairs %d\n", qMax);


    fscanf(fp, "D_Sz_All: %d\n", &dSize);
    ErrorCode(dSize);
    //printf("Show pairs %d\n", dSize);
}


/*
 * Function that calculates the length of the anti-diagonal.
 * We have three seperate cases.
 * -Case 1:Elements are increasing(until we reach max length)
 * -Case 2:Elements are decreasing when we reach the end of the table
 * -Case 3:Elements are stable(Happens only once if our table is symmetrical which is not our case so we have some more)
 */
long long int antiDiagLength(long long int currAnti){
    long int minFinder = D_len < Q_len ? D_len : Q_len; //Our bottleneck
    long int maxFinder = D_len < Q_len ? Q_len : D_len;

    if(currAnti < D_len && currAnti < Q_len){
        //CASE 1
        return currAnti;
    }else if(currAnti > D_len && currAnti > Q_len ){
        //CASE 2
        return minFinder  - (currAnti - maxFinder  ) ;
    }else{
        //CASE 3
        return  minFinder;
    }
}

long long int * firstAntiElement(long long int currAnti, _td* tData){

//    _td *t;
//    tData = (_td *)threadArg;
    //static long long int coor[2];

    //long long tempX, tempY;
    if(currAnti <= Q_len){
        tData->xCalc = currAnti;
        tData->yCalc = 1;
//        tempY = 1;
//        tempX = currAnti;
        //coor[1] = currAnti;
        //coor[2] = 1;
    }else{
        tData->xCalc = Q_len;
        tData->yCalc = currAnti - Q_len + 1;
        //coor[1] = Q_len;
        //coor[2] = currAnti - Q_len + 1;
//        tempY = currAnti - Q_len + 1;
//        tempX = Q_len;
    }
    //printf("X %lld AND Y %lld\n",tempX,tempY);
    //return coor;
}

void updateScore(long long int x, long long int y){

    //printf("x %lld, y %lld\t", x, y);

    long long int up, left, diag;

    up = ScoreTable[x-1][y] + GAP;

    left = ScoreTable[x][y-1] + GAP; //we want same row and one column left

    long long int tempDiag,tempMax;


    //TODO wrong while printing.Fix that


    tempDiag = (D[y - 1] == Q[x - 1]) ? MATCH : MISMATCH;


    diag = ScoreTable[x - 1][y - 1] + tempDiag;

    //printf("D %c, Q %c up %lld, left %lld, diag %lld\n", D[x], Q[y], up, left, diag );

    if( up <= 0 && left <= 0 && diag <= 0 ){    //case we have max as negative value
        tempMax = 0;
    }else{
        if( up > left && up > diag){
            tempMax = up;
        }else if (left >= up && left > diag){
            tempMax = left;
        }else{
            tempMax = diag;
        }
        //_cellVal++;
    }

    ScoreTable[x][y] = tempMax;
    //printf("MAX %lld\n",tempMax);
}



//void *calcSimilarity(void *threadArg){
//    _td *tData;
//     tData = (_td *)threadArg;
//     //long long int *coor;
//
//    //printf("Hello Houston from thread %d\n\n\n",tData->thread_id);
//
//    //TODO FIX VARIABLE DECLARATION
//    //Don't start counting from the beginning. Eliminate elements from first line and column
//    printf("D size %d, Q size %d anti %lld\n", D_len, Q_len, antiDiagNum);
//
//    for (long long int i = 1; i<=antiDiagNum; i++){
//        //printf("Antidiagonal length %lld\n\n", antiDiagLength(i));
//        tData[tData->thread_id].antiLen = antiDiagLength(i);
//        //firstAntiElement(i, &tData[tData->thread_id]);
//
//        if(i <= Q_len){
//            tData[tData->thread_id].xCalc = i;
//            tData[tData->thread_id].yCalc = 1;
////        tempY = 1;
////        tempX = currAnti;
//            //coor[1] = currAnti;
//            //coor[2] = 1;
//        }else{
//            tData[tData->thread_id].xCalc = Q_len;
//            tData[tData->thread_id].yCalc = i - Q_len + 1;
//            //coor[1] = Q_len;
//            //coor[2] = currAnti - Q_len + 1;
////        tempY = currAnti - Q_len + 1;
////        tempX = Q_len;
//        }
//        //tData->xCalc = coor[1];
//        //tData->yCalc = coor[2];
//        //printf("Hello Houston from thread %d and len %lld\n",tData->thread_id, tData->antiLen);
//        //printf("X %lld, Y %lld\n",tData[i].xCalc, tData[i].yCalc);
//
//
//        //SerialFlag = 0;
//        //balanceCount = 1;
//
//        long long int lenDiv = tData[tData->thread_id].antiLen / THREADS;
//        long long int modDiv = tData[tData->thread_id].antiLen % THREADS;
//
//
//        long long int finEleX = (tData[tData->thread_id].xCalc-1) * lenDiv;
//        long long int finEleY = (tData[tData->thread_id].yCalc+1) * lenDiv;
//        //printf("calc X %lld,calc Y %lld and thread %d\n",tData[i].xCalc,tData[i].yCalc,tData->thread_id);
//        printf("FIN X %lld,FIN Y %lld and thread %d\n",finEleX,finEleY,tData->thread_id);
//        printf("DIV %lld, MOD %lld\n",lenDiv,modDiv);
//
//        for(long long j=0;j<lenDiv;j++){
//
//            long long int currX, currY;
//
//
//            currX = tData[tData->thread_id].xCalc - tData->thread_id;
//            currY = tData[tData->thread_id].yCalc + tData->thread_id;
//            //printf("TSO KAI LO %lld , %lld\n",currX,currY);
//
//            tData[tData->thread_id].xCalc = (tData[tData->thread_id].xCalc-1) * lenDiv;
//            tData[tData->thread_id].yCalc = (tData[tData->thread_id].yCalc+1) * lenDiv;
//            updateScore(currX, currY);
//
//        }
//
//        if(lenDiv == 0) {   //TODO works fine
//            pthread_mutex_lock(&lock);
//            SerialFlag++;
//            if (SerialFlag == 1) {
//                for (int j = 0; j < tData->antiLen; j++) {
//
//                    long long int currX, currY;
//                    currX = tData[tData->thread_id].xCalc - j;
//                    currY = tData[tData->thread_id].yCalc + j;
//
//                    //printf("Hello case X %lld, Y %lld and thread %d\n", currX, currY, tData->thread_id);
//                    updateScore(currX, currY);
//
//                }
//            }
//
////            SerialFlag = 1;
//
//            pthread_mutex_unlock(&lock);
//        }else if( modDiv!=0){
//            pthread_mutex_lock(&lock);
//
//            if(balanceCount <= modDiv){
//                long long int currX, currY;
//                currX = finEleX - balanceCount;
//                currY = finEleY + balanceCount;
//
//                printf("Balance case X %lld, Y %lld and thread %d\n",currX,currY,tData->thread_id);
//
//                updateScore(currX, currY);
//            }
//
//
//            balanceCount++;
//            pthread_mutex_unlock(&lock);
//        }
//        printf("WAIT thread %d\n",tData->thread_id);
//        pthread_barrier_wait(&barrier);
//        SerialFlag = 0;
//        balanceCount = 1;
//        printf("\n\n");
//
//    }
//
//    printf("WATI BAR\n\n");
//    pthread_barrier_wait(&barrier);
//    return ((void *) 0);
//}

void *calcSimilarity(void *threadArg){
//    _td *tData;
//    tData = (_td *)threadArg;
    long tid;
    tid = (long)threadArg;
    long long int diaLen,initX,initY;

    long long int currX, currY;
    //long long int *coor;

    //printf("Hello Houston from thread %d\n\n\n",tData->thread_id);

    //TODO FIX VARIABLE DECLARATION
    //Don't start counting from the beginning. Eliminate elements from first line and column
    //printf("D size %d, Q size %d anti %lld\n", D_len, Q_len, antiDiagNum);

    for (long long int i = 1; i<=antiDiagNum; i++){
        SerialFlag = 0;
        balanceCount = 0;
        //printf("Antidiagonal length %lld\n\n", antiDiagLength(i));

        //tData[tData->thread_id].antiLen = antiDiagLength(i);
        diaLen = antiDiagLength(i);

        //printf("DIA %lld\n\n", diaLen);
        //firstAntiElement(i, &tData[tData->thread_id]);

        if(i <= Q_len){
            //tData[tData->thread_id].xCalc = i;
            //tData[tData->thread_id].yCalc = 1;
            initX = i;
            initY = 1;
//        tempY = 1;
//        tempX = currAnti;
            //coor[1] = currAnti;
            //coor[2] = 1;
        }else{
            initX = Q_len;
            initY = i - Q_len + 1;
            //tData[tData->thread_id].xCalc = Q_len;
            //tData[tData->thread_id].yCalc = i - Q_len + 1;
            //coor[1] = Q_len;
            //coor[2] = currAnti - Q_len + 1;
//        tempY = currAnti - Q_len + 1;
//        tempX = Q_len;
        }

        //printf("DIA %lld, x %lld, y %lld on thread %ld \n\n", diaLen, initX, initY, tid);
        //tData->xCalc = coor[1];
        //tData->yCalc = coor[2];
        //printf("Hello Houston from thread %d and len %lld\n",tData->thread_id, tData->antiLen);
        //printf("X %lld, Y %lld\n",tData[i].xCalc, tData[i].yCalc);


        //SerialFlag = 0;
        //balanceCount = 1;

//        long long int lenDiv = tData[tData->thread_id].antiLen / THREADS;
//        long long int modDiv = tData[tData->thread_id].antiLen % THREADS;

        long long int lenDiv = diaLen / THREADS;
        long long int modDiv = diaLen % THREADS;

//        long long int finEleX = (tData[tData->thread_id].xCalc-1) * lenDiv;
//        long long int finEleY = (tData[tData->thread_id].yCalc+1) * lenDiv;

        long long int finEleX = initX - ((THREADS *lenDiv)-1);
        long long int finEleY = initY + ((THREADS *lenDiv)-1);


        //printf("calc X %lld,calc Y %lld and thread %d\n",tData[i].xCalc,tData[i].yCalc,tData->thread_id);
        //printf("FIN X %lld,FIN Y %lld and thread %d\n",finEleX,finEleY,tData->thread_id);
        printf("DIV %lld, MOD %lld and %ld\n",lenDiv,modDiv,tid);
        //printf("FIN X %lld,FIN Y %lld init X %lld initY %lld and thread %ld\n",finEleX,finEleY,initX,initY,tid);

        for(long long j=1;j<=lenDiv;j++){


            currX = initX - tid;
            currY = initY + tid;
            //currX = initX - ((tid *j));
            //currY = initY + ((tid *j));

            //printf("TSO KAI LO %lld , %lld\n",currX,currY);

            initX = initX - (THREADS);
            initY = initY + (THREADS);
            //initY = (initY) * (lenDiv+1);
            printf("MAIN, %lld, %lld, %ld",currX,currY,tid);
            updateScore(currX, currY);


        }

        if(lenDiv == 0) {   //TODO works fine

            pthread_mutex_lock(&lock);
            SerialFlag+=1;
            pthread_mutex_unlock(&lock);
            //printf("SERAL %d\n",SerialFlag);
            if (SerialFlag == 1) {
                //for (int j = 0; j < tData->antiLen; j++) {
                for (int j = 0; j < diaLen; j++) {
                    long long int currXa, currYa = 0;
                    currXa = initX - j;
                    currYa = initY + j;

                    //printf("Hello case X %lld, Y %lld and thread %ld\n", currXa, currYa, tid);
                    printf("DIV0, %lld, %lld, %ld",currXa,currYa,tid);
                    updateScore(currXa, currYa);

                }
            }

            //SerialFlag = 1;
        //}

        }else if( modDiv!=0 && tid < modDiv){
//            pthread_mutex_lock(&lock);
//            balanceCount+=1;
//            pthread_mutex_unlock(&lock);

            //if(balanceCount <= modDiv){
                long long int currXb, currYb =0;
                currXb = finEleX - (tid +1);
                currYb = finEleY + (tid+1);

                //printf("Balance case X %lld, Y %lld and thread %ld bal %d\n",currXb,currYb,tid, balanceCount);
                printf("MOD, %lld, %lld, %ld",currXb,currYb,tid);
                updateScore(currXb, currYb);
            //}



        }
        //printf("WAIT thread %d\n",tData->thread_id);
        pthread_barrier_wait(&barrier);
        //printf("WAIT thread %ld\n",tid);
//        pthread_mutex_lock(&swmutex[tid]);
//        barrierCounter++;
//        pthread_mutex_lock(&swmutex[tid]);
//        pthread_mutex_lock(&lock);
//        SerialFlag = 0;
//        balanceCount = 0;
//        pthread_mutex_unlock(&lock);
        //printf("NEXT \n\n");
        printf("\n\n");
    }

    pthread_barrier_wait(&barrier);

    //printf("WATI BAR\n\n");
    //pthread_barrier_wait(&barrier);
    pthread_exit(NULL);
}

/*
 * This function is where the algorithm and the calculations take place.
 * It takes as parameters the 2 strings that we need to apply the algorithm
 * and the game begins
 */

void dataParser(){

    //At first write as output the strings as they are
    //TODO ENABLE
    //fprintf(finFile, "\n\nQ: \t%s", Q);
    //fprintf(finFile, "\nD: \t%s\n", D);

    Q_len = strlen(Q);
    D_len = strlen(D);

    /*
     * We use and initialize.The first one is responsible to mark the scores from the
     * calculations.The second contains backtracking information.
     */
//    int *ScoreTable[Q_len+1];
//    for(int i = 0; i< Q_len+1; i++){
//        ScoreTable[i] = malloc((D_len+1)* sizeof(*ScoreTable[i]));
//        if(ScoreTable[i] == NULL){
//            printf("Could not allocate memory for matrix ScoreTable.Terminating");
//            exit(-2);
//        }
//    }
    ScoreTable = malloc((Q_len+1)* sizeof(int *));
    if(ScoreTable ==NULL){
        printf("Could not allocate memory for matrix ScoreTable.Terminating");
        exit(-2);
    }
    for(int i = 0; i< Q_len+1; i++){
        ScoreTable[i] = malloc((D_len+1)* sizeof(int));
        if(ScoreTable[i] == NULL){
            printf("Could not allocate memory for matrix ScoreTable.Terminating");
            exit(-2);
        }
    }

//    ScoreTable = malloc((Q_len+1) * sizeof(int *));
//    for(int i =0 ; Q_len + 1; i++){
//        ScoreTable[i] = malloc((D_len+1)* sizeof(int));
//    }

//    for(int i = 0; i< Q_len+1; i++){
//        ScoreTable[i] = malloc((D_len+1)* sizeof(*ScoreTable[i]));
//        if(ScoreTable[i] == NULL){
//            printf("Could not allocate memory for matrix ScoreTable.Terminating");
//            exit(-2);
//        }
//    }
    printf("Memory Allocated successfully for ScoreTable matrix.\n");

    //Initialize 2D array (We have one additional column and row)

    for(int i=0; i < Q_len + 1; i++){
        ScoreTable[i][0] = 0;
    }

    for(int i=0; i < D_len + 1; i++){
        ScoreTable[0][i] = 0;
    }






    //TODO HERE WE GO

    /*
     * Mutex Init should return zero. Otherwise print error message and exit
     */

    pthread_mutexattr_init(&mutexattr);
    pthread_mutexattr_settype(&mutexattr,PTHREAD_MUTEX_NORMAL);

    if(pthread_mutex_init(&lock,&mutexattr) != 0){
        printf("Error occured while trying to init mutexes. Exitiing...");
        exit(-1);
    }
//    swmutex = malloc(sizeof(pthread_mutex_t)*(THREADS+1));
//    for(int i = 0;i<THREADS+1; i++){
//        if(pthread_mutex_init(&swmutex[i],NULL) != 0){
//            printf("Error occured while trying to init mutexes. Exitiing...");
//            exit(-1);
//        }
//    }


    pthread_t mthread[THREADS];

    pthread_barrier_init(&barrier, NULL,THREADS);

    //_td *tData = malloc(THREADS * sizeof(*tData));
    //memset(tData, 0, sizeof(_td));
    //_td tData[THREADS];

    antiDiagNum = D_len + Q_len  - 1;

    /*
     * Create the demanded number of threads on the beginning
     */
    //_td *tData;
    int threadChecker;
    for(long i = 0; i<THREADS; i++){
        //TODO LAST PARAMETER CAN PASS ARGUMENTS.CHECK IF WE CAN USE THAT
        /*
         * NULL for attributes as we want to use default
         */
        //printf("MAIN : Creating thread %ld\n",i);
        //tData = (_td *)malloc(THREADS * sizeof(tData));
//        tData[i].thread_id = i;
//        tData[i].yCalc = 0;
//        tData[i].xCalc = 0;
//        tData[i].antiLen = 1;
        //threadChecker = pthread_create(&mthread[i], NULL, calcSimilarity, (void *)&tData[i]);
        threadChecker = pthread_create(&mthread[i], NULL, calcSimilarity, (void *)i);
        if(threadChecker){
            printf("Thread creation failed.Exiting");
            exit(-1);
        }


    }

    /*
     * After the parallel part is over join all the threads and move on
     */

    for(int i = 0; i<THREADS; i++){
       // printf("MAIN : completed join with thread %d \n\n",tData[i].thread_id);
        pthread_join(mthread[i], NULL);

    }

    pthread_barrier_destroy(&barrier);
    pthread_mutex_destroy(&lock);
//    for(int i = 0;i<THREADS+1; i++){
//        pthread_mutex_destroy(&swmutex[i]);
//    }

    //free(mthread);
//    for(int i = 0; i< THREADS; i++) {
//        free(tData);
//    }

    //free(tData);




    printf("START SIMILARITY MATRIX\n");

    /*
     * Similarity matrix calculations. We also find how many
     * max values we have in each matrix
     */
//    int _counterMax = 0;
//
//    double cellTimeInit = gettime();
//    _MAX_SIMILARITY=0;
//
//    for(int i = 1; i < Q_len + 1; i++){
//        for(int j = 1; j < D_len + 1; j++){
//            int up, left, diag;
//
//            up = ScoreTable[i-1][j] + GAP;  //we want same column and one row up
//            left = ScoreTable[i][j-1] + GAP; //we want same row and one column left
//
//            //for diagonal check if we have match or mismatch
//            uint tempDiag;
//
//            tempDiag = (D[j - 1] == Q[i - 1]) ? MATCH : MISMATCH;
//
//            diag = ScoreTable[i - 1][j - 1] + tempDiag;
//
//            uint tempMax = 0;
//
//            if( up <= 0 && left <= 0 && diag <= 0 ){    //case we have max as negative value
//                tempMax = 0;
//            }else{
//                if( up > left && up > diag){
//                    tempMax = up;
//                }else if (left >= up && left > diag){
//                    tempMax = left;
//                }else{
//                    tempMax = diag;
//                }
//                _cellVal++;
//            }
//
//            ScoreTable[i][j] = tempMax;
//
//            if(tempMax > _MAX_SIMILARITY) {
//                _counterMax = 0;
//                _MAX_SIMILARITY = tempMax;
//            }else if(tempMax == _MAX_SIMILARITY && tempMax!=0){
//                _counterMax++;
//            }
//        }
//    }
//
//
//    double cellTimeFin = gettime();
//    _totalCellTime+= (cellTimeFin-cellTimeInit);

    printf("Similarity matrix done. Start backtracking\n");









    //printf("MAX %d\n",_MAX_SIMILARITY);

    /*
     * FIND ALL MAX VALUES FOR BACKTRACKING
     */
    //TODO IN THE END ENABLE BACKTRACKING AGAIN

//    long xMax[_counterMax+1];
//    long yMax[_counterMax+1];
//    long _endKeeper[_counterMax+1];
//
//    int tempCount=0;
//    for (int i = 0; i < Q_len + 1; i++){
//        for (int j = 0; j < D_len + 1; j++){
//            if(_MAX_SIMILARITY == ScoreTable[i][j]){
//                _endKeeper[tempCount] = j;
//                xMax[tempCount] = i;
//                yMax[tempCount] = j;
//                tempCount++;
//            }
//        }
//    }
//
//
//    /*
//     * BACKTRACKING. We continue with the loop until it hits zero
//     */
//
//    double traceTimeInit = gettime();
//    for(int i = 0; i< _counterMax+1; i++){
//
//        long currXpos = xMax[i];
//        long currYpos = yMax[i];
//        char xElem = Q[currXpos];
//        char yElem = D[currYpos];
//        char currNode = ScoreTable[currXpos][currYpos];
//
//        char *_qOut = NULL;
//
//        _qOut = (char *)malloc((xMax[i]+1)*(yMax[i]+1)*sizeof(char));
//        if(_qOut==NULL){
//            printf("Error occured while trying to allocate memory for traceback.Terminating....");
//            exit(-1);
//        }
//
//        char *_dOut = NULL;
//        _dOut = (char *)malloc((xMax[i]+1)*(yMax[i]+1)* sizeof(char));
//        if(_dOut==NULL){
//            printf("Error occured while trying to allocate memory for traceback.Terminating....");
//            exit(-1);
//        }
//
//        int lengthCount = 0;
//
//        //We initialize traceflag to zero when we want to find a new max value
//        //Traceback until we reach zero. Added some cases to make it even more strict
//        u_int8_t traceFlag = 0;
//        while(traceFlag != 1){
//            int up, left, diag;
//
//            up = ScoreTable[currXpos-1][currYpos];  //we want same column and one row up
//            left = ScoreTable[currXpos][currYpos-1]; //we want same row and one column left
//            diag = ScoreTable[currXpos - 1][currYpos - 1];
//
//            if(diag == 0 && ScoreTable[currXpos][currYpos] - MATCH == 0){
//                traceFlag = 1;
//                currYpos--;
//                currXpos--;
//                xElem = Q[currXpos];
//                yElem = D[currYpos];
//            }else{
//                if( up > left && up > diag){
//                    currXpos--;
//                    xElem = Q[currXpos];
//                    yElem = '-';
//                }else if (left > up && left > diag){
//                    currYpos--;
//                    xElem = '-';
//                    yElem = D[currYpos];
//                }else if(diag >= up && diag >= left){
//                    if(Q[currXpos-1] == D[currYpos-1] || (diag > up && diag > left) ){
//                        currYpos--;
//                        currXpos--;
//                        xElem = Q[currXpos];
//                        yElem = D[currYpos];
//                    }else if(Q[currXpos-1] == D[currYpos]){
//                        currXpos--;
//                        xElem = Q[currXpos];
//                        yElem = '-';
//                    }else{
//                        currYpos--;
//                        xElem = '-';
//                        yElem = D[currYpos];
//                    }
//
//                }else{
//                    if(Q[currXpos] == D[currYpos-1]){
//                        currYpos--;
//                        xElem = '-';
//                        yElem = D[currYpos];
//                    }else{
//                        currXpos--;
//                        xElem = Q[currXpos];
//                        yElem = '-';
//                    }
//                }
//
//            }
//            _traceSteps++;
//            _qOut[lengthCount] = xElem;
//            _dOut[lengthCount] = yElem;
//            lengthCount++;
//
//        }
//
//        _qOut[lengthCount] = '\0';
//        _dOut[lengthCount] = '\0';
//
//        _dOut = reverseArr(_dOut, strlen(_dOut));
//        _qOut = reverseArr(_qOut, strlen(_qOut));
//
//        //printf("\nQ reversed %s and \nD reversed %s",_qOut,_dOut);
//        /*
//         * Write all the info demanded on an output file
//         */
//        //TODO ENABLE IT AT THE END
//        fprintf(finFile, "\nMATCH %d [SCORE: %d,START: %ld,STOP: %ld]\n\tD: %s\n\tQ: %s\n", i+1, _MAX_SIMILARITY, (_endKeeper[i]-lengthCount) , _endKeeper[i]-1, _dOut, _qOut);
//
//        free(_qOut);
//        free(_dOut);

        //printf("\n\n");
//    }
//
//    double traceTimeFin = gettime();
//    _totalTraceTime+=(traceTimeFin - traceTimeInit);

    //TODO DISABLE IT ON PRINT


    //PRINTER OF THE  ARRAYS. USED FOR DEBBUGING PURPOSES ONLY
    for (int i = 0; i < Q_len + 1; i++){
        printf("\n");
        for (int j = 0; j < D_len + 1; j++){
            printf("%d\t",ScoreTable[i][j]);
        }

    }

    for(int i = 0; i< Q_len+1; i++) {
        free(ScoreTable[i]);
    }


}
/*
 * This function is related with file parser. We find the point where the buffer
 * hits the 2nd Q that find and we return that position which will continue in the next loop
 */
long fillDataBuffer(char * buf, long bytereader ,int compFlag){

    long index = 0;
    long qIndex = 0;
    long dIndex = 0;
    //TODO check that memory was allocated successfully


    Q = (char *)malloc(bytereader* sizeof(char));
    if(Q==NULL){
        printf("Error occured while trying to allocate memory for buffer.Terminating....");
        exit(-1);
    }

    D = (char *)malloc(bytereader*sizeof(char));
    if(D==NULL){
        printf("Error occured while trying to allocate memory for buffer.Terminating....");
        exit(-1);
    }

    uint8_t dFlag=0;
    uint8_t qFlag=0;
    uint8_t qCount=0;

    size_t bufLen = strlen(buf);

    if(compFlag == 1) return -1;

    for(index;index < bufLen; index++){
        if(buf[index] == ':' || buf[index] == '\\' || buf[index] == '\t'|| buf[index] == '\n' || buf[index] == '\r'){
            continue;
        }

        if ( buf[index] == 'Q' ){
            qFlag = 1;
            qCount++;
            if(qCount>1){

                Q[qIndex] = '\0';
                D[dIndex] = '\0';

                //dataParser(Q, D);
                dataParser();

                free(Q);
                free(D);
                return index;
            }
            continue;
        }else if(buf[index] == 'D'){
            qFlag = 0;
            dFlag = 1;
            continue;
        }

        if(qFlag == 1){
            Q[qIndex] = buf[index];
            qIndex++;
        }else if(dFlag == 1){
            D[dIndex] = buf[index];
            dIndex++;
        }
    }

    Q[qIndex] = '\0';
    D[dIndex] = '\0';

    //dataParser(Q, D);
    dataParser();

    free(Q);
    free(D);
    return bufLen;
}

/*
 * This function is used to read input from file. As files are too large,
 * c functions such as sscanf, fgets cannot be applied. For that reason we use
 * buffer that we fill until we find our escape in which case we return. The next
 * buffer starts from that point until EOF. More info about that can be found on the
 * link : http://www.softwareprojects.com/resources/programming/t-processing-large-files-in-c-1636.html?fbclid=IwAR3H8v6N2xoVIFXxTzUswMNSsNPzYrm439KkLcUMCxmtb2mrMWRPEvOeF7w
 */
void fileParser(FILE *fp){
    //TODO WHAT IS GOING ON WITH BUFSIZE AND D8
    int BUFSIZE =  dSize+ qMax + 10000;
    long numbytes = 0;
    char buf[BUFSIZE];
    int  sizeLeftover=0;
    int isCompleted = 0;
    long pos = 0;
    int pairCounter = 0;

    do
    {

        // Read next block from file and save into buf, right after the
        // "left over" buffer
        numbytes = fread(buf+sizeLeftover, 1, sizeof(buf)-1-sizeLeftover, fp);
        if (numbytes<1)
        {
            // Turn on 'loop completed' flag so that we know to exit at the bottom
            // Still need to process any block we currently have in the
            // leftover buffer
            isCompleted = 1;
            numbytes  = 0;
        }

        // Add NULL terminator at the end of our buffer
        buf[numbytes+sizeLeftover] = 0;

        // Process data - Replace with your function
        //
        // Function should return the position in the file or -1 if failed
        //
        // We are also passing bLoopCompleted to let ProcessData know whether this is
        // the last record (in which case - if no end-of-record separator,
        // use eof and process anyway)

        pos = fillDataBuffer(buf, numbytes+sizeLeftover,isCompleted);
        pairCounter++;
        // If error occured, bail
        if (pos<1)
        {
            isCompleted = 1;
            pos      = 0;
        }

        sizeLeftover = numbytes+sizeLeftover-pos;
        // Extra protection - should never happen but you can never be too safe
        if (sizeLeftover<1) sizeLeftover=0;

        // If we have a leftover unprocessed buffer, move it to the beginning of
        // read buffer so that when reading the next block, it will connect to the
        // current leftover and together complete a full readable line
        if (pos!=0 && sizeLeftover!=0)
            memmove(buf, buf+pos, sizeLeftover);

    } while(!isCompleted && pairCounter != pairs);
}


//TODO finalize input from terminal and directories
int main(int argc, char * argv[]) {

    double timeInitTotal = gettime();

    commandChecker(argc, argv);

    printf("STARTING EXECUTION....\n");

    FILE *fp;

    fp = fopen("/home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Smith-Waterman/MyDocs/D1.txt","r");

    if(fp == NULL){
        printf("Error opening file\n");
        exit(-9);
    }

    //TODO ENABLE
    //finFile = fopen("/home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Smith-Waterman/MyDocs/FINAL.txt","a");
    printf("THREADS : %d\n",THREADS);

//    if(finFile == NULL){
//        printf("Error while opening write file!\n");
//        exit(1);
//    }

    fileHeaderValues(fp);
    fileParser(fp);


    fclose(fp);
    //TODO ENABLE
    //fclose(finFile);

    double timeFinTotal = gettime();
    _totalTime = timeFinTotal - timeInitTotal;

    //terminalPrinter();

    return 0;

}


