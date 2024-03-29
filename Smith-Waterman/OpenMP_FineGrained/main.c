#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <memory.h>
#include <stdbool.h>
#include <stdint.h>
#include <unistd.h>
#include <omp.h>


u_int64_t _cellVal = 0;
u_int64_t _traceSteps = 0;
double _totalTime = 0;

double _totalCellTime = 0;
double _totalTraceTime = 0;

int _MAX_SIMILARITY;
int _counterMax;

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
char *INPUT;
char *REPORT;
char threadC[2];
char cwd[512];


int antiDiagNum = 1 ;

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
 * This function is responsible to create final path for output file.
 * Gets current directory and concats it with value from command line
 */
void finFileDir(){
    //printf("INPUT %s\n",INPUT);

    if(getcwd(cwd, sizeof(cwd))==NULL){

        printf("Error while trying to get current directory");
        exit(-1);
    }
    /*
     * Here concat any strings you want to create your directory
     */
    sprintf(threadC,"%d",THREADS);
    strcat(cwd,"/Report_");
    strcat(cwd,REPORT);
    strcat(cwd,"_OMP_");
    strcat(cwd,threadC);
    strcat(cwd,".txt");

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
            REPORT = argv[++i];
            nameBool = true;
        }
        if(strcmp(argv[i],inVariable[1]) == 0){
            INPUT = argv[++i];

            /*
             * We check that our input path end in .txt
             */
            if(strcmp(&INPUT[strlen(INPUT)-4],".txt")!=0){
                printf("Your input path does not end in .txt . Exiting.\n");
                exit(-9);
            }
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

    fscanf(fp, "Q_Sz_Min: %d\n", &qMin);
    ErrorCode(qMin);

    fscanf(fp, "Q_Sz_Max: %d\n", &qMax);
    ErrorCode(qMax);


    fscanf(fp, "D_Sz_All: %d\n", &dSize);
    ErrorCode(dSize);
}


/*
 * Function that calculates the length of the anti-diagonal.
 * We have three seperate cases.
 * -Case 1:Elements are increasing(until we reach max length)
 * -Case 2:Elements are decreasing when we reach the end of the table
 * -Case 3:Elements are stable(Happens only once if our table is symmetrical which is not our case so we have some more)
 */
int antiDiagLength(int currAnti){
    int minFinder = D_len < Q_len ? D_len : Q_len; //Our bottleneck
    int maxFinder = D_len < Q_len ? Q_len : D_len;

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

/*
 * Function especially designed for threds and open mp parallel program that does all the calculations needed to apply
 * Smith-Waterman algorithm
 */
int updateScore(long long int x, long long int y){


    int up, left, diag;

    up = ScoreTable[x-1][y] + GAP;

    left = ScoreTable[x][y-1] + GAP; //we want same row and one column left

    int tempDiag,tempMax;
    int tempCell=0;


    tempDiag = (D[y - 1] == Q[x - 1]) ? MATCH : MISMATCH;


    diag = ScoreTable[x - 1][y - 1] + tempDiag;


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

        tempCell++;
    }

    ScoreTable[x][y] = tempMax;


    if(tempMax > _MAX_SIMILARITY) {
        _counterMax = 0;
        _MAX_SIMILARITY = tempMax;
    }else if(tempMax == _MAX_SIMILARITY && tempMax!=0){
        _counterMax++;
    }

    return tempCell;
}


/*
 * This function is where the algorithm and the calculations take place.
 * It takes as parameters the 2 strings that we need to apply the algorithm
 * and the game begins
 */

void dataParser(){

    //At first write as output the strings as they are
    //TODO ENABLE
    fprintf(finFile, "\n\nQ: \t%s", Q);
    fprintf(finFile, "\nD: \t%s\n", D);

    Q_len = strlen(Q);
    D_len = strlen(D);

    /*
     * We use and initialize.The first one is responsible to mark the scores from the
     * calculations.The second contains backtracking information.
     */

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


    //printf("Memory Allocated successfully for ScoreTable matrix.\n");

    //Initialize 2D array (We have one additional column and row)

    for(int i=0; i < Q_len + 1; i++){
        ScoreTable[i][0] = 0;
    }

    for(int i=0; i < D_len + 1; i++){
        ScoreTable[0][i] = 0;
    }

    _counterMax = 0;
    _MAX_SIMILARITY=0;


    antiDiagNum = D_len + Q_len  - 1;

    double cellTimeInit = gettime();


    //printf("START SIMILARITY MATRIX\n");

    /*
     * This is where omp parallelization starts.
     */


        int diaLen,initX, initY,sum;



        for ( int i = 1; i<=antiDiagNum; i++){
            sum = 0;
            diaLen = antiDiagLength(i);
            if(i <= Q_len){
                initX = i;
                initY = 1;
            }else{
                initX = Q_len;
                initY = i - Q_len + 1;

            }


            //TODO SEE WHAT WE CAN DO WITH THAT.BOTTLENECK
//            if(diaLen < Q_len){
//                for(int k =0;k< diaLen;k++){
//                    _cellVal +=updateScore(initX - k,initY + k);
//                }
//                continue;
//            }

            int j=0;
#pragma omp parallel for num_threads(THREADS) default(none)\
            shared(initX,initY,diaLen,THREADS) private(j) reduction(+:sum)
                for(j =0;j< diaLen;j++){

                        sum +=updateScore(initX - j,initY + j);
                }
                _cellVal+=sum;
            }


    double cellTimeFin = gettime();
    _totalCellTime+= (cellTimeFin-cellTimeInit);



    //printf("Similarity matrix done. Start backtracking\n");



    //printf("MAX %d\n",_MAX_SIMILARITY);

    /*
     * FIND ALL MAX VALUES FOR BACKTRACKING
     */
    //TODO IN THE END ENABLE BACKTRACKING AGAIN

    long xMax[_counterMax+1];
    long yMax[_counterMax+1];
    long _endKeeper[_counterMax+1];

    int tempCount=0;
    for (int i = 0; i < Q_len + 1; i++){
        for (int j = 0; j < D_len + 1; j++){
            if(_MAX_SIMILARITY == ScoreTable[i][j]){
                _endKeeper[tempCount] = j;
                xMax[tempCount] = i;
                yMax[tempCount] = j;
                tempCount++;
            }
        }
    }


    /*
     * BACKTRACKING. We continue with the loop until it hits zero
     */

    double traceTimeInit = gettime();
    for(int i = 0; i< _counterMax+1; i++){

        long currXpos = xMax[i];
        long currYpos = yMax[i];
        char xElem = Q[currXpos];
        char yElem = D[currYpos];
        char currNode = ScoreTable[currXpos][currYpos];

        char *_qOut = NULL;

        _qOut = (char *)malloc((xMax[i]+1)*(yMax[i]+1)*sizeof(char));
        if(_qOut==NULL){
            printf("Error occured while trying to allocate memory for traceback.Terminating....");
            exit(-1);
        }

        char *_dOut = NULL;
        _dOut = (char *)malloc((xMax[i]+1)*(yMax[i]+1)* sizeof(char));
        if(_dOut==NULL){
            printf("Error occured while trying to allocate memory for traceback.Terminating....");
            exit(-1);
        }

        int lengthCount = 0;

        //We initialize traceflag to zero when we want to find a new max value
        //Traceback until we reach zero. Added some cases to make it even more strict
        u_int8_t traceFlag = 0;
        while(traceFlag != 1){
            int up, left, diag;

            up = ScoreTable[currXpos-1][currYpos];  //we want same column and one row up
            left = ScoreTable[currXpos][currYpos-1]; //we want same row and one column left
            diag = ScoreTable[currXpos - 1][currYpos - 1];

            if(diag == 0 && ScoreTable[currXpos][currYpos] - MATCH == 0){
                traceFlag = 1;
                currYpos--;
                currXpos--;
                xElem = Q[currXpos];
                yElem = D[currYpos];
            }else{
                if( up > left && up > diag){
                    currXpos--;
                    xElem = Q[currXpos];
                    yElem = '-';
                }else if (left > up && left > diag){
                    currYpos--;
                    xElem = '-';
                    yElem = D[currYpos];
                }else if(diag >= up && diag >= left){
                    if(Q[currXpos-1] == D[currYpos-1] || (diag > up && diag > left) ){
                        currYpos--;
                        currXpos--;
                        xElem = Q[currXpos];
                        yElem = D[currYpos];
                    }else if(Q[currXpos-1] == D[currYpos]){
                        currXpos--;
                        xElem = Q[currXpos];
                        yElem = '-';
                    }else{
                        currYpos--;
                        xElem = '-';
                        yElem = D[currYpos];
                    }

                }else{
                    if(Q[currXpos] == D[currYpos-1]){
                        currYpos--;
                        xElem = '-';
                        yElem = D[currYpos];
                    }else{
                        currXpos--;
                        xElem = Q[currXpos];
                        yElem = '-';
                    }
                }

            }
            _traceSteps++;
            _qOut[lengthCount] = xElem;
            _dOut[lengthCount] = yElem;
            lengthCount++;

        }

        _qOut[lengthCount] = '\0';
        _dOut[lengthCount] = '\0';

        _dOut = reverseArr(_dOut, strlen(_dOut));
        _qOut = reverseArr(_qOut, strlen(_qOut));

        /*
         * Write all the info demanded on an output file
         */
        //TODO ENABLE IT AT THE END
        fprintf(finFile, "\nMATCH %d [SCORE: %d,START: %ld,STOP: %ld]\n\tD: %s\n\tQ: %s\n", i+1, _MAX_SIMILARITY, (_endKeeper[i]-lengthCount) , _endKeeper[i]-1, _dOut, _qOut);

        free(_qOut);
        free(_dOut);

        //printf("\n\n");
    }

    double traceTimeFin = gettime();
    _totalTraceTime+=(traceTimeFin - traceTimeInit);




    //PRINTER OF THE  ARRAYS. USED FOR DEBBUGING PURPOSES ONLY
//    for (int i = 0; i < Q_len + 1; i++){
//        printf("\n");
//        for (int j = 0; j < D_len + 1; j++){
//            printf("%d\t",ScoreTable[i][j]);
//        }
//
//    }


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
    //printf("GO\n");

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
    //printf("RETURN FROM PARSING\n\n");
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
    finFileDir();

    FILE *fp;

    fp = fopen(INPUT,"r");

    if(fp == NULL){
        printf("Error opening file\n");
        exit(-9);
    }

    //TODO ENABLE
    finFile = fopen(cwd,"a");

    if(finFile == NULL){
        printf("Error while opening write file!\n");
        exit(1);
    }

    fileHeaderValues(fp);
    fileParser(fp);


    fclose(fp);
    //TODO ENABLE
    fclose(finFile);

    double timeFinTotal = gettime();
    _totalTime = timeFinTotal - timeInitTotal;

    terminalPrinter();

    return 0;

}


