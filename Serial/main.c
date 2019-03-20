#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <memory.h>
#include <stdbool.h>
//#include <w32api/rpcndr.h>


int _cellVal = 0;
int _traceSteps = 0;
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

int Q_len, D_len;

int MATCH;
int MISMATCH;
int GAP;

bool nameBool, inputBool, matchBool, misBool, gapBool = false;

int pairs= -1;
int qMin = -1;
int qMax = -1;
int dSize = -1;


FILE *finFile;

/*
 * Add more input variables if you wish and update commandChecker function
 */
const char* inVariable[] = {"-name", "-input", "-match", "-mismatch", "-gap" };

void terminalPrinter(){
    printf("\n\n");
    printf("NUMBER OF Q-D PAIRS: %d\n", pairs);
    printf("TOTAL NUMBER OF CELLS UPDATED: %d\n",_cellVal);
    printf("TRACEBACK STEPS: %d\n",_traceSteps);
    printf("TOTAL TIME : %lf\n", _totalTime);
    printf("TOTAL CELL TIME : %lf\n", _totalCellTime);
    printf("TOTAL TRACEBACK TIME : %lf\n", _totalTraceTime);
    printf("CUPS TOTAL TIME: %.3lf\n", _cellVal/_totalTime);
    printf("CUPS CELL TIME: %.3lf\n", _cellVal/_totalCellTime);

}

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
        //printf("Argument %s\n",argv[i]);
        if(strcmp(argv[i],inVariable[0]) == 0){
            //printf("NAME FOUND %s\n",argv[++i]);
            nameBool = true;
        }
        if(strcmp(argv[i],inVariable[1]) == 0){
            //printf("INPUT FOUND %s\n",argv[++i]);
            inputBool = true;
        }
        if(strcmp(argv[i],inVariable[2]) == 0) {
            //printf("MATCH FOUND %s\n", argv[++i]);
            MATCH = atoi(argv[++i]);
            matchBool = true;
        }
        if(strcmp(argv[i],inVariable[3]) == 0) {
           // printf("MISMATCH FOUND %s\n", argv[++i]);
            MISMATCH = atoi(argv[++i]);
            misBool = true;
        }
        if(strcmp(argv[i],inVariable[4]) == 0) {
            //printf("GAP FOUND %s\n", argv[++i]);
            GAP = atoi(argv[++i]);
            gapBool = true;
            //GAP = atoi(argv[++i]);
        }
    }

    if( !nameBool | !inputBool | !matchBool | !misBool | !gapBool){
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

char* reverseArr(char *str, size_t len) {
    size_t i = 0;
    while (len > i) {
        char tmp = str[--len];
        str[len] = str[i];
        str[i++] = tmp;
    }
    return str;
}
//TODO FIX THE PROBLEM WITH /001
//char * stringFixer(char *str) {
//    char *loc = strchr(str,'\001');
//    char *pDest = str;
//    if(loc != NULL){
//        while(*str){
//            if(*str != '\001')
//                *pDest++ = *str;
//            str++;
//        }
//        *pDest = '\0';
//    }
//    return *pDest;
//}

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

void dataParser(char * Q, char * D){

    //At first write as output the strings as they are
    fprintf(finFile, "\n\nQ: \t%s", Q);
    fprintf(finFile, "\nD: \t%s\n", D);

    Q_len = strlen(Q);
    D_len = strlen(D);

    int *ScoreTable[Q_len+1];
    for(int i = 0; i< Q_len+1; i++){
        ScoreTable[i] = (int *)malloc((D_len+1)* sizeof(int));
        if(ScoreTable[i] == NULL){
            printf("Could not allocate memory for matrix ScoreTable.Terminating");
            exit(-2);
        }
    }
    printf("Memory Allocated successfully for ScoreTable matrix.\n");


    char *TraceTable[Q_len+1];
    for(int i = 0; i< Q_len+1; i++){
        TraceTable[i] = (char *)malloc((D_len+1)* sizeof(char));
        if(TraceTable[i] == NULL){
            printf("Could not allocate memory for matrix TraceTable.Terminating");
            exit(-2);
        }
    }

    printf("Memory Allocated successfully for TraceTable matrix.\n");
    //Initialize 2D array (We have one additional column and row)

    for(int i=0; i < Q_len + 1; i++){
        ScoreTable[i][0] = 0;
        //TraceTable[i][0] = 0;
        TraceTable[i][0] = 'Z';
    }


    for(int i=0; i < D_len + 1; i++){
        ScoreTable[0][i] = 0;
        //TraceTable[0][i] = 0;
        TraceTable[0][i] = 'Z';
    }


    printf("START SIMILARITY MATRIX\n");

    int _counterMax = 0;

    double cellTimeInit = gettime();
    _MAX_SIMILARITY=0;

    for(int i = 1; i < Q_len + 1; i++){
        for(int j = 1; j < D_len + 1; j++){
            int up, left, diag;

            up = ScoreTable[i-1][j] + GAP;  //we want same column and one row up
            left = ScoreTable[i][j-1] + GAP; //we want same row and one column left

            //for diagonal check if we have match or mismatch
            uint tempDiag;

            tempDiag = (D[j - 1] == Q[i - 1]) ? MATCH : MISMATCH;

            diag = ScoreTable[i - 1][j - 1] + tempDiag;

            //printf("up %d, left %d, diag %d \n", up, left, diag);

            uint tempMax = 0;

            //TODO find smarter way
            if( up <= 0 && left <= 0 && diag <= 0 ){    //case we have max as negative value
                tempMax = 0;
                TraceTable[i][j] = ZEROCODE;
            }else{
                if( up > left && up >= diag){
                    tempMax = up;
                    TraceTable[i][j] = UPCODE;
                }else if (left >= up && left >= diag){
                    tempMax = left;
                    TraceTable[i][j] = LEFTCODE;
                }else{
                    tempMax = diag;
                    TraceTable[i][j] = DIAGCODE;
                }
                _cellVal++;
            }

            //printf("check per loop %d \n", tempMax);
            ScoreTable[i][j] = tempMax;

            if(tempMax > _MAX_SIMILARITY) {
                _counterMax = 0;
                _MAX_SIMILARITY = tempMax;
            }else if(tempMax == _MAX_SIMILARITY && tempMax!=0){
                _counterMax++;
            }
        }
    }

    for(int i = 0; i< Q_len+1; i++) {
        free(ScoreTable[i]);
    }

    double cellTimeFin = gettime();
    _totalCellTime+= (cellTimeFin-cellTimeInit);

    printf("Similarity matrix done. Start backtracking\n");

    /*
     * We find all max value we need to backtrack
     */
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

    //printf("SIMILARITY MAX %d and COUNTER MAX %d\n", _MAX_SIMILARITY, _counterMax);

    double traceTimeInit = gettime();
    for(int i = 0; i< _counterMax+1; i++){

        long currXpos = xMax[i];
        long currYpos = yMax[i];
        char xElem = Q[currXpos];
        char yElem = D[currYpos];
        char currNode = TraceTable[currXpos][currYpos];
        //char qOut[xMax[i]+1];
        //char dOut[yMax[i]+1];

        //TODO also works with dynamic allocation. Check if it can work with large files
        char *_qOut;

        //_qOut = (char *)calloc((xMax[i]+1)*(yMax[i]+1), sizeof(char));
        _qOut = (char *)malloc((xMax[i]+1)*(yMax[i]+1)*sizeof(char));
        if(_qOut==NULL){
            printf("Error occured while trying to allocate memory for traceback.Terminating....");
            exit(-1);
        }

        char *_dOut;
        //_dOut = (char *)calloc((xMax[i]+1)*(yMax[i]+1), sizeof(char));
        _dOut = (char *)malloc((xMax[i]+1)*(yMax[i]+1)* sizeof(char));
        if(_dOut==NULL){
            printf("Error occured while trying to allocate memory for traceback.Terminating....");
            exit(-1);
        }
        //char *qOut;

        //THIS SEEMS TO WORK TODO CHECK
        //qOut = (char *)malloc(sizeof(char) * (maxFinder(xMax[i],yMax[i])+1));
        int lengthCount = 0;

        while(currNode != ZEROCODE ){

            if(currNode == UPCODE){
                currXpos--;
                xElem = Q[currXpos];
                yElem = '-';

                currNode = TraceTable[currXpos][currYpos];

            }else if(currNode == LEFTCODE){
                currYpos--;
                xElem = '-';
                yElem = D[currYpos];

                currNode = TraceTable[currXpos][currYpos];

            }else if(currNode == DIAGCODE){
                currYpos--;
                currXpos--;
                xElem = Q[currXpos];
                yElem = D[currYpos];

                currNode = TraceTable[currXpos][currYpos];

            }
            _traceSteps++;
            //printf("TRACE %d, TESTD %c, TESTQ %c\n",ScoreTable[currXpos][currYpos], TestD[currYpos], TestQ[currXpos]);
//            printf("TestD %c , TESTQ %c\n",yElem, xElem);
//            qOut[i][destCounter] = xElem;
//            destCounter++;
            _qOut[lengthCount] = xElem;
            _dOut[lengthCount] = yElem;

            lengthCount++;

        }


        //TODO check what is going on here
//        if(strlen(_dOut)!= strlen(_qOut)){
//            if(strlen(_dOut) < strlen(_qOut) ){
//                _qOut = (char *)stringFixer(_qOut);
//            }else{
//                _dOut = (char *)stringFixer(_dOut);
//            }
//        }
        _dOut = reverseArr(_dOut, strlen(_dOut));
        _qOut = reverseArr(_qOut, strlen(_qOut));


        //CHECKED OK
        //Here we are writing the other demanded info(score, start, stop)
        fprintf(finFile, "\nMATCH %d [SCORE: %d,START: %d,STOP: %d]\n\tD: %s\n\tQ: %s\n", i+1, _MAX_SIMILARITY, (_endKeeper[i]-lengthCount) , _endKeeper[i]-1, _dOut, _qOut);

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
//
//    printf("\n\n\n");
//    for (int i = 0; i < Q_len + 1; i++){
//        printf("\n");
//        for (int j = 0; j < D_len + 1; j++){
//            printf("%c\t",TraceTable[i][j]);
//        }
//
//    }


}

//TODO finishing touches and check what is needed and not
long fillDataBuffer(char * buf, long bytereader ,int compFlag){
    //Needed for testing. Will use it later.
//    FILE *fwrite = fopen("D:\\TUC_PROJECT\\TUC_Parallel_Computer_Architecture\\MyDocs\\FINAL.txt","a");
//    if(fwrite == NULL){
//        printf("Error while opening write file!\n");
//        exit(1);
//    }

    long index = 0;
    long qIndex = 0;
    long dIndex = 0;
    //TODO check that memory was allocated successfully

    char * Q;
    char * D;

    Q = (char *)calloc(bytereader, sizeof(char));
    if(Q==NULL){
        printf("Error occured while trying to allocate memory for buffer.Terminating....");
        exit(-1);
    }

    D = (char *)calloc(bytereader, sizeof(char));
    if(D==NULL){
        printf("Error occured while trying to allocate memory for buffer.Terminating....");
        exit(-1);
    }
    //char D[100000];
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
                //printf ( "Q inedx %d\n", index );
                //printf ( "Q %s\n", Q );
                //printf ( "D is %s\n", D );
                dataParser(Q, D);
                //fprintf(fwrite, "\n\nQ: \t%s\n", Q);
                //fprintf(fwrite, "\n\nD: \t%s\n", D);
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


    //printf ( "HI Q %s\n", Q );
    //printf ( "HI D is %s\n", D );
    dataParser(Q, D);
//    fprintf(fwrite, "\n\nQ: \t%s\n", Q);
//    fprintf(fwrite, "\n\nD: \t%s\n", D);
    free(Q);
    free(D);
    return bufLen;
}

//TODO finishing touches and check what is needed and not
void fileParser(FILE *fp){
    long long MAXLINELENGTH = 100000000000000;
    int BUFSIZE =  dSize+ qMax + 100;
    long            bytesread;
    char            buf[BUFSIZE];
    int              sizeLeftover=0;
    int              bLoopCompleted = 0;
    long        pos = 0;
    int pairCounter = 0;

    do
    {

        // Read next block from file and save into buf, right after the
        // "left over" buffer
        bytesread = fread(buf+sizeLeftover, 1, sizeof(buf)-1-sizeLeftover, fp);
        if (bytesread<1)
        {
            // Turn on 'loop completed' flag so that we know to exit at the bottom
            // Still need to process any block we currently have in the
            // leftover buffer
            bLoopCompleted = 1;
            bytesread  = 0;
        }

        // Add NULL terminator at the end of our buffer
        buf[bytesread+sizeLeftover] = 0;

        // Process data - Replace with your function
        //
        // Function should return the position in the file or -1 if failed
        //
        // We are also passing bLoopCompleted to let ProcessData know whether this is
        // the last record (in which case - if no end-of-record separator,
        // use eof and process anyway)

        pos = fillDataBuffer(buf, bytesread+sizeLeftover,
                          bLoopCompleted);
        pairCounter++;
        // If error occured, bail
        if (pos<1)
        {
            bLoopCompleted = 1;
            pos      = 0;
        }

        // Set Left over buffer size to
        //
        //  * The remaining unprocessed buffer that was not processed
        //  by ProcessData (because it couldn't find end-of-line)
        //
        // For protection if the remaining unprocessed buffer is too big
        // to leave sufficient room for a new line (MAXLINELENGTH), cap it
        // at maximumsize - MAXLINELENGTH
        //sizeLeftover = mymin(bytesread+sizeLeftover-pos, sizeof(buf)-MAXLINELENGTH);
        if(bytesread+sizeLeftover-pos > sizeof(buf)-MAXLINELENGTH){
            sizeLeftover = sizeof(buf)-MAXLINELENGTH;
        }else{
            sizeLeftover = bytesread+sizeLeftover-pos;
        }
        // Extra protection - should never happen but you can never be too safe
        if (sizeLeftover<1) sizeLeftover=0;

        // If we have a leftover unprocessed buffer, move it to the beginning of
        // read buffer so that when reading the next block, it will connect to the
        // current leftover and together complete a full readable line
        if (pos!=0 && sizeLeftover!=0)
            memmove(buf, buf+pos, sizeLeftover);

    } while(!bLoopCompleted && pairCounter != pairs);
}


int main(int argc, char * argv[]) {

    double timeInitTotal = gettime();

    commandChecker(argc, argv);

    //printf("MATCH  %d, MISMATCH %d, GAP %d\n", MATCH, MISMATCH, GAP);
    printf("STARTING EXECUTION....\n");

    FILE *fp;

    fp = fopen("D:\\TUC_PROJECT\\TUC_Parallel_Computer_Architecture\\MyDocs\\D8.txt","r");

    if(fp == NULL){
        printf("Error opening file\n");
        exit(-9);
    }

    finFile = fopen("D:\\TUC_PROJECT\\TUC_Parallel_Computer_Architecture\\MyDocs\\FINAL.txt","a");

    if(finFile == NULL){
        printf("Error while opening write file!\n");
        exit(1);
    }

    fileHeaderValues(fp);
    fileParser(fp);


    fclose(fp);
    fclose(finFile);

    double timeFinTotal = gettime();
    _totalTime = timeFinTotal - timeInitTotal;

    terminalPrinter();

    return 0;

}


