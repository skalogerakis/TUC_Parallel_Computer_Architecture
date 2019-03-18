#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <memory.h>
#include <stdbool.h>
//#include <w32api/rpcndr.h>

double gettime(void)
{
    struct timeval ttime;
    gettimeofday(&ttime, NULL);
    return ttime.tv_sec+ttime.tv_usec * 0.000001;
}

//char* Q, D;
int *similarity, *backtracking;

//int UPCODE = 1;
//int LEFTCODE = 2;
//int DIAGCODE = 3;

static char const UPCODE = 'U';
static char const LEFTCODE = 'L';
static char const DIAGCODE = 'D';
static char const ZEROCODE = 'Z';


int MAX_SIMILARITY=0;

char TestQ[] = "abc";
char TestD[] = "xxxabxcxxxaabbcc";

//char TestQ[] = "GGTTGACTA";
//char TestD[] = "TGTTACGG";

int Q_len, D_len;

int MATCH;
int MISMATCH;
int GAP;

bool nameBool, inputBool, matchBool, misBool, gapBool = false;


int pairs= -1;
int qMin = -1;
int qMax = -1;
int dSize = -1;

void testers(){

//    while(pairs>0){
//        fscanf(fp,)
//        pairs--;
//    }

//    char line[1000];
//    while((fscanf(fp, "%[^\n]", line)) != EOF){
//        fgetc(fp);
//        printf("Hello\n");
//        printf("Line = %s \n", line);
//    }
    char * line = NULL;
    //char * Q = (char *)malloc(sizeof(char)*qMax);



//    while(pairs>0){
//        //char  Q[1000];
//        //char D[1000];
//        char iniQ[qMax+10];
//        char iniD[dSize+10];
////        char *Q = (char *)malloc(sizeof(char)*(qMax+1));
////        char *D = (char *)malloc(sizeof(char)*(dSize+1));
//
//        fgets(iniQ, sizeof(iniQ) , fp);
//        printf("INI Q %d\n", strlen(iniQ));
////        if( sscanf(iniQ, "Q: %s\n", Q) == 1){
////            printf("Match Q %d\n", strlen(Q));
////        }
//        fgets(iniD, sizeof(iniD), fp);
//        printf("INI D %d\n", strlen(iniD));
////        if( sscanf(iniD, "D: %s\n", D) == 1){
////            printf("Match D %d\n", strlen(D));
////        }
//        //printf("Match D %s\n", D);
////        free(Q);
////        free(D);
//        free(iniD);
//        free(iniQ);
//        pairs--;
//        printf("\n\n\n");
//    }

//    ssize_t lineRead;
//    size_t lineLen =0;  //this will get reallocated
//    while ((lineRead = getline(&line, &lineLen, fp)) != -1) {
//        printf("Retrieved line of length %zu:\n", lineRead);
//
//        char  Q[100000];
//        char D[100000];
//        char tline[1000000];
//
//        //sscanf returns 1 when matches
//        if( sscanf(line, "Q: %s\n", Q) == 1){
//            fgets(tline, 1000000, fp);
//            //char *test = (char *)malloc(sizeof(char)*strlen(Q));
//            //strcpy(test,Q);
//            printf("Match Q %s\n", tline);
//            //Q = (char *)malloc(sizeof(char)*lineRead);
//
//        }else{
//            printf("NO match Q at %s\n", line);
//        }
//
//        if( sscanf(line, "D: %s\n", D) == 1){
//            printf("Match D %s\n", D);
//            //Q = (char *)malloc(sizeof(char)*lineRead);
//        }else{
//            printf("NO match D at %s\n", line);
//        }
//
//        //printf("%s\n", Q);
//        free(Q);
//    }
    //LOOKS WORKS
//    unsigned long count = 0;
//    size_t num;
//    unsigned long BUFELEM = 100;
//    unsigned char buf[BUFELEM];
//    do
//    {
//        num = fread( buf, 1, BUFELEM, fp );
//        ++count;
//        printf("%s buf %ld\n ",buf,count);
///* do sth. with buffer contents here [1] */
//    }
//    while ( num == BUFELEM );
}

/*
 * Add more input variables if you wish and update commandChecker function
 */
const char* inVariable[] = {"-name", "-input", "-match", "-mismatch", "-gap" };

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
    for (int i = 0; i < argc; i++) {
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

/*
 * This function is responsible to update all header file value.
 * In case the format is not correct then we exit the program
 */
void fileHeaderValues(FILE *fp){
    fscanf(fp, "Pairs: %d\n", &pairs);
    ErrorCode(pairs);
    //printf("Show pairs %d\n", pairs);

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

long fillDataBuffer(char * buf, long bytereader ,int compFlag){

    FILE *fwrite = fopen("D:\\TUC_PROJECT\\TUC_Parallel_Computer_Architecture\\MyDocs\\FINAL.txt","a");
    if(fwrite == NULL){
        printf("Error while opening write file!\n");
        exit(1);
    }

    long index = 0;
    long qIndex = 0;
    long dIndex = 0;
    //TODO change that
    //char Q[100000];
    char * Q = (char *)calloc(bytereader, sizeof(char));
    char * D = (char *)calloc(bytereader, sizeof(char));

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
                printf ( "Q %s\n", Q );
                printf ( "D is %s\n", D );
                fprintf(fwrite, "\n\nQ: \t%s\n", Q);
                fprintf(fwrite, "\n\nD: \t%s\n", D);
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


    //printf("%zu\n", bufLen);

//    if(bufLen<bytereader){
//        //printf ( "Q %s\n", Q );
//        //printf ( "D %s\n", D );
//        //free(Q);
//        //free(D);
//        return -1;
//    }


    printf ( "HI Q %s\n", Q );
    printf ( "HI D is %s\n", D );
    fprintf(fwrite, "\n\nQ: \t%s\n", Q);
    fprintf(fwrite, "\n\nD: \t%s\n", D);
    free(Q);
    free(D);
    return bufLen;

    //while ( buf != EOF ) {
//    for(index;index < bytereader; index++){
//
//        if ( buf[index] != '\n' ){
//            if(qFlag == 0){
//                Q[qIndex] = buf[index];
//                qIndex++;
//            }else{
//                D[dIndex] = buf[index];
//                dIndex++;
//            }
//
//        }else{
//            if(qFlag==0){
//                qFlag = 1;
//                printf ( "%s\n", Q );
//            }else{
//                printf ( "%s\n", D );
//                if(compFlag == 1) return -1;
//                return index;
//            }
//
//        }

//    }
//
   // return bytereader;
}

void fileParser(FILE *fp){
    long long MAXLINELENGTH = 100000000000000;
    int BUFSIZE =  dSize+ qMax +10;
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

    //TODO DONE
    commandChecker(argc, argv);

    printf("MATCH  %d, MISMATCH %d, GAP %d\n", MATCH, MISMATCH, GAP);

    FILE *fp;

    fp = fopen("D:\\TUC_PROJECT\\TUC_Parallel_Computer_Architecture\\MyDocs\\D6.txt","r");

    if(fp == NULL){
        printf("Error opening file\n");
        exit(-9);
    }
    //printf("Success\n");

    //TODO DONE
    fileHeaderValues(fp);
    fileParser(fp);

    fclose(fp);
    exit(1);


    Q_len = strlen(TestQ);
    D_len = strlen(TestD);
    printf("Length Q %d\n", Q_len);
    printf("Length D %d\n", D_len);



    //TODO try allocate 2d array

    int *ScoreTable[Q_len+1];
    for(int i = 0; i< Q_len+1; i++){
        ScoreTable[i] = (int *)malloc((D_len+1)* sizeof(int));
        if(ScoreTable[i] == NULL){
            printf("Could not allocate memory for matrix ScoreTable.Terminating");
            exit(-2);
        }
    }
    printf("Memory Allocated successfully for ScoreTable matrix.\n");

    //DONOT DELETE YET. INITIAL APPROACH WITH INT
//    int *TraceTable[Q_len+1];
//    for(int i = 0; i< Q_len+1; i++){
//        TraceTable[i] = (int *)malloc((D_len+1)* sizeof(int));
//        if(TraceTable[i] == NULL){
//            printf("Could not allocate memory for matrix TraceTable.Terminating");
//            exit(-2);
//        }
//    }

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
        //printf("Init check row %d\n",i);
    }

    printf("check done\n");


    for(int i=0; i < D_len + 1; i++){
        ScoreTable[0][i] = 0;
        //TraceTable[0][i] = 0;
        TraceTable[0][i] = 'Z';
        //printf("Init check column %d\n",i);
    }


    printf("START SIMILARITY MATRIX\n");

    int counterMax = 0;


    for(int i = 1; i < Q_len + 1; i++){
        for(int j = 1; j < D_len + 1; j++){
            //TableCalc(QDArray[i][j],Tr_QDArray[i][j],i,j);
            int up, left, diag;

            up = ScoreTable[i-1][j] + GAP;  //we want same column and one row up
            left = ScoreTable[i][j-1] + GAP; //we want same row and one column left

            //for diagonal check if we have match or mismatch
            uint tempDiag;

            //printf("Compare %c, %c\n",TestD[j - 1], TestQ[i - 1]);

            tempDiag = (TestD[j - 1] == TestQ[i - 1]) ? MATCH : MISMATCH;

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
            }

            //printf("check per loop %d \n", tempMax);
            ScoreTable[i][j] = tempMax;

            if(tempMax > MAX_SIMILARITY) {
                counterMax = 0;
                MAX_SIMILARITY = tempMax;
            }else if(tempMax == MAX_SIMILARITY && tempMax!=0){
                counterMax++;
            }
        }
    }


    printf("SIMILARITY MAX %d and COUNTER MAX %d\n", MAX_SIMILARITY, counterMax);

    /*
     * We find all max value we need to backtrack
     */
    int xMax[counterMax+1];
    int yMax[counterMax+1];

    int tempCount=0;
    for (int i = 0; i < Q_len + 1; i++){
        for (int j = 0; j < D_len + 1; j++){
            if(MAX_SIMILARITY == ScoreTable[i][j]){
                printf("i %d and j %d\n", i, j);

                xMax[tempCount] = i;
                yMax[tempCount] = j;
                tempCount++;
            }
        }
    }

    printf("SIMILARITY MAX %d and COUNTER MAX %d\n", MAX_SIMILARITY, counterMax);

    //char **qOut = (char **)malloc((counterMax+1)* sizeof(char *));

    for(int i = 0; i< counterMax+1; i++){

        int currXpos = xMax[i];
        int currYpos = yMax[i];
        char xElem = TestQ[currXpos];
        char yElem = TestD[currYpos];
        char currNode = TraceTable[currXpos][currYpos];
        int score = ScoreTable[xMax[i]][yMax[i]];
        int destCounter=0;
        //qOut[i] = (char *)malloc()
        char *qOut;


        //THIS SEEMS TO WORK TODO CHECK
        qOut = (char *)malloc(sizeof(char)*10000);
        int lengthCount=0;
        //char qOut[1000];
        //qOut[0] = 'Q';
        //strcpy(qOut,'Q');
        //if(xMax[i]>yMax[i]){
            //qOut[i] = (char *)malloc(sizeof(char)*xMax[i]);
        //}else{
            //qOut[i] = (char *)malloc(sizeof(char)*1000);
        //}

        while(currNode != ZEROCODE ){

            if(currNode == UPCODE){
                //printf("TRACE %d, - , TESTQ %c\n",ScoreTable[currXpos][currYpos], TestQ[currXpos]);
                currXpos--;
                xElem = TestQ[currXpos];
                yElem = '-';

                currNode = TraceTable[currXpos][currYpos];

            }else if(currNode == LEFTCODE){
                //printf("TRACE %d, TESTD %c, - \n",ScoreTable[currXpos][currYpos], TestD[currYpos]);
                currYpos--;
                xElem = '-';
                yElem = TestD[currYpos];

                currNode = TraceTable[currXpos][currYpos];

            }else if(currNode == DIAGCODE){
                //printf("TRACE %d, TESTD %c, TESTQ %c\n",ScoreTable[currXpos][currYpos], TestD[currYpos], TestQ[currXpos]);
                currYpos--;
                currXpos--;
                xElem = TestQ[currXpos];
                yElem = TestD[currYpos];

                currNode = TraceTable[currXpos][currYpos];

            }

            //printf("TRACE %d, TESTD %c, TESTQ %c\n",ScoreTable[currXpos][currYpos], TestD[currYpos], TestQ[currXpos]);
//            printf("TestD %c , TESTQ %c\n",yElem, xElem);
//            qOut[i][destCounter] = xElem;
//            destCounter++;
            qOut[lengthCount] = xElem;
            lengthCount++;
            printf("TestD %c , TESTQ %c\n",yElem, xElem);

        }
        free(qOut);
        //printf("TestD %c , TESTQ %c\n",TestD[currYpos], TestQ[currXpos]);
        //printf("TestD %c , TESTQ %c\n",yElem, xElem);
        printf("\n\n");
    }

    for (int i = 0; i < Q_len + 1; i++){
        printf("\n");
        for (int j = 0; j < D_len + 1; j++){
            printf("%d\t",ScoreTable[i][j]);
        }

    }

        printf("\n\n\n");
    for (int i = 0; i < Q_len + 1; i++){
        printf("\n");
        for (int j = 0; j < D_len + 1; j++){
            printf("%c\t",TraceTable[i][j]);
        }

    }

    //while

    //TODO add calc time in the end


//    double time0 = gettime();
//    printf("Hello, World!\n");
//    double time1 = gettime();
//
//    double elapsedTime = time1 - time0;
//
//
//    printf("Hello, World! %lf\n", &elapsedTime);
    return 0;
}


