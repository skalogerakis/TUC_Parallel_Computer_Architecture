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

char Q[];
char D[];
int MATCH;
int MISMATCH;
int GAP;
bool nameBool, inputBool, matchBool, misBool, gapBool = false;

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
            printf("NAME FOUND %s\n",argv[++i]);
            nameBool = true;
        }
        if(strcmp(argv[i],inVariable[1]) == 0){
            printf("INPUT FOUND %s\n",argv[++i]);
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


int main(int argc, char * argv[]) {


    commandChecker(argc, argv);

    printf("MATCH  %d, MISMATCH %d, GAP %d\n", MATCH, MISMATCH, GAP);

    //exit(1);

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

            if(tempMax > MAX_SIMILARITY){
                MAX_SIMILARITY = tempMax;
            }
        }
    }

    printf("SIMILARITY MAX%d \n", MAX_SIMILARITY);

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


