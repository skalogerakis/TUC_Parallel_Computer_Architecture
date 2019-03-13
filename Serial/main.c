#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <memory.h>
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
int MATCH = 3;
int MISMATCH = -1;
int GAP = -1;


//void initializeMatrixes(int row_index, int column_index){
//
//    for(int i = 0; i < row_index +1 ; i++){
//        for(int j = 0; j < column_index + 1; j++){
//
//        }
//    }
//
//
//}
//
//void similarityCalc(int curr_index, int row_index, int column_index){
//
//    //for initialization purposes only. We know that first column and row are filled with zero
//    if(row_index == 0 || column_index == 0){
//        similarity[curr_index] = 0;
//        backtracking[curr_index] = 0;
//        return;
//    }
//

//}

//void TableCalc(int ScoreTable[D_len][Q_len] ,int TraceTable[D_len][Q_len],int rowIndex, int columnIndex){
//        printf("DIDI\n");
//        int up, left, diag;
//
//        up = ScoreTable[rowIndex-1][columnIndex] + GAP;  //we want same column and one row up
//        left = ScoreTable[rowIndex][columnIndex-1] + GAP; //we want same row and one column left
//
//        //for diagonal check if we have match or mismatch
//        uint tempDiag;
//
//        tempDiag = (TestD[rowIndex - 1] == TestQ[columnIndex - 1]) ? MATCH : MISMATCH;
//
//        diag = ScoreTable[rowIndex - 1][columnIndex - 1] + tempDiag;
//
//        uint tempMax = 0;
//
//        if( up <= 0 && left <= 0 && diag <= 0 ){    //case we have max as negative value
//            tempMax = 0;
//        }else{
//            if( up > left && up >= diag){
//                tempMax = up;
//                TraceTable[rowIndex][columnIndex] = UPCODE;
//            }else if (left >= up && left >= diag){
//                tempMax = left;
//                TraceTable[rowIndex][columnIndex] = LEFTCODE;
//            }else{
//                tempMax = diag;
//                TraceTable[rowIndex][columnIndex] = DIAGCODE;
//            }
//        }
//
//        ScoreTable[rowIndex][columnIndex] = tempMax;
//
//        if(tempMax > MAX_SIMILARITY){
//            MAX_SIMILARITY = tempMax;
//        }
//
//        return;
//
//
//
//}

int main() {

    /*TODO remove static values and append dynamically
     * Input variables:
     * -match
     * -mismatch
     * -gap
     */

    int match = 3;
    int mismatch = -1;
    int gap = -1;
    int current_index =0;




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
        printf("Init check row %d\n",i);
    }

    printf("check done\n");


    for(int i=0; i < D_len + 1; i++){
        ScoreTable[0][i] = 0;
        //TraceTable[0][i] = 0;
        TraceTable[0][i] = 'Z';
        printf("Init check column %d\n",i);
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


