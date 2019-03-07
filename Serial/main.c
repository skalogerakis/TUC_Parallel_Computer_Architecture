#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <memory.h>

double gettime(void)
{
    struct timeval ttime;
    gettimeofday(&ttime, NULL);
    return ttime.tv_sec+ttime.tv_usec * 0.000001;
}

char* Q, D;
int *similarity, *backtracking;

const char* TestQ = "abc";
const char* TestD = "xxxabxcxxxaabbcc";

void initializeMatrixes(int row_index, int column_index){

    for(int i = 0; i < row_index +1 ; i++){
        for(int j = 0; j < column_index + 1; j++){

        }
    }


}

void similarityCalc(int curr_index, int row_index, int column_index){

    //for initialization purposes only. We know that first column and row are filled with zero
    if(row_index == 0 || column_index == 0){
        similarity[curr_index] = 0;
        backtracking[curr_index] = 0;
        return;
    }


}

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


    size_t Q_len = strlen(TestQ);
    size_t D_len = strlen(TestD);

    printf("check %d %d\n", Q_len, D_len);



    //TODO check if needs that checks after malloc
    similarity = (int *)malloc( Q_len * D_len * sizeof(int));
    if(similarity == NULL){
        fprintf(stderr, "Fatal: failed to allocate %zu bytes.\n", Q_len * D_len * sizeof(int));
        abort();
    }

    backtracking = malloc( Q_len * D_len * sizeof(int));
    if(backtracking == NULL){
        fprintf(stderr, "Fatal: failed to allocate %zu bytes.\n", Q_len * D_len * sizeof(int));
        abort();
    }



    //printf("check %d\n", sizeof(similarity));

    //FILL SIMILARITY MATRIX

    for(int i = 0; i < D_len + 1; i++){
        for(int j = 0; j < Q_len + 1; j++){
            current_index++;
            similarityCalc(current_index, i,j);
        }
    }


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


