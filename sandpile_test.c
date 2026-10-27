#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int compare_files(const char *file_1, const char *file_2){
    FILE *fp1, *fp2;
    int c1, c2;
    int res = 1;
    if((fp1 = fopen(file_1, "rb")) != NULL && (fp2 = fopen(file_2, "rb")) != NULL){
        do{
            c1 = getc(fp1);
            c2 = getc(fp2);
            if(c1 != c2){
                res = 0;
                break;
            }
        }while(c1 != EOF && c2 != EOF);

        fclose(fp1);
        fclose(fp2);
    }

    return res;
}

int main(int argc, char *argv[]){
    if(argc != 3){
        fprintf(stderr, "usage: sandpile_test <file_name_1> <file_name_2>");
        return 1;
    }

    if(compare_files(argv[1], argv[2])){
        printf("Test passed!\n");
    }else{
        printf("Test failed!\n");
    }
    return 0;
}