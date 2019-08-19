#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <errno.h>
#include <string.h>

#define MAXLENGTH 2048

/*
Check that the reads of the paired fastq files are sorted in paired order
and the only difference between paired read ids are '1'>'2' or '2'>'1' only.
*/
int fastq_read_pair_id_check(const char * filename1, const char * filename2){
    char line1[MAXLENGTH];
    char line2[MAXLENGTH];
    int result = 1;
    unsigned long int line_count1 = 0;
    unsigned long int line_count2 = 0;
    gzFile file1 = gzopen(filename1, "rb");
    gzFile file2 = gzopen(filename2, "rb");
    if (file1 == NULL || file2 == NULL){
        perror("File open failed");
    }
    while (gzgets(file1, line1, MAXLENGTH -1) && gzgets(file2, line2, MAXLENGTH -1) && result){
        if ((line_count1 % 4 == 0) && (line_count2 % 4 ==0)) {
            if (strlen(line1) != strlen(line2)){
                result = 0;
            } else {
                for (unsigned int i = 0; i < strlen(line1); i++){
                    if (line1[i] != line2[i]){
                        if (line1[i] == '1'){
                            if (line2[i] != '2'){
                                result = 0;
                            }
                        }else if (line1[i] == '2') {
                            if (line2[i] != '1'){
                                result = 0;
                            }
                        }
                    }
                }
            }
        }
        line_count1++;
        line_count2++;
    }
    // Count remaining lines
    while (gzgets(file1, line1, MAXLENGTH -1)){
        line_count1++;
    }
    while (gzgets(file2, line2, MAXLENGTH -1)){
        line_count2++;
    }
    if (line_count1 != line_count2){
        result = 0;
    }
    gzclose(file1);
    gzclose(file2);
    return result;
}

// flip the boolean output since 0 signal success for unix processes.
int main(int argc, char *argv[]){
    return !fastq_read_pair_id_check(argv[1], argv[2]);
}