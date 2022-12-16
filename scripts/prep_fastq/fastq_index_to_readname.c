// Authors: Ben Parks

// Change log:
// 11/29/22 - Add --output-cmd option
// 9/27/22 - Created

// To compile: gcc -O3 -o fastq_index_to_readname fastq_index_to_readname.c

// Add index reads to the name field for fastqs output from bcl2fastq
// Input: R1, R2, I1, I2 files
// Output: R1, R2 files

// Output read names will be like:
// @A00509:608:HNT5JDSX3:1:2154:1063:1016 1:N:0:NCGGACGATCATG+TTTCNAGC
// Where the last sequence bits are I1 then I2

// It's possible to output compressed files by passing the last argument as
// --output-cmd 'gzip -c > $FILE'. This argument can be any command that will
// compress stdin to the path $FILE.

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#define BUF_SIZE 2048

// Used for temporary read/writes
char line_buf[BUF_SIZE];

// Read up to size bytes of the sequence field from a fastq file into str
void read_sequence(FILE *f, int size, char * str);

// Write a full fastq record from in to out, but add `prefix` then `index`
// to the end of the read name
void write_record(FILE *in, FILE *out, char* prefix, char* index);

int main(int argc, char **argv) {
    bool filter_output = false;
    FILE *R1_out, *R2_out;
    if (argc == 7) {
        filter_output = false;
        R1_out = fopen(argv[5], "w");
        R2_out = fopen(argv[6], "w");
    } else if (argc == 9 && strcmp(argv[7], "--output-cmd") == 0) {
        // Use popen to open the output files with a compression command
        if (setenv("FILE", argv[5], 1) != 0) {
            perror("ERROR");
            exit(1);
        }
        R1_out = popen(argv[8], "w");
        if (setenv("FILE", argv[6], 1) != 0) {
            perror("ERROR");
            exit(1);
        }
        R2_out = popen(argv[8], "w");
    } else {
        printf("Usage: %s R1_in R2_in I1_in I2_in R1_out R2_out [--output-cmd 'gzip -c > $FILE']\n", argv[0]);
        return 1;
    }
    
    char index_buf[BUF_SIZE];
    index_buf[BUF_SIZE - 1] = '*'; // sentinel
    line_buf[BUF_SIZE - 1] = '*'; // sentinel

    FILE *R1_in = fopen(argv[1], "r");
    FILE *R2_in = fopen(argv[2], "r");
    FILE *I1_in = fopen(argv[3], "r");
    FILE *I2_in = fopen(argv[4], "r");
    
    if (R1_in == NULL || R2_in == NULL || I1_in == NULL || I2_in == NULL || R1_out == NULL || R2_out == NULL) {
        fprintf(stderr, "ERROR: could not open all input and output files");
    }
    
    bool finished_reading = false;
    while (!feof(I1_in) && !feof(I2_in) && !feof(R1_in) && !feof(R2_in)) {
        // Read I1 + I2
        read_sequence(I1_in, BUF_SIZE, index_buf);
        int index_len = strlen(index_buf);
        index_buf[index_len - 1] = '+';
        read_sequence(I2_in, BUF_SIZE - index_len, index_buf + index_len);

        // Read + print R1
        write_record(R1_in, R1_out, "1:N:0:", index_buf);
        write_record(R2_in, R2_out, "2:N:0", index_buf);
        
        if (index_buf[BUF_SIZE - 1] == '\0' || line_buf[BUF_SIZE - 1] == '\0') {
            fprintf(stderr, "ERROR: fastq line length > %d characters\n", BUF_SIZE-1);
            exit(1);
        }

        if (ferror(R1_in) || ferror(R2_in) || ferror(I1_in) || ferror(I2_in) || ferror(R1_out) || ferror(R2_out)) {
            perror("ERROR");
            exit(1);
        }
    }

    // Check that all files finished reading at the same time
    if (!feof(I1_in) || !feof(I2_in) || !feof(R1_in) || !feof(R2_in)) {
        fprintf(stderr, "ERROR: input files not all of same length\n");
        exit(1);
    }

    // Clean up files
    fclose(I1_in); fclose(I2_in); fclose(R1_in); fclose(R2_in);
    if (filter_output) {
        if (pclose(R1_out) != 0 || pclose(R2_out) != 0) {
            fprintf(stderr, "ERROR: output filter commands had non-zero exit codes");
        }
    } else {
        fclose(R1_out); fclose(R2_out);
    }

    return 0;
}

// Read up to size bytes of the sequence field from a fastq file into str
void read_sequence(FILE *f, int size, char* index) {
    if (NULL == fgets(line_buf, BUF_SIZE, f)) return; // Skip, or return if nothing to read
    fgets(index, size, f); // Keep sequence
    fgets(line_buf, BUF_SIZE, f); // Skip
    if (feof(f)) {
        fprintf(stderr, "ERROR: truncated input file\n");
        exit(1);
    }
    fgets(line_buf, BUF_SIZE, f); // Skip
}

// Write a full fastq record from in to out, but add `prefix` then `index`
// to the end of the read name
void write_record(FILE *in, FILE *out, char* prefix, char* index) {
    // Parse first part of name
    if (NULL == fgets(line_buf, BUF_SIZE, in)) return; // return if nothing to read
    int name_length = 0;
    while(line_buf[name_length] != ' ' && line_buf[name_length] != '\0') {
        name_length += 1;
    }
    if (line_buf[name_length] == '\0') {
        fprintf(stderr, "ERROR: unexpected read name format - ' ' character not found\n");
        exit(1);
    }
    // Write first part of name
    fwrite(line_buf, sizeof(char), name_length + 1, out);
    // Write index sequence
    fputs(prefix, out);
    fputs(index, out);

    // Copy rest of the fastq lines
    fgets(line_buf, BUF_SIZE, in);
    fputs(line_buf, out);
    fgets(line_buf, BUF_SIZE, in);
    fputs(line_buf, out);
    if (feof(in)) {
        printf("ERROR: truncated input file\n");
        exit(1);
    }
    fgets(line_buf, BUF_SIZE, in);
    fputs(line_buf, out);
}