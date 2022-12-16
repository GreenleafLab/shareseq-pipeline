// Authors: Ben Parks
// Last updated: 11/21/22

// To compile: gcc -O3 -o count_unique count_unique.c

// Deduplicate and count sorted lines.
// This is similar to `uniq -c`, but the count is added to the end of the line after a tab,
// and the command can also add together counts after merging sorted files.
//
// This can be used to count duplicates of ATAC fragments and RNA UMIs/genes
// 
// Input is sorted lines where duplicate entries will be on consecutive lines. If `merge`
// is passed as the only argument, each line will end in a tab character followed by a duplicate count.
// Output is deduplicated entries with summed duplicate counts

// Usage: 
//   - Deduplicate and add counts to the end of each line:
//     count_unique < input.tsv > output.tsv
//   - Deduplicate a sorted + merged output from dedup_fragments.py:
//     sort --merge input1.tsv input2.tsv | count_unique merge > output.tsv


// Reads from stdin, writes to stdout
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define BUF_SIZE 2048

int main(int argc, char **argv) {
    bool merge = false;
    if (argc == 2 && strcmp(argv[1], "merge") == 0) {
        merge = true;
    }
    if (argc != 1 && !merge) {
        fprintf(stderr, "Usage: %s [merge] < inputfile > outputfile\n", argv[0]);
        exit(1);
    }
    // Used for temporary read/writes
    char cur_line[BUF_SIZE];
    char prev_line[BUF_SIZE];
    int prev_count = 0;

    cur_line[BUF_SIZE - 1] = '*'; // sentinel
    prev_line[BUF_SIZE - 1] = '*'; // sentinel
    prev_line[0] = '\0';

    while (fgets(cur_line, BUF_SIZE, stdin)) {
        if (cur_line[BUF_SIZE - 1] == '\0') {
            fprintf(stderr, "ERROR: input line length > %d characters\n", BUF_SIZE-1);
            exit(1);
        }

        int count = 1;
        if (merge) {
            // Parse count number after tab
            char *dup_count = strrchr(cur_line, '\t');
            if (dup_count == NULL) {
                fprintf(stderr, "No tab character found in input line\n");
                exit(1);
            }
            count = atoi(dup_count + 1);
            *dup_count = '\0';
        } else {
            // Remove the newline from the string
            char *newline = strchr(cur_line, '\n');
            if (newline != NULL) {
                *newline = '\0';
            }
        }

        // Update state
        if (strcmp(prev_line, cur_line) == 0) {
            prev_count += count;
        } else {
            if (prev_line[0] != '\0') printf("%s\t%d\n", prev_line, prev_count);
            prev_count = count;
            strcpy(prev_line, cur_line);
        }
    }
    if (prev_line[0] != '\0') printf("%s\t%d\n", prev_line, prev_count);
}
