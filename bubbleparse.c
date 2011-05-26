/*----------------------------------------------------------------------*
 * File:    bubbleparse.c                                               *
 *                                                                      *
 * Purpose: Parse files output by cortex bubble finding and produce     *
 *          ranked list.                                                *
 *                                                                      *
 * Author:  Richard Leggett                                             *
 *          The Genome Analysis Centre                                  *
 *          Norwich Research Park, Colney, Norwich, NR4 7UH, UK         *
 *          richard.leggett@bbsrc.ac.uk                                 * 
 *                                                                      *
 * History: 06-Jul-10: RML: Created                                     *
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <execinfo.h>
#include <signal.h>
#include <inttypes.h>
#include "binary_kmer.h"
#include "flags.h"
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "logger.h"

/*----------------------------------------------------------------------*
 * Constants                                                            *
 *----------------------------------------------------------------------*/
#define MAX_PATHS_PER_MATCH 8 
#define MAX_COLOURS 2
#define TRUE 1
#define FALSE 0

/*----------------------------------------------------------------------*
 * Flags                                                                *
 *----------------------------------------------------------------------*/
#define FLAG_IGNORE                    1
#define FLAG_NO_MATCH_LOW              2
#define FLAG_NO_MATCH_HIGH             4
#define FLAG_BELOW_MINIMUM_CONTIG_SIZE 8
#define FLAG_REPEAT                    16

#define FLAG_PATHS_SAME_LENGTH         1
#define FLAG_COVERAGE_WITHIN_BOUNDS    2

/*----------------------------------------------------------------------*
 * Structures to hold data                                              *
 *----------------------------------------------------------------------*/
typedef struct {
    int* coverage;
} CoverageArray;

typedef struct {
    char type[16];
    long double cpnp;
    short flags;
    int q_total;
    long double p_ratio[2];
    long double p_value[MAX_PATHS_PER_MATCH][MAX_COLOURS];
    short coverage_complete[MAX_PATHS_PER_MATCH][MAX_COLOURS];
    double coverage_av[MAX_PATHS_PER_MATCH][MAX_COLOURS];
    double coverage_pc[MAX_PATHS_PER_MATCH][MAX_COLOURS];
    double coverage_difference[MAX_COLOURS];
    double combined_coverage;
} RankingStats;

typedef struct {
    long contig_file_offset;    
    long coverage_file_offset;
    char* contig_pre;
    char* contig_mid;
    char* contig_post;
    int contig_length;
    int pre;
    int mid;
    int post;
    short flags;
    short first_bubble_kmer_offset;
    char *first_bubble_kmer;
    char *first_contig_kmer;
    char *last_contig_kmer;
    CoverageArray* colour_coverage[MAX_COLOURS];
} MatchPath;

struct _Match {
    int match_number;
    int number_of_paths;
    int longest_contig;
    short flags;
    //struct _Match* linked_match;
    RankingStats statistics;
    MatchPath* paths[MAX_PATHS_PER_MATCH];
};

typedef struct _Match Match;

typedef struct {
    int length;
    int pre;
    int mid;
    int post;
    long file_offset;
    char* seq;
} PreprocessedPath;

typedef struct {
    int match_number;
    int number_of_paths;
    PreprocessedPath paths[MAX_PATHS_PER_MATCH];
} PreprocessedMatch;

/*----------------------------------------------------------------------*
 * Global variables                                                     *
 *----------------------------------------------------------------------*/
int max_matches = 300000;
char fasta_filename[4096];
char coverage_filename[4096];
Match** matches = 0;
int number_of_matches = 0;
int number_of_colours = 2;
int kmer_size = 0;
int no_match_low = 0;
int no_match_high = 0;
double expected_coverage_pc[MAX_PATHS_PER_MATCH][MAX_COLOURS];
double coverage_pc_tolerance[MAX_COLOURS];
HashTable * hash_table = NULL;
char* base_filename;
char* file_of_filenames;
char* expected_filename = 0;
char* options_filename = 0;
char* match_selection_filename = 0;
char* rank_table_filename = 0;
char* rank_contig_filename = 0;
char* rank_csv_filename = 0;
char* rank_log_filename = 0;
char* debug_log_filename = 0;
char selected_matches_filename[4096];
int minimum_contig_size = 0;
int quality_offset = 64;
int n_length_adjusted = 0;
int n_kmers_already = 0;
int deduplicate = 0;
int use_quality_scores = 0;
int n_deduplicates = 0;
int max_line_length = 8*1024*1024;

/*----------------------------------------------------------------------*
 * Table column widths                                                  *
 *----------------------------------------------------------------------*/
//                     0  1  2  3   4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24  25  26  27  28  29  30 
int column_widths[] = {6, 6, 3, 10, 5, 3, 3, 3, 6, 8, 8, 8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 12, 10, 10, 10, 10, 10, 10};
int space_between_columns = 1;
#ifdef EXTRA_STATS
int number_of_columns = 31;
#else
int number_of_columns = 24;
#endif

/*----------------------------------------------------------------------*
 * Flags                                                                *
 *----------------------------------------------------------------------*/
void set_flag(short* s, short f)
{
    *s |= f;
}

void unset_flag(short*s, short f)
{
    *s ^= f;
}

int flag_is_set(short s, short f)
{
    return s & f;
}

int flag_is_unset(short s, short f)
{
    return !(s & f);
}

/*----------------------------------------------------------------------*
 * Function: make_reverse_compliment                                    *
 * Purpose:  make a reverse compliment string                           *
 * Params:   in -> input string                                         *
 *           out -> buffer for reverse compliment string                *
 * Returns:  Pointer to output string                                   *
 *----------------------------------------------------------------------*/
char* make_reverse_compliment(char* in, char* out)
{
    int i;
    for (i=0; i<strlen(in); i++) {
        switch(in[strlen(in)-1-i]) {
            case 'T':
                out[i] = 'A';
                break;
            case 'G':
                out[i] = 'C';
                break;
            case 'C':
                out[i] = 'G';
                break;
            case 'A':
                out[i] = 'T';
                break;
        }
    }
    
    out[strlen(in)] = 0;
    
    return out;
}

/*----------------------------------------------------------------------*
 * Function: make_time_string                                           *
 * Purpose:  make a time string of the form hh:mm:ss.                   *
 * Params:   s -> where to store string                                 *
 * Returns:  Pointer to string                                          *
 *----------------------------------------------------------------------*/
char* make_time_string(char* s)
{
    time_t t = time(NULL);
    struct tm* loctime = localtime(&t);
    sprintf(s, "%2d:%02d:%02d", loctime->tm_hour, loctime->tm_min, loctime->tm_sec);
    return s;
}

/*----------------------------------------------------------------------*
 * Function: clean_line                                                 *
 * Purpose:  remove redundant non-printable characters from the end of  *
 *           lines.                                                     *
 * Params:   line -> pointer to line to change.                         *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void clean_line(char* line)
{
    int i;
    for (i=strlen(line)-1; i>=0; i--) {
        if (line[i] <= ' ') {
            line[i] = 0;
        } else {
            break;
        }
    }
}

char* make_upper_case(char* line)
{
    int i;
    for (i=0; i<strlen(line); i++) {
        line[i] = toupper(line[i]);
    }
    
    return line;
}

char* make_lower_case_copy(char* copy, char* original)
{
    int i;
    for (i=0; i<strlen(original); i++) {
        copy[i] = tolower(original[i]);
    }
    copy[i] = 0;
    
    return copy;
}

void strcpy_upper(char* dest, char* src, int n) {
    int i=0;
    
    while ((i<strlen(src)) && (i<n)) {
        dest[i] = toupper(src[i]);
        i++;
    }
    dest[i] = 0;
}

/*----------------------------------------------------------------------*
 * Function: clean_and_check_sequence                                   *
 * Purpose:  remove non-printable terminating characters and change all *
 *           to upper case.                                             *
 * Params:   line -> pointer to line to change.                         *
 * Returns:  0 for failure, 1 for success                               *
 *----------------------------------------------------------------------*/
int clean_and_check_sequence(char* line)
{
    int i;
    int ok = 1;
    
    // Remove any problem characters at the end of line - eg. CR, LF
    clean_line(line);
    
    // Convert case and check
    for (i=0; i<strlen(line); i++) {
        // Convert case
        if ((line[i] >= 'a') && (line[i] <= 'z')) {
            line[i] = toupper(line[i]);
        }
        
        // Check T, G, A or C
        if ((line[i] != 'T') && (line[i] != 'G') && (line[i] != 'A') && (line[i] != 'C')) {
            ok = 0;
            break;
        }
    }
    
    return ok;
}

/*----------------------------------------------------------------------*
 * Function: get_next_line                                              *
 * Purpose:  Get next line of length > 1 from file.                     *
 * Params:   line -> pointer to character buffer to store line          *
 *           max_size = maximum length of buffer                        *
 *           fp -> file pointer                                         *
 * Returns:  Pointer to line                                            *
 *----------------------------------------------------------------------*/
char* get_next_line(char* line, int max_size, FILE* fp)
{
    char* r = 0;
    
    while (!r) {
        if (!fgets(line, max_size, fp)) {
            break;
        } else {
            if (strlen(line) > 1) {
                r = line;
            }
        }
    }
    
    return r;
}

/*----------------------------------------------------------------------*
 * Function: get_first_and_last_kmers                                   *
 * Purpose:  Get and store the first kmer in the bubble and the first   *
 *           and last kmers from a path                                 *
 * Params:   line -> string containing contig                           *
 *           m = match number                                           *
 *           p = path number                                            *
 *           pre = size of prefix flanking                              *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void get_first_and_last_kmers(char* line, int m, int p)
{
    int pre = matches[m]->paths[p]->pre;
    int k_offset = pre - kmer_size + 1;
    int i;
    
    // First 
    if (k_offset < 0) {
        printf("Error: match %i path %i, pre too small (%i). Is kmer size correct?\n", m, p, pre);
        exit(1);
    }
    
    if (strlen(line) < (k_offset + kmer_size)) {
        printf("Error: match %i path %i, line too short (%i)\n", m, p, k_offset + kmer_size);
        exit(1);
    }
    
    matches[m]->paths[p]->first_bubble_kmer = calloc(kmer_size+1, sizeof(char));
    matches[m]->paths[p]->first_contig_kmer = calloc(kmer_size+1, sizeof(char));
    matches[m]->paths[p]->last_contig_kmer = calloc(kmer_size+1, sizeof(char));

   if (!matches[m]->paths[p]->first_bubble_kmer ||
       !matches[m]->paths[p]->first_contig_kmer ||
       !matches[m]->paths[p]->last_contig_kmer) {
        printf("Error: can't get memory to store first/last kmer\n");
        exit(1);
    }
    
    strncpy(matches[m]->paths[p]->first_bubble_kmer, line + k_offset, kmer_size);
    matches[m]->paths[p]->first_bubble_kmer[kmer_size] = 0;                 
    matches[m]->paths[p]->first_bubble_kmer_offset = kmer_size-1;
    
    for (i=0; i<strlen(matches[m]->paths[p]->first_bubble_kmer); i++) {
        char c = toupper(matches[m]->paths[p]->first_bubble_kmer[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            printf("Error: bad character in match %i path %i (%c)\n", m, p, c);
            exit(1);
        }
    }
        
    strncpy(matches[m]->paths[p]->first_contig_kmer, line, kmer_size);
    strncpy(matches[m]->paths[p]->last_contig_kmer, line + strlen(line) - kmer_size, kmer_size);
    matches[m]->paths[p]->first_contig_kmer[kmer_size] = 0;
    matches[m]->paths[p]->last_contig_kmer[kmer_size] = 0;
}

/*----------------------------------------------------------------------*
 * Function: preprocess_match                                           *
 * Purpose:  Preprocesses a match structure - currently, this means     *
 *           verifying/correcting the flanking lengths, as there is     *
 *           a cortex_bub problem that results in these occasionally    *
 *           being out by 1. In the future, this will be corrected.     *
 * Params:   match -> a PreprocessedMatch to check                      *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void preprocess_match(PreprocessedMatch *match)
{
    int p;
    int n_prefix;
    int n_suffix;
    char current_kmer[kmer_size+1];
    char c;
    int found_difference;
    int reached_end;
    
    n_prefix = -1;
    found_difference = FALSE;
    reached_end = FALSE;
    do {
        n_prefix++;
        c = 0;
        for (p=0; p<match->number_of_paths; p++) {
            if (n_prefix >= strlen(match->paths[p].seq)) {
                reached_end = TRUE;
            } else if (c == 0) {
                c = match->paths[p].seq[n_prefix];
            } else if (match->paths[p].seq[n_prefix] != c) {
                found_difference = TRUE;
            }
        }
    } while (!found_difference && !reached_end);

    if (reached_end) {
        printf("Error: Reached end of string while preprocessing match %d.\n", match->match_number);
        exit(1);
    }    
    
    n_suffix = -1;
    found_difference = FALSE;
    reached_end = FALSE;
    do {
        n_suffix++;
        current_kmer[0] = 0;
        for (p=0; p<match->number_of_paths; p++) {
            int i_suffix = strlen(match->paths[p].seq) - n_suffix - kmer_size;
            if (n_suffix >= strlen(match->paths[p].seq)) {
                reached_end = TRUE;
            } else if (current_kmer[0] == 0) {
                strncpy(current_kmer, match->paths[p].seq+i_suffix, kmer_size);
                current_kmer[kmer_size] = 0;
            } else {
                if (strncmp(match->paths[p].seq+i_suffix, current_kmer, kmer_size) != 0) {
                    found_difference = TRUE;
                }
            }
        }
    } while (!found_difference && !reached_end);
    
    if (reached_end) {
        printf("Error: Reached end of string while preprocessing match %d.\n", match->match_number);
        exit(1);
    }
    
    for (p=0; p<match->number_of_paths; p++) {                
        if ((match->paths[p].pre != n_prefix) || (match->paths[p].post != n_suffix)) {
            log_printf("Warning: Flanking ");

            if ((match->paths[p].pre != n_prefix) && (match->paths[p].post != n_suffix)) {
                log_printf("(prefix and suffix) ");
            } else if (match->paths[p].pre != n_prefix) {
                log_printf("(prefix) ");
            } else if (match->paths[p].post != n_suffix) {
                log_printf("(suffix) ");
            }
            
            log_printf("length differs (%d, %d, %d, %d) ",
                       match->paths[p].length, match->paths[p].pre, match->paths[p].mid, match->paths[p].post);
                       
            match->paths[p].pre = n_prefix;
            match->paths[p].post = n_suffix;
            match->paths[p].mid = match->paths[p].length - n_prefix - n_suffix;

            log_printf("to calculated (%d, %d, %d, %d) for match %d path %d. Correcting.\n",
                       match->paths[p].length, match->paths[p].pre, match->paths[p].mid, match->paths[p].post, match->match_number, p);
            
            if (match->paths[p].mid < 0) {
                printf("Error: calculated mid length wrong for match %d path %d. Is kmer size correct?\n", match->match_number, p);
                exit(1);
            }
            
            n_length_adjusted++;
        }
    }
}

/*----------------------------------------------------------------------*
 * Function: read_next_match                                            *
 * Purpose:  Read the next match from the FASTA file and put into a     *
 *           PreprocessedMatch structure.                               *
 * Params:   fp -> file pointer                                         *
 *           match -> a PreprocessedMatch to check                      *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void read_next_match(FILE *fp, PreprocessedMatch *match)
{
    long position;
    char line[1024];
    char* s;
    int m;
    int p;
    int l;
    int end_of_match = 0;
    int ignore_path;

    match->match_number = -1;
    match->number_of_paths = 0;
    
    while (!feof(fp) && (!end_of_match)) {
        position = ftell(fp);
        if (fgets(line, 1024, fp)) {
            if (line[0] == '>') {
                // Reset ignore_path
                ignore_path = 0;
                
                // Read match and path number
                sscanf(line, ">match_%d_path_%d length:%d ", &m, &p, &l);
                
                // Check if we've got a new match
                if (match->match_number == -1) {
                    match->match_number = m;
                } else if (m != match->match_number) {
                    end_of_match = 1;
                    fseek(fp, position, SEEK_SET);
                }
                
                // Check path number (indexed from 0)
                if (!end_of_match) {
                    if (p >= MAX_PATHS_PER_MATCH) {
                        log_printf("Warning: Match %d path index %d too high (Max %d paths allowed). Path ignored.\n", m, p, MAX_PATHS_PER_MATCH);
                        ignore_path = 1;
                    } else if (p != match->number_of_paths) {
                        log_and_screen_printf("Error: nonsequential path number (%d found, expected %d) on match %d.\n", p, match->number_of_paths, m);
                        fclose(fp);
                        exit(1);
                    }
                }
                                
                if ((!end_of_match) && (!ignore_path)) {
                    // Store length
                    match->paths[p].length = l;
                    match->number_of_paths = p+1;
                    match->paths[p].seq = malloc(l+1);
                    
                    // Read pre, mid and post lengths
                    s = strstr(line, "pre_length");
                    if (s) {
                        sscanf(s, "pre_length:%d mid_length:%d post_length:%d", &(match->paths[p].pre), &(match->paths[p].mid), &(match->paths[p].post)); 
                        match->paths[p].file_offset = ftell(fp);
                        if (!fgets(match->paths[p].seq, l+1, fp)) {
                            log_and_screen_printf("Error: couldn't read sequence for match %d path %d.\n", m, p);
                            fclose(fp);
                            exit(1);
                        }

                        if (!clean_and_check_sequence(match->paths[p].seq)) {
                            log_and_screen_printf("Error: Strange contig data for match %d path %d.\n", m, p);
                            fclose(fp);
                            exit(1);
                        }                        
                    } else {
                        log_and_screen_printf("Error: couldn't extract header fields for match %d path %d.\n", m, p);
                        fclose(fp);
                        exit(1);
                    }
                }
            }
        }
    }
    
    if (match->match_number != -1) {
        preprocess_match(match);
    }
}

/*----------------------------------------------------------------------*
 * Function: read_fasta_file                                            *
 * Purpose:  Read the FASTA file and store data.                        *
 * Params:   filename -> filename of FASTA file                         *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void read_fasta_file(char* filename)
{
    FILE* fp;
    int m, p;
    PreprocessedMatch ppmatch;
    
    log_newline();
    log_write_timestamp(0);
    log_printf(" START Read FASTA file\n");
    
    printf("Reading FASTA file...\n");
    
    fp = fopen(filename, "r");
    if (!fp) {
        printf("Error: Can't open FASTA file.\n");
        exit(1);
    }
    
    while (!feof(fp)) {
        read_next_match(fp, &ppmatch);
        m = ppmatch.match_number;

        if (m != -1) {
            if (m != number_of_matches) {
                printf("Error: Non-sequential match numbers: Found %d expected %d.\n", m, number_of_matches);
                fclose(fp);
                exit(1);
            }
            
            // Do we need to allocate more memory to store matches?
            if (m > (max_matches-1)) {              
                int new_max_matches = max_matches + 1000;
                Match** new_matches = realloc(matches, new_max_matches*sizeof(Match*));
                int i;
                
                if (!new_matches) {
                    printf("Error: Out of memory trying to store match %d.\n", m);
                    fclose(fp);
                    exit(1);
                }
                
                for (i=max_matches; i<new_max_matches; i++) {
                    new_matches[i] = 0;
                }
                
                matches = new_matches;
                max_matches = new_max_matches;
            }
            
            // Get space for match
            if (matches[m] == 0) {
                int i;
                
                matches[m] = calloc(1, sizeof(Match));
                
                if (!matches[m]) {
                    printf("Error: couldn't get memory for Match structure.\n");
                    fclose(fp);
                    exit(1);
                }
                
                for (i=0; i<MAX_PATHS_PER_MATCH; i++) {
                    matches[m]->paths[i] = 0;
                }
            }

            number_of_matches = m+1;            
            matches[m]->match_number = m;
            matches[m]->number_of_paths = ppmatch.number_of_paths;
            matches[m]->flags = 0;
                        
            // Go through paths
            for (p=0; p<ppmatch.number_of_paths; p++) {
                int pre = ppmatch.paths[p].pre;
                int mid = ppmatch.paths[p].mid;
                int post = ppmatch.paths[p].post;
                int l = ppmatch.paths[p].length;
                char* line = ppmatch.paths[p].seq;
                
                // Get space for path
                if (matches[m]->paths[p] == 0) {
                    matches[m]->paths[p] = malloc(sizeof(MatchPath));
                    if (!matches[m]->paths[p]) {
                        printf("Error: couldn't get memory for MatchPath structure.\n");
                        fclose(fp);
                        exit(1);
                    }
                }
                        
                #ifndef DEBUG_MAC
                if (l != (pre+mid+post)) {
                    printf("  Match %d path %d: pre, mid, post (%d, %d, %d) don't add up to length (%d)\n", m, p, pre, mid, post, l);
                }                   
                #endif

                matches[m]->paths[p]->contig_length = l;
                matches[m]->paths[p]->contig_file_offset = ppmatch.paths[p].file_offset;
                matches[m]->paths[p]->pre = pre;
                matches[m]->paths[p]->mid = mid;
                matches[m]->paths[p]->post = post;
                matches[m]->paths[p]->flags = 0;                                
                                
                matches[m]->paths[p]->contig_pre = malloc(pre+1);
                matches[m]->paths[p]->contig_mid = malloc(mid+1);
                matches[m]->paths[p]->contig_post = malloc(post+1);
                if ((!matches[m]->paths[p]->contig_pre) ||
                    (!matches[m]->paths[p]->contig_mid) ||
                    (!matches[m]->paths[p]->contig_post)) {
                    printf("Error: couldn't allocate memory for contig (lengths %d, %d, %d)!\n", pre, mid, post);
                    exit(1);
                }
                strcpy_upper(matches[m]->paths[p]->contig_pre, line, pre);                  
                strcpy_upper(matches[m]->paths[p]->contig_mid, line+pre, mid);                  
                strcpy_upper(matches[m]->paths[p]->contig_post, line+pre+mid, post);                    

                get_first_and_last_kmers(line, m, p);
            }            
        }
    }
    
    fclose(fp);

    log_printf("%d paths length adjusted.\n", n_length_adjusted);

    log_write_timestamp(0);
    log_printf(" ENDED Read FASTA file\n");
}


#ifdef OLD_READ_FASTA_FILE
/*----------------------------------------------------------------------*
 * Function: read_fasta_file                                            *
 * Purpose:  Read the FASTA file and store data.                        *
 * Params:   filename -> filename of FASTA file                         *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void read_fasta_file(char* filename)
{
    FILE* fp;
    char* s;
    int m,p,l;
    int pre, mid, post;
    long pos;
    int current_match = -1;
    int current_path = -1;
    int ignore_path = 0;
    char *line;
        
    log_newline();
    log_write_timestamp(0);
    log_printf(" START Read FASTA file\n");

    printf("Reading FASTA file...\n");

    line = malloc(max_line_length);
    if (!line) {
        printf("Error: Insufficient memory to allocate space to store line.\n");
        exit(1);
    }   
    
    fp = fopen(filename, "r");
    if (!fp) {
        printf("Error: Can't open FASTA file.\n");
        exit(1);
    }
    
    while (!feof(fp)) {
        ignore_path = 0;
        if (fgets(line, max_line_length, fp)) {
            if (line[0] == '>') {               
                // Read match and path number
                sscanf(line, ">match_%d_path_%d length:%d ", &m, &p, &l);

                // Check for sequential match numbers
                if ((m != current_match) && (m != (current_match+1))) {
                    printf("Error: Missing match data? Went from %d to %d.\n", current_match, m);
                    fclose(fp);
                    exit(1);
                }
                
                // Do we need to allocate more memory to store matches?
                if (m > (max_matches-1)) {              
                    int new_max_matches = max_matches + 1000;
                    Match** new_matches = realloc(matches, new_max_matches*sizeof(Match*));
                    int i;
                                        
                    if (!new_matches) {
                        printf("Error: Out of memory trying to store match %d.\n", m);
                        fclose(fp);
                        exit(1);
                    }
                    
                    for (i=max_matches; i<new_max_matches; i++) new_matches[i] = 0;
                    matches = new_matches;
                    max_matches = new_max_matches;
                }

                // Check for too many paths
                if (p >= MAX_PATHS_PER_MATCH) {
                    log_printf("Warning: Match %d path index %d too high (Max %d paths allowed). Path ignored. \n", m, p, MAX_PATHS_PER_MATCH);
                    ignore_path = 1;
                }
                                
                // Check for missing path 0.
                if ((m != current_match) && (p != 0)) {
                    printf("Error: found path %d before path 0 on match %d.\n", p, m);
                    fclose(fp);
                    exit(1);
                }
                
                // Check for sequential path numbers
                if ((m == current_match) && (p != (current_path+1))) {
                    printf("Error: found path %d before path %d on match %d.\n", p, current_path+1, m);
                    fclose(fp);
                    exit(1);
                }
                
                // Read pre, mid and post lengths
                s = strstr(line, "pre_length");
                if (s) {
                    sscanf(s, "pre_length:%d mid_length:%d post_length:%d", &pre, &mid, &post); 
                }
                
                pos = ftell(fp);                
                if (fgets(line, max_line_length, fp)) {
                    if ((strlen(line) == max_line_length) && (line[max_line_length-1] != '\n')) {
                        printf("Error: line length for match %d path %d too long. Try increasing with -y [int kb].\n", m, p);
                        exit(1);
                    }
                    
                    if (!clean_and_check_sequence(line)) {
                        printf("Error: Strange contig data for match %d path %d.\n", m, p);
                        fclose(fp);
                        exit(1);
                    }
                    
                    if (ignore_path == 0) {                    
                        if (matches[m] == 0) {
                            int i;
                            matches[m] = calloc(1, sizeof(Match));
                            current_path = -1;
                            if (!matches[m]) {
                                printf("Error: couldn't get memory for Match structure.\n");
                                fclose(fp);
                                exit(1);
                            }
                            for (i=0; i<MAX_PATHS_PER_MATCH; i++) {
                                matches[m]->paths[i] = 0;
                            }
                        }
                        
                        if (matches[m]->paths[p] == 0) {
                            matches[m]->paths[p] = malloc(sizeof(MatchPath));
                            if (!matches[m]->paths[p]) {
                                printf("Error: couldn't get memory for MatchPath structure.\n");
                                fclose(fp);
                                exit(1);
                            }
                        }
                        
                        #ifndef DEBUG_MAC
                        if (l != (pre+mid+post)) {
                            printf("  Match %d path %d: pre, mid, post (%d, %d, %d) don't add up to length (%d)\n", m, p, pre, mid, post, l);
                        }                   
                        #endif
                        
                        matches[m]->paths[p]->contig_pre = malloc(pre+1);
                        matches[m]->paths[p]->contig_mid = malloc(mid+1);
                        matches[m]->paths[p]->contig_post = malloc(post+1);
                        if ((!matches[m]->paths[p]->contig_pre) ||
                            (!matches[m]->paths[p]->contig_mid) ||
                            (!matches[m]->paths[p]->contig_post)) {
                            printf("Error: couldn't allocate memory for contig!\n");
                            exit(1);
                        }
                        strcpy_upper(matches[m]->paths[p]->contig_pre, line, pre);                  
                        strcpy_upper(matches[m]->paths[p]->contig_mid, line+pre, mid);                  
                        strcpy_upper(matches[m]->paths[p]->contig_post, line+pre+mid, post);                    
                    }
                    
                    number_of_matches = m+1;
                    current_match = m;
                    current_path = p;
                    
                    if (ignore_path == 0) {
                        matches[m]->match_number = m;
                        matches[m]->paths[p]->contig_length = l;
                        matches[m]->paths[p]->contig_file_offset = pos;
                        matches[m]->number_of_paths = p+1;
                        matches[m]->paths[p]->pre = pre;
                        matches[m]->paths[p]->mid = mid;
                        matches[m]->paths[p]->post = post;
                        matches[m]->flags = 0;
                        matches[m]->paths[p]->flags = 0;
                        //strncpy(matches[m]->paths[p]->first_bubble_kmer, line+pre-kmer_size+1, kmer_size);
                        //matches[m]->paths[p]->first_bubble_kmer[kmer_size] = 0;   
                        //matches[m]->paths[p]->first_bubble_kmer_offset = kmer_size-1;
                        get_first_and_last_kmers(line, m, p);
                    }
                } else {
                    printf("Error: couldn't find a contig for match %d path %d.\n", m, p);
                    fclose(fp);
                    exit(1);
                }
            }
        }
    }
        
    fclose(fp);
    log_printf(" ENDED Read FASTA file\n");
}
#endif

/*----------------------------------------------------------------------*
 * Function: parse_coverage_line                                        *
 * Purpose:  Parse a line from the coverage file.                       *
 * Params:   m = match number                                           *
 *           p = path number                                            *
 *           c = colour number                                          *
 *           line -> pointer to line to change.                         *
 * Returns:  0 for failure, 1 for success                               *
 *----------------------------------------------------------------------*/
int parse_coverage_line(int m, int p, int c, char* line)
{
    int n = 0;
    char* s;
    CoverageArray* cov;
    
    cov = malloc(sizeof(CoverageArray));
    if (!cov) {
        printf("Error: couldn't allocate memory to store coverage array.\n");
        return 0;
    }
    matches[m]->paths[p]->colour_coverage[c] = cov;
    
    cov->coverage = calloc(matches[m]->paths[p]->contig_length, sizeof(int));
    if (!cov->coverage) {
        printf("Error: couldn't allocate memory to store coverage.\n");
        return 0;
    }
    
    s = strtok(line, " ");
    while (s) {
        if (s) {
            if (n < matches[m]->paths[p]->contig_length) {
                cov->coverage[n++] = atoi(s);
            } else {
                printf("Error: too many coverage values for match %d path %d colour %d.\n", m, p, c);
                printf("%d %d", n, matches[m]->paths[p]->contig_length);
                return 0;
            }
        }
        s = strtok(0, " ");
    }
            
    if (n != matches[m]->paths[p]->contig_length) {
        printf("Error: not enough coverage values for match %d path %d colour %d. Expected %d, found %d.\n", m, p, c, matches[m]->paths[p]->contig_length, n);
        return 0;
    }
        
    return 1;
}

/*----------------------------------------------------------------------*
 * Function: read_coverage_file                                         *
 * Purpose:  Read the coverage file and store data.                     *
 * Params:   filename -> filename of FASTA file                         *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void read_coverage_file(char* filename)
{
    FILE *fp;
    int m,p,c;
    int current_match = -1;
    int current_path = -1;
    int ignore_path = 0;
    long pos;
    char *line;
        
    log_newline();
    log_write_timestamp(0);
    log_printf(" START Read coverage file\n");
    
    printf("\nReading coverage file...\n");

    line = malloc(max_line_length);
    if (!line) {
        printf("Error: Insufficient memory to allocate space to store line.\n");
        exit(1);
    }        
    
    fp = fopen(filename, "r");
    if (!fp) {
        printf("Error: Can't open coverage file.\n");
        exit(1);
    }

    while (!feof(fp)) {
        ignore_path = 0;
        if (fgets(line, max_line_length, fp)) {
            if (line[0] == '>') {
                sscanf(line, ">match_%d_path_%d", &m, &p);
                
                if (m > (number_of_matches-1)) {
                    printf("Error: Match number %d too high.\n", m);
                    fclose(fp);
                    exit(1);
                }

                if (matches[m] == 0) {
                    printf("Error: Match number %d doesn't exist in FASTA file.\n", m);
                    fclose(fp);
                    exit(1);
                }

                if (p >= MAX_PATHS_PER_MATCH) {
                    ignore_path = 1;
                }
                
                if ((ignore_path == 0) && (p >= matches[m]->number_of_paths)) {
                    printf("Error: Match %d path %d doesn't exist in FASTA file.\n", m, p);
                    fclose(fp);
                    exit(1);
                }
                                                
                if ((m != current_match) && (m != (current_match+1))) {
                    printf("Error: Missing match data? Went from %d to %d.\n", current_match, m);
                    fclose(fp);
                    exit(1);
                }
                
                if ((m != current_match) && (p != 0)) {
                    printf("Error: found path %d before path 0 on match %d.\n", p, m);
                    fclose(fp);
                    exit(1);
                }
                
                if ((m == current_match) && (p != (current_path+1))) {
                    printf("Error: found path %d before path %d on match %d.\n", p, current_path+1, m);
                    fclose(fp);
                    exit(1);
                }
                
                pos = ftell(fp);
                for (c=0; c<number_of_colours; c++) {
                    if (get_next_line(line, max_line_length, fp)) {                        
                        if ((strlen(line) == max_line_length-1) && (line[max_line_length-1] != '\n')) {
                            printf("Error: line length for match %d path %d too long. Try increasing with -y [int kb].\n", m, p);
                            exit(1);
                        }
                        
                        clean_line(line);
                        current_match = m;
                        current_path = p;
                        if (ignore_path == 0) {
                            matches[m]->paths[p]->coverage_file_offset = pos;
                            if (!parse_coverage_line(m, p, c, line)) {
                                fclose(fp);
                                exit(1);
                            }
                        }
                    } else {
                        printf("Error: couldn't find a contig for match %d path %d.\n", m, p);
                        fclose(fp);
                        exit(1);
                    }
                }
                
                // Add kmer to hash table
                if (ignore_path == 0) {
                    if (strlen(matches[m]->paths[p]->first_bubble_kmer) >= kmer_size) {
                        BinaryKmer b;
                        BinaryKmer tmp_kmer;
                        seq_to_binary_kmer(matches[m]->paths[p]->first_bubble_kmer, kmer_size, &b);
                        
                        // Check if it's there...
                        Element *e = hash_table_find(element_get_key(&b, kmer_size, &tmp_kmer), hash_table);
                        if (e != NULL) {
                            // If already there...
                            log_printf("Warning: Match %d path %d kmer (%s) already in table.\n", m, p, matches[m]->paths[p]->first_bubble_kmer);
                            n_kmers_already++;
                            //set_flag(&matches[m]->flags, FLAG_IGNORE);
                        } else {                                    
                            // If not found, insert it
                            Element * current_entry = hash_table_insert(element_get_key(&b, kmer_size, &tmp_kmer), hash_table); 
                            if (current_entry == NULL) {
                                printf("Error: couldn't add kmer to hash table.\n");
                                fclose(fp);
                                exit(1);
                            }
                            for (c=0; c<number_of_colours; c++) {
                                element_preallocate_quality_strings(current_entry, c, matches[m]->paths[p]->colour_coverage[c]->coverage[matches[m]->paths[p]->pre]);
                            }
                        }
                    }               
                }
            }
        }
    }
    
    fclose(fp);
    
    log_printf("%d warnings for kmers already in table.\n", n_kmers_already);
    log_write_timestamp(0);
    log_printf(" ENDED Read coverage file\n");
}

/*----------------------------------------------------------------------*
 * Function: get_bubble_length_if_equal                                 *
 * Purpose:  Get the bubble length of this match, if all paths are of   *
 *           equal length.                                              *
 * Params:   m = match number                                           *
 * Returns:  Bubble length, or 0 if not all the same length.            *
 *----------------------------------------------------------------------*/
int get_bubble_length_if_equal(int m) {
    int length_matches = matches[m]->paths[0]->mid;
    int j;
    
    for (j=1; j<matches[m]->number_of_paths; j++) {
        if (matches[m]->paths[j]->mid != matches[m]->paths[0]->mid) {
            length_matches = 0;
        }
    }   
    
    return length_matches;
}

/*----------------------------------------------------------------------*
 * Function: generate_overall_stats                                     *
 * Purpose:  Output high level stats about the files.                   *
 * Params:   None.                                                      *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void generate_overall_stats(void)
{
    int i;
    int n_kmer_length = 0;
    int n_equal_length = 0;
    int n_different_length = 0;
    int n_more_than_two_paths = 0;
    int n_more_than_two_and_kmer = 0;
    int n_more_than_two_and_equal = 0;
    int length_matches;
    
    log_newline();
    log_write_timestamp(0);
    log_printf(" START Generate overall stats\n");
    
    for (i=0; i<number_of_matches; i++) {       
        if (!(matches[i]->flags & FLAG_IGNORE)) {
            if (matches[i]->number_of_paths > 2) {
                n_more_than_two_paths++;            
            }
            
            length_matches = get_bubble_length_if_equal(i);
                    
            if (length_matches) {
                if (matches[i]->paths[0]->mid == kmer_size) {
                    n_kmer_length++;
                    if (matches[i]->number_of_paths > 2) {
                        n_more_than_two_and_kmer++;
                    }
                } else {
                    n_equal_length++;
                    if (matches[i]->number_of_paths > 2) {
                        n_more_than_two_and_equal++;
                    }
                }
            } else {
                n_different_length++;
            }
        }
    }
    
    printf("\n");
    log_and_screen_printf("Total matches: %d\n", number_of_matches);
    log_and_screen_printf("- Number with bubble path length %d: %d\n", kmer_size, n_kmer_length);
    log_and_screen_printf("- Others with equal bubble path length: %d\n", n_equal_length);
    log_and_screen_printf("- Number with differing bubble path length: %d\n", n_different_length);
    printf("\n");
    log_and_screen_printf("Number with more than two paths through bubble: %d out of %d\n", n_more_than_two_paths, number_of_matches);
    log_and_screen_printf("- Number with bubble path length %d: %d\n", kmer_size, n_more_than_two_and_kmer);
    log_and_screen_printf("- Others with equal bubble path length: %d\n", n_more_than_two_and_equal);

    log_write_timestamp(0);
    log_printf(" ENDED Generate overall stats\n");
}

/*----------------------------------------------------------------------*
 * Function: get_average_bubble_coverage                                *
 * Purpose:  Get the average coverage for route through bubble.         *
 * Params:   m = match number                                           *
 *           p = path number                                            *
 *           c = colour number                                          *
 * Returns:  Average coverage                                           *
 *----------------------------------------------------------------------*/
double get_average_bubble_coverage(int m, int p, int c, short* coverage_complete)
{
    CoverageArray* ca = matches[m]->paths[p]->colour_coverage[c];
    int i;
    int total = 0;
    int start = matches[m]->paths[p]->pre;
    int end = start + matches[m]->paths[p]->mid;
    double average = 0.0;

    *coverage_complete = 1;
    
    for (i=start; i<end; i++) {
        total = total + ca->coverage[i];
        
        if (ca->coverage[i] == 0) {
            *coverage_complete = 0;
        }
    }
    
    if (total == 0) {
        return 0;
    }
    
    //if (*coverage_complete == 1) {
        average = (double)total/(double)(end-start);
    //}
    
    return average;
}

/*----------------------------------------------------------------------*
 * Function: get_quality_scores_string                                  *
 * Purpose:  Get quality scores string for last nucleotide of given     *
 *           kmer in given colour. Produces string as sequence of       *
 *           numbers or as sequence of chars.                           * 
 * Params:   kmer_string -> kmer                                        *
 *           output_string -> storage for output string, assumed to be  *
 *                            long enough!                              *
 *           offset = offset to nucleotide (kmer_length-1 usually)      *
 *           c = colour to get quality scores for                       *
 *           max_number = maximum number of scores to fetch             *
 *           as_numbers = 1 to return numbers, or 0 to return chars     *
 * Returns:  Pointer to output_string                                   *
 *----------------------------------------------------------------------*/
char* get_quality_scores_string(char* kmer_string, char* output_string, int offset, int c, int max_number, int as_numbers)
{
    BinaryKmer b;
    BinaryKmer tmp_kmer;
    int i, n;
    char tmp_string[16];
    
    strcpy(output_string, "");
    
    seq_to_binary_kmer(kmer_string, kmer_size, &b); 
    Element *e = hash_table_find(element_get_key(&b, kmer_size, &tmp_kmer), hash_table);
    if (e == NULL) {
        printf("Error: can't find kmer %s\n", kmer_string);
        exit(1);
    }
    
    QualityStringArray *qa = &e->quality_string_arrays[c];
        
    n = qa->number_of_strings;
    if (n > max_number) {
        n = max_number;
    }
    
    for (i=0; i<n; i++) {
        char* qs = qa->quality_strings[i].quality;
        if ((as_numbers) && (i > 0)) {
            strcat(output_string, " ");
        }
        
        if (as_numbers) {
            sprintf(tmp_string, "%d", qs[offset]-quality_offset);
        } else {
            sprintf(tmp_string, "%c", qs[offset]);
        }
        strcat(output_string, tmp_string);
    }
    
    if (!as_numbers) {
        for (i=n; i<max_number; i++) {
            strcat(output_string, " ");
        }
    }
    
    return output_string;
}

/*----------------------------------------------------------------------*
 * Function: get_all_colour_quality_scores_strings                      *
 * Purpose:  Get quality scores string for last nucleotide of given     *
 *           kmer in all colours. Produces string as sequence of        *
 *           numbers or as sequence of chars. Each colour begins with   *
 *           c0(...) etc.                                               *
 * Params:   kmer_string -> kmer                                        *
 *           output_string -> storage for output string, assumed to be  *
 *                            long enough!                              *
 *           offset = offset to nucleotide (kmer_length-1 usually)      *
 *           max_number = maximum number of scores to fetch             *
 *           as_numbers = 1 to return numbers, or 0 to return chars     *
 * Returns:  Pointer to output_string                                   *
 *----------------------------------------------------------------------*/
char* get_all_colour_quality_scores_strings(char* kmer_string, char* output_string, int offset, int max_number, int as_numbers)
{
    int c;
    char tmp_string[1024];
    char c_string[64];
    
    strcpy(output_string, "");
    
    for (c=0; c<number_of_colours; c++) {
        get_quality_scores_string(kmer_string, tmp_string, offset, c, max_number, as_numbers);
        if (c > 0) {
            sprintf(c_string, " ");
        }
        sprintf(c_string, "c%d(", c);
        strcat(output_string, c_string);
        strcat(output_string, tmp_string);
        strcat(output_string, ")");
    }
    
    return output_string;
}

/*----------------------------------------------------------------------*
 * Function: p_to_q, p_to_q_char, q_to_p, char_q_to_p                   *
 * Purpose:  Conversion between p-value and Q for quality scores.       *
 * Params:   N/A                                                        *
 * Returns:  N/A                                                        *
 *----------------------------------------------------------------------*/
short p_to_q(double p)
{
    return -10 * log10(p);
}

char p_to_q_char(double p)
{
    return p_to_q(p)+quality_offset;
}

double q_to_p(double q)
{
    return pow(10, q/-10.0);
}

double char_q_to_p(char c)
{
    double q = c - quality_offset;
    return q_to_p(q);
}

/*----------------------------------------------------------------------*
 * Function: check_all_paths_at_least_kmer_size                         *
 * Purpose:  Check that all paths for a match are at least kmer size in *
 *           length.                                                    *
 * Params:   m = match number to calculate cpnp for.                    *
 * Returns:  1 if all are at least kmer size, 0 otherwise               *
 *----------------------------------------------------------------------*/
int check_all_paths_at_least_kmer_size(int m)
{
    int p;
    int rc = 1;
    
    for (p=0; p<matches[m]->number_of_paths; p++) {
        if (matches[m]->paths[p]->mid < kmer_size) {
            rc = 0;
            break;
        }
    }
    
    return rc;
}

/*----------------------------------------------------------------------*
 * Function: calculate_probabilities                                    *
 * Purpose:  Calculate ranking probabilities                            *
 * Params:   m = match number to calculate for                          *
 * Returns:  cpnp value                                                 *
 *----------------------------------------------------------------------*/
void calculate_probabilities(int m)
{
    BinaryKmer b;
    BinaryKmer tmp_kmer;
    int c, i, p;
    int first_node;
    short total_coverage;
    long double combined_p;
    long double cpnp = 1.0;
    long double pn;
    long double num, den;
    char* kmer_string;

    // Check all paths are kmer size...
    if (!check_all_paths_at_least_kmer_size(m)) {
        return;
    }
    
    matches[m]->statistics.q_total = 0;
    
    // Build up an array of p_values
    for (p=0; p<matches[m]->number_of_paths; p++)
    {
        // Get the index to the first node of the bubble
        first_node = matches[m]->paths[p]->pre;

        // Get the quality string for the first kmer in this path
        kmer_string = matches[m]->paths[p]->first_bubble_kmer;

        //log_printf(" [%s]", matches[m]->paths[p]->first_bubble_kmer);
        
        // Find first kmer in hash table
        seq_to_binary_kmer(kmer_string, kmer_size, &b); 
        Element *e = hash_table_find(element_get_key(&b, kmer_size, &tmp_kmer), hash_table);
        if (e == NULL) {
            printf("Error: can't find kmer %s\n", kmer_string);
            exit(1);
        }
        
        // Find total coverage across all colours for first node of this path
        total_coverage = 0;
        for (c=0; c<number_of_colours; c++) {
            total_coverage += matches[m]->paths[p]->colour_coverage[c]->coverage[first_node];
        }
        
        // Now go through colours, calculating combined p_value
        for (c=0; c<number_of_colours; c++) {
            QualityStringArray *qa = &e->quality_string_arrays[c];
            int q_total = 0;
            int number_of_quality_scores = qa->number_of_strings;
            
            // Warnings if number of quality scores found in file different to coverage
            if ((matches[m]->paths[p]->colour_coverage[c]->coverage[first_node] > 0) &&
                (number_of_quality_scores != matches[m]->paths[p]->colour_coverage[c]->coverage[first_node])) {
                log_printf("Warning: Number of quality scores read (%d) doesn't match expected coverage (%d) for kmer (%s) on match %d, path %d, colour %d\n",
                           number_of_quality_scores,
                           matches[m]->paths[p]->colour_coverage[c]->coverage[first_node],
                           kmer_string,
                           m,
                           p,
                           c);
                if (number_of_quality_scores > matches[m]->paths[p]->colour_coverage[c]->coverage[first_node]) {
                    no_match_high++;
                    set_flag(&matches[m]->flags, FLAG_NO_MATCH_HIGH);
                    set_flag(&matches[m]->paths[p]->flags, FLAG_NO_MATCH_HIGH);
                } else {
                    no_match_low++;
                    set_flag(&matches[m]->flags, FLAG_NO_MATCH_LOW);
                    set_flag(&matches[m]->paths[p]->flags, FLAG_NO_MATCH_LOW);
                }
            }
            
            // CPNP calculation

            // Sum up quality scores - these become the exponent
            for (i=0; i<qa->number_of_strings; i++) {
                char* qs = qa->quality_strings[i].quality;
                q_total += qs[matches[m]->paths[p]->first_bubble_kmer_offset]-quality_offset;
            }
                        
            // Workout combined (combination of quality scores) p-value for this path and colour
            combined_p = powl(10.0, (q_total/-10));
            combined_p *= matches[m]->paths[p]->colour_coverage[c]->coverage[first_node]; 
            combined_p /= total_coverage;
                        
            // Now, we want to build a probability that this is correct, so the probablity that
            // this path and colour is correct is 1-combined_p. The CPNP is obtained by
            // multiplying all these values together.
            cpnp *= (1-combined_p);
            
            // Alternative to CPNP - total quality score
            
            // Update Q total
            matches[m]->statistics.q_total += q_total;

            // Another alternative - sum up the probabilities
            pn = 0;
            for (i=0; i<qa->number_of_strings; i++) {
                char* qs = qa->quality_strings[i].quality;
                int qi = qs[matches[m]->paths[p]->first_bubble_kmer_offset]-quality_offset;
                pn = pn + (1 - powl(10.0, (qi/-10)));
            }
            if (matches[m]->statistics.coverage_complete[p][c]) {
                matches[m]->statistics.p_value[p][c] = pn;
            } else {
                matches[m]->statistics.p_value[p][c] = 0;
            }
        }
    }
    
    // Now work out the p ratios
    for (c=0; c<number_of_colours; c++) {
        if (matches[m]->number_of_paths == 2) {
            if (matches[m]->statistics.p_value[0][c] > matches[m]->statistics.p_value[1][c]) {
                num = matches[m]->statistics.p_value[0][c];
                den = matches[m]->statistics.p_value[1][c];
            } else {
                num = matches[m]->statistics.p_value[1][c];
                den = matches[m]->statistics.p_value[0][c];
            }
            
            if (den > 0) {
                matches[m]->statistics.p_ratio[c] = num/den;
            } else {
                matches[m]->statistics.p_ratio[c] = num/0.1;
            }               
        } else {
            matches[m]->statistics.p_ratio[0] = 0;
            matches[m]->statistics.p_ratio[1] = 0;
        }
    }
    
    matches[m]->statistics.cpnp = cpnp;
}

/*----------------------------------------------------------------------*
 * Function: calculate_coverage_difference                              *
 * Purpose:  Calculate difference between coverage % and expected %.    *
 * Params:   m = match number to calculate for                          *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void calculate_coverage_difference(int m)
{
    int c = 0;
    
    for (c=0; c<number_of_colours; c++) {
        if ((matches[m]->statistics.coverage_pc[0][c] == 0) && (matches[m]->statistics.coverage_pc[1][c] == 0)) {
            matches[m]->statistics.coverage_difference[c] = expected_coverage_pc[0][c] > expected_coverage_pc[1][c] ? expected_coverage_pc[0][c] : expected_coverage_pc[1][c];
        } else if (fabs(expected_coverage_pc[0][c] - matches[m]->statistics.coverage_pc[0][c]) <
            fabs(expected_coverage_pc[1][c] - matches[m]->statistics.coverage_pc[0][c])) {
            matches[m]->statistics.coverage_difference[c] = fabs(expected_coverage_pc[0][c] - matches[m]->statistics.coverage_pc[0][c]);
        } else {
            matches[m]->statistics.coverage_difference[c] = fabs(expected_coverage_pc[1][c] - matches[m]->statistics.coverage_pc[0][c]);
        }
    }
}

/* Sorting function for find_type */
int int_cmp(const void *a, const void *b)
{
    const int *ia = (const int *)a;
    const int *ib = (const int *)b;
    return *ib  - *ia; 
}

/*----------------------------------------------------------------------*
 * Function: find_type                                                  *
 * Purpose:  Work out a type string.                                    *
 * Params:   m = match number to calculate for                          *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void find_type(int m)
{
    int c, p;
    int c_on_p[MAX_PATHS_PER_MATCH];
    char type_string[64];
    char temp[16];
    int n_paths;
    
    // Clear array
    for (p=0; p<MAX_PATHS_PER_MATCH; p++) {
        c_on_p[p] = 0;
    }
    
    // Work out how many colours on each path
    for (p=0; p<matches[m]->number_of_paths; p++) {
        for (c=0; c<number_of_colours; c++) {
            if (matches[m]->statistics.coverage_complete[p][c]) {
                c_on_p[p]++;
            }
        }
    }
    
    // Sort array - lowest first
    qsort(c_on_p, MAX_PATHS_PER_MATCH, sizeof(int), int_cmp);
    
    type_string[0] = 0;
    n_paths = 0;
    for (p=0; p<MAX_PATHS_PER_MATCH; p++) {
        if (c_on_p[p] > 0) {
            n_paths++;
            
            if (p > 0) {
                sprintf(temp, ",%d", c_on_p[p]);
            } else {
                sprintf(temp, "%d", c_on_p[p]);
            }
                        
            strcat(type_string, temp);
        }
    }
    
    if (n_paths == 0) {
        strcpy(matches[m]->statistics.type, "z");
    } else if (n_paths == 1) {
        sprintf(matches[m]->statistics.type, "y %s", type_string);
    } else if (matches[m]->statistics.flags & FLAG_PATHS_SAME_LENGTH) {
        sprintf(matches[m]->statistics.type, "S %s", type_string);
    } else {
        sprintf(matches[m]->statistics.type, "i %s", type_string);
    }
}

/*----------------------------------------------------------------------*
 * Function: Convert type_for_csv                                       *
 * Purpose:  Converts type to one suitable for CSV file, by replacing   *
 *           commas and spaces with underscore character.               *
 * Params:   type -> type string                                        *
 *           converted -> buffer for converted string                   *
 * Returns:  Pointer to converted string.                               *
 *----------------------------------------------------------------------*/
char* convert_type_for_csv(char* type, char* converted)
{
    int i;
    
    for (i=0; i<strlen(type); i++) {
        if ((type[i] == ',') || (type[i] == ' ')) {
            converted[i] = '_';
        } else {
            converted[i] = type[i];
        }
    }
    
    converted[i] = 0;
    
    return converted;
}

/*----------------------------------------------------------------------*
 * Function: score_matches                                              *
 * Purpose:  Score matches (calculate stats) prior to ranking           *
 * Params:   filename -> filename of log file                           *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void score_matches(void)
{
    int m, p, c;
    double total;
    double lowest;
    FILE *fp = 0;
    char qs_string[1024];
    
    log_newline();
    log_write_timestamp(0);
    log_printf(" START Score matches\n");

    if (rank_log_filename) {
        fp = fopen(rank_log_filename, "w");
        if (!fp) {
            printf("Error: Can't open scoring log.\n");
            exit(1);
        }
    }
    
    printf("\nScoring matches...\n");
    
    for (m=0; m<number_of_matches; m++) {        
        if (((m % 100000) == 0) || (m == number_of_matches-1)) {
            log_printf("Reached match %d\n", m);
        }
        
        //log_printf("Match %d...", m);
        matches[m]->statistics.combined_coverage = 0;
        if (!(matches[m]->flags & FLAG_IGNORE)) {
            // Clear coverage matrix
            for (p=0; p<MAX_PATHS_PER_MATCH; p++) {
                for (c=0; c<number_of_colours; c++) {
                    matches[m]->statistics.coverage_av[p][c] = 0;
                    matches[m]->statistics.coverage_pc[p][c] = 0;
                }
            }
            
            //log_printf(" 1");

            // Fill coverage matrix
            for (c=0; c<number_of_colours; c++) {
                for (p=0; p<matches[m]->number_of_paths; p++) {
                    matches[m]->statistics.coverage_av[p][c] = get_average_bubble_coverage(m, p, c, &(matches[m]->statistics.coverage_complete[p][c]));
                    matches[m]->statistics.combined_coverage += matches[m]->statistics.coverage_av[p][c];
                }
            }

            //log_printf(" 2");
            
            // Check if bubble paths are the same length
            set_flag(&matches[m]->statistics.flags, FLAG_PATHS_SAME_LENGTH);
            matches[m]->longest_contig = 0;
            for (p=1; p<matches[m]->number_of_paths; p++) {
                if (matches[m]->paths[p]->mid != matches[m]->paths[0]->mid) {
                    unset_flag(&matches[m]->statistics.flags, FLAG_PATHS_SAME_LENGTH);
                }
                if (matches[m]->paths[p]->contig_length > matches[m]->longest_contig) {
                    matches[m]->longest_contig = matches[m]->paths[p]->contig_length;
                }
            }

            //log_printf(" 3");
            
            // Make it a percentage
            lowest = 100.0;
            for (c=0; c<number_of_colours; c++) {
                total = 0.0;
                for (p=0; p<matches[m]->number_of_paths; p++) {
                    if (matches[m]->statistics.coverage_complete[p][c]) {
                        total += matches[m]->statistics.coverage_av[p][c];
                    }
                }
                if (total > 0) {
                    for (p=0; p<matches[m]->number_of_paths; p++) {
                        if ((matches[m]->statistics.coverage_av[p][c] == 0) ||
                            (matches[m]->statistics.coverage_complete[p][c] == 0)) {
                            matches[m]->statistics.coverage_pc[p][c] = 0.0;
                        } else {
                            matches[m]->statistics.coverage_pc[p][c] = (100.0*matches[m]->statistics.coverage_av[p][c]) / total;
                        }
                        
                        if (matches[m]->statistics.coverage_pc[p][c] < lowest) {
                            lowest = matches[m]->statistics.coverage_pc[p][c];
                        }
                    }
                } else {
                    lowest = 0;
                }
            }

            //log_printf(" 4");
            
            // Find type
            find_type(m);

            //log_printf(" 5");
            
            // Check if coverage percentage within tolerance
            calculate_coverage_difference(m);
            set_flag(&matches[m]->statistics.flags, FLAG_COVERAGE_WITHIN_BOUNDS);
            for (c=0; c<number_of_colours; c++) {
                if (matches[m]->statistics.coverage_difference[c] > coverage_pc_tolerance[c]) {
                    unset_flag(&matches[m]->statistics.flags, FLAG_COVERAGE_WITHIN_BOUNDS);
                    break;
                }
            }

            //log_printf(" 6");
            
            // Calculate probabilities for ranking
            calculate_probabilities(m);
                        
            //log_printf(" 7");
            
            // Output to log file
            if (rank_log_filename) {
                fprintf(fp, "Match %d\n", m);
                fprintf(fp, "  Type: %s\n", matches[m]->statistics.type);

    #ifdef USE_CPNP_RANKING
                fprintf(fp, "  CPNP: %Le\n", matches[m]->statistics.cpnp);
    #endif
                
                for (p=0; p<matches[m]->number_of_paths; p++) {
                    if (strlen(matches[m]->paths[p]->first_bubble_kmer) != kmer_size) {
                        printf("Error: match %d path %d first kmer (%s) not equal to kmer_size\n", m, p, matches[m]->paths[p]->first_bubble_kmer);
                        exit(1);
                    }
                    
                    fprintf(fp, "  Path %d\n", p);
                    fprintf(fp, "                    First kmer: %s\n", matches[m]->paths[p]->first_bubble_kmer);
                    fflush(fp);
                    for (c=0; c<number_of_colours; c++) {
                        fprintf(fp, "          Coverage in colour %d: %d\n", c, matches[m]->paths[p]->colour_coverage[c]->coverage[matches[m]->paths[p]->pre]);
                        fprintf(fp, "    Coverage mean for colour %d: %.2f\n", c, matches[m]->statistics.coverage_av[p][c]);
                        fprintf(fp, "       Coverage %% for colour %d: %.2f\n", c, matches[m]->statistics.coverage_pc[p][c]);
                        fprintf(fp, "       Quality scores colour %d: %s\n", c, get_quality_scores_string(matches[m]->paths[p]->first_bubble_kmer, qs_string, matches[m]->paths[p]->first_bubble_kmer_offset, c, 256, TRUE));
                        fprintf(fp, "                 as characters: %s\n", get_quality_scores_string(matches[m]->paths[p]->first_bubble_kmer, qs_string, matches[m]->paths[p]->first_bubble_kmer_offset, c, 256, FALSE));
                        fprintf(fp, "    Coverage diff for colour %d: %.2f\n", c, matches[m]->statistics.coverage_difference[c]);
                        fflush(fp);
                    }
                }           
            }
        }
        
        //log_printf(" 8\n");

        if (rank_log_filename) {
            fprintf(fp, "\n");
        }
    }
    
    printf("Kmers with no. quality scores > coverage: %d\n", no_match_high);
    printf("Kmers with no. quality scores < coverage: %d\n", no_match_low);

    if ((rank_log_filename) && (fp)) {
        fclose(fp);
    }

    log_write_timestamp(0);
    log_printf(" ENDED Score matches\n");
}

/*----------------------------------------------------------------------*
 * Function: match_alpha_cmp                                            *
 * Purpose:  qsort matching function for Match* array                   *
 * Params:   a -> first item to compare                                 *
 *           b -> second item to compare                                *
 * Returns:  -ve value if a should be nearer the top of the list than b *
 *           +ve value if a should be lower down the list than b        *
 *           0 if a and b are equal                                     *
 *----------------------------------------------------------------------*/
int match_alpha_cmp(const void *a, const void *b)
{
    const Match **pma = (const Match **)a;
    const Match **pmb = (const Match **)b;
    Match *ma = (Match*)*pma;
    Match *mb = (Match*)*pmb;
    
    return strcmp(ma->paths[0]->first_contig_kmer, mb->paths[0]->first_contig_kmer);
}    

/*----------------------------------------------------------------------*
 * Function: match_rank                                                 *
 * Purpose:  qsort matching function for Match* array                   *
 * Params:   a -> first item to compare                                 *
 *           b -> second item to compare                                *
 * Returns:  -ve value if a should be nearer the top of the list than b *
 *           +ve value if a should be lower down the list than b        *
 *           0 if a and b are equal                                     *
 *----------------------------------------------------------------------*/
int match_rank(const void *a, const void *b)
{
    const Match **pma = (const Match **)a;
    const Match **pmb = (const Match **)b;
    Match *ma = (Match*)*pma;
    Match *mb = (Match*)*pmb;   
    int i;
        
    // We favour one where the ignore flag isn't set!
    if ((!(ma->flags & FLAG_IGNORE)) && (mb->flags & FLAG_IGNORE)) {
        return -1;
    } else if ((ma->flags & FLAG_IGNORE) && (!(mb->flags & FLAG_IGNORE))) {
        return 1;
    }
    
    // Rank on type
    i = strcmp(ma->statistics.type, mb->statistics.type);
    if (i != 0)
        return i;
    
    // Put repeats at the bottom
    if ((!(ma->flags & FLAG_REPEAT)) && (mb->flags & FLAG_REPEAT)) {
        return -1;
    } else if ((ma->flags & FLAG_REPEAT) && (!(mb->flags & FLAG_REPEAT))) {
        return 1;
    }
    
    // Must be at least minimum_contig_size
    if ((ma->longest_contig < minimum_contig_size) && (mb->longest_contig >= minimum_contig_size)) {
        return 1;
    } else if ((ma->longest_contig >= minimum_contig_size) && (mb->longest_contig < minimum_contig_size)) {
        return -1;
    }
        
    // We favour one where bubble paths are the same length
    if ((!(ma->statistics.flags & FLAG_PATHS_SAME_LENGTH)) && (mb->statistics.flags & FLAG_PATHS_SAME_LENGTH)) {
        return 1;
    } else if ((ma->statistics.flags & FLAG_PATHS_SAME_LENGTH) && (!(mb->statistics.flags & FLAG_PATHS_SAME_LENGTH))) {
        return -1;
    }
    
    // We favour SNPs (where bubble path = kmer size)
    if ((ma->statistics.flags & FLAG_PATHS_SAME_LENGTH) && (mb->statistics.flags & FLAG_PATHS_SAME_LENGTH)) {
        if ((ma->paths[0]->mid != kmer_size) && (mb->paths[0]->mid == kmer_size)) {
            return 1;
        } else if ((ma->paths[0]->mid == kmer_size) && (mb->paths[0]->mid != kmer_size)) {
            return -1;
        }
    }
    
    // We favour coverage within the tolerance bounds   
    if ((!(ma->statistics.flags & FLAG_COVERAGE_WITHIN_BOUNDS)) && (mb->statistics.flags & FLAG_COVERAGE_WITHIN_BOUNDS)) {
        return 1;
    } else if ((ma->statistics.flags & FLAG_COVERAGE_WITHIN_BOUNDS) && (!(mb->statistics.flags & FLAG_COVERAGE_WITHIN_BOUNDS))) {
        return -1;
    }

#ifdef USE_CPNP_RANKING
    // On CPNP
    if (ma->statistics.cpnp > mb->statistics.cpnp) {
        return -1;
    } else if (ma->statistics.cpnp < mb->statistics.cpnp) {
        return 1;
    }
#endif
    
    // On total Q
    if (mb->statistics.q_total != ma->statistics.q_total)
        return mb->statistics.q_total - ma->statistics.q_total;
    
    // On combined coverage
    return mb->statistics.combined_coverage - ma->statistics.combined_coverage;
}

/*----------------------------------------------------------------------*
 * Function: rank_matches                                               *
 * Purpose:  Sort matches.                                              *
 * Params:   None.                                                      *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void rank_matches(void)
{
    log_newline();
    log_write_timestamp(0);
    log_printf(" START Rank matches\n");

    printf("\nSorting matches...\n");
    
    qsort(matches, number_of_matches, sizeof(Match*), match_rank);

    log_write_timestamp(0);
    log_printf(" ENDED Rank matches\n");
}

/*----------------------------------------------------------------------*
 * Function: add_to_table_with_span                                     *
 * Purpose:  Add item to table and make it span multiple columns.       *
 * Params:   row -> the row we're building up                           *
 *           item -> item to add (a string)                             *
 *           column -> column counter                                   *
 *           span = number of columns to span                           *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void add_to_table_with_span(char *row, char* item, int* column, int span)
{
    int width = 0;
    int spare = 0;
    int left = 0;
    int right = 0;
    int i;
    
    for (i=0; i<span; i++) {
        width += column_widths[(*column)+i];
        
        if (i > 0) {
            width += space_between_columns;
        }
    }

    spare = width - strlen(item);
    left = spare / 2;
    right = spare - left;    
    
    if (strlen(row) > 0) {
        for (i=0; i<space_between_columns; i++) {
            strcat(row, " ");
        }
    }
    
    if (left < 0) {
        left = 0;
    }
    
    if (right < 0) {
        right = 0;
    }
    
    for (i=0; i<left; i++) {
        strcat(row, " ");
    }
    
    strcat(row, item);
    
    for (i=0; i<right; i++) {
        strcat(row, " ");
    }
    
    *column = (*column)+span;
}

void add_to_table(char* row, char* item, int* column)
{
    add_to_table_with_span(row, item, column, 1);
}

void add_int_to_table(char* row, char* format, int i, int* column)
{
    char temp[256];
    
    sprintf(temp, format, i);
    add_to_table(row, temp, column);
}

void add_double_to_table(char* row, char* format, double d, int* column)
{
    char temp[256];
    
    sprintf(temp, format, d);
    add_to_table(row, temp, column);
}

void add_long_double_to_table(char* row, char* format, long double d, int* column)
{
    char temp[256];
    
    sprintf(temp, format, d);
    add_to_table(row, temp, column);
}


/*----------------------------------------------------------------------*
 * Function: output_header_line                                         *
 * Purpose:  Print the rank table header line to screen and file.       *
 * Params:   fp -> handle of file to print to                           *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void output_header_line(FILE *fp)
{
    char line1[1024];
    char line2[1024];
    int i;
    int totalwidth = 0;
    int column1 = 0;
    int column2 = 0;
            
    fputs("\n", fp); 
    
    column1 = 0;
    column2 = 0;
    line1[0] = 0;
    line2[0] = 0;
    
    add_to_table(line1, "", &column1);
    add_to_table(line2, "Rank", &column2);    // 0
    
    add_to_table(line1, "", &column1);
    add_to_table(line2, "Match", &column2);   // 1

    add_to_table(line1, "Num", &column1);
    add_to_table(line2, "Pth", &column2);     // 2

    add_to_table(line1, "", &column1);
    add_to_table(line2, "Type", &column2);    // 3

    add_to_table(line1, "Lngst", &column1);
    add_to_table(line2, "Cntig", &column2);   // 4

    add_to_table_with_span(line1, "Path Len", &column1, 3);
    add_to_table(line2, "P0", &column2);      // 5
    add_to_table(line2, "P1", &column2);      // 6
    add_to_table(line2, "P2", &column2);      // 7

    add_to_table(line1, "", &column1);
    add_to_table(line2, "Flags", &column2);   // 8

    add_to_table_with_span(line1, "c0 Coverage", &column1, 3);
    add_to_table(line2, "P0", &column2);      // 9
    add_to_table(line2, "P1", &column2);      // 10
    add_to_table(line2, "P2", &column2);      // 11

    add_to_table_with_span(line1, "c1 Coverage", &column1, 3);
    add_to_table(line2, "P0", &column2);      // 12
    add_to_table(line2, "P1", &column2);      // 13
    add_to_table(line2, "P2", &column2);      // 14

    add_to_table_with_span(line1, "c0 Coverage %", &column1, 3);
    add_to_table(line2, "P0", &column2);      // 15
    add_to_table(line2, "P1", &column2);      // 16
    add_to_table(line2, "P2", &column2);      // 17

    add_to_table_with_span(line1, "c1 Coverage %", &column1, 3);
    add_to_table(line2, "P0", &column2);      // 18
    add_to_table(line2, "P1", &column2);      // 19
    add_to_table(line2, "P2", &column2);      // 20
 
    add_to_table_with_span(line1, "Difference", &column1, 2);
    add_to_table(line2, "c0", &column2);      // 21
    add_to_table(line2, "c1", &column2);      // 22
    
    add_to_table(line1, "", &column1);
    add_to_table(line2, "QTotal", &column2);  // 23

#ifdef EXTRA_STATS
    add_to_table(line1, "", &column1);
    add_to_table(line2, "CPNP", &column2);    // 24

    add_to_table_with_span(line1, "p ratio", &column1, 2);
    add_to_table(line2, "c0", &column2);      // 25
    add_to_table(line2, "c1", &column2);      // 26
    
    add_to_table_with_span(line1, "c0 p-value", &column1, 2);
    add_to_table(line2, "P0", &column2);      // 27
    add_to_table(line2, "P1", &column2);      // 28
    
    add_to_table_with_span(line1, "c1 p-value", &column1, 2);
    add_to_table(line2, "P0", &column2);      // 29
    add_to_table(line2, "P1", &column2);      // 30
#endif
    
    fputs(line1, fp);
    fputs("\n", fp);
    fputs(line2, fp);
    fputs("\n", fp);
    
    for (i=0; i<number_of_columns; i++) {
        if (i > 0) {
            totalwidth += space_between_columns;
        }
        totalwidth += column_widths[i];
    }
    
    for (i=0; i<totalwidth; i++) {
        fputs("-", fp);
    }
    fputs("\n", fp);
}

/*----------------------------------------------------------------------*
 * Function: makes_flag_string                                          *
 * Purpose:  Make a string showing which flags of a match are on.       *
 * Params:   tmp_string -> Place to store string.                       *
 *           m = match number.                                          *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void make_flags_string(char* tmp_string, int m)
{
    tmp_string[0] = 0;
    if (matches[m]->paths[0]->flags & FLAG_NO_MATCH_LOW) {
        strcat(tmp_string, "L");
    } else if (matches[m]->paths[0]->flags & FLAG_NO_MATCH_HIGH) {
        strcat(tmp_string, "H");
    } else {
        strcat(tmp_string, "-");
    }
    
    if (matches[m]->paths[1]->flags & FLAG_NO_MATCH_LOW) {
        strcat(tmp_string, "L");
    } else if (matches[m]->paths[1]->flags & FLAG_NO_MATCH_HIGH) {
        strcat(tmp_string, "H");
    } else {
        strcat(tmp_string, "-");
    }
    
    if (matches[m]->statistics.flags & FLAG_PATHS_SAME_LENGTH) {
        strcat(tmp_string, "S");
    } else {
        strcat(tmp_string, "-");
    }
    
    if (matches[m]->statistics.flags & FLAG_COVERAGE_WITHIN_BOUNDS) {
        strcat(tmp_string, "W");
    } else {
        strcat(tmp_string, "-");
    }
    
    if (matches[m]->flags & FLAG_IGNORE) {
        strcat(tmp_string, "I");
    } else {
        strcat(tmp_string, "-");
    }                               
    
    if (matches[m]->flags & FLAG_REPEAT) {
        strcat(tmp_string, "R");
    } else {
        strcat(tmp_string, "-");
    }    
}

/*----------------------------------------------------------------------*
 * Function: output_rank_table                                          *
 * Purpose:  Output ranked list                                         *
 * Params:   None.                                                      *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void output_rank_table(void)
{
    int n_output_tab = 0;
    int n_output_csv = 0;
    int m, p, pto;
    FILE* table_fp = 0;
    FILE* csv_fp = 0;
    char table_line[1024];
    char csv_line[1024];
    char tmp_string[64];
    char previous_type[64];
    int rank_count = 1;
    int max_paths_to_show = 3;
    int column = 0;
    
    log_newline();
    log_write_timestamp(0);
    log_printf(" START Output rank table\n");

    previous_type[0] = 0;
    
    printf("\nWriting rank file...\n");
    
    if (rank_table_filename) {
        table_fp = fopen(rank_table_filename, "w");
        if (!table_fp) {
            printf("Error: Can't open table file %s\n", rank_table_filename);
        }
    }

    if (rank_csv_filename) {
        csv_fp = fopen(rank_csv_filename, "w");
        if (!csv_fp) {
            printf("Error: Can't open table file %s\n", rank_csv_filename);
        }
    }        
    
    if (csv_fp) {
        fputs("Rank,Match,NumPth,Type,LngstCntig,LenP0,LenP1,CovC0P0,CovC0P1,CovC1P0,CovC1P1,PcC0P0,PcC0P1,PcC1P0,PcC1P1,C0Dif,C1Dif,QTotal\n", csv_fp);
    }
    
    strcpy(previous_type, matches[0]->statistics.type);

    if (table_fp) {
        output_header_line(table_fp);
    }
    
    for (m=0; m<number_of_matches; m++) {
        column = 0;
        table_line[0] = 0;
        
        if (strcmp(previous_type, matches[m]->statistics.type) != 0) {
            strcpy(previous_type, matches[m]->statistics.type);
            if (table_fp) {
                output_header_line(table_fp);
            }
            rank_count = 1;
        }
        
        if (matches[m]->number_of_paths > 0) {
            add_int_to_table(table_line, "%d", rank_count, &column);
            add_int_to_table(table_line, "%d", matches[m]->match_number, &column);          
            add_int_to_table(table_line, "%d", matches[m]->number_of_paths, &column);

            strcpy(tmp_string, matches[m]->statistics.type);
            make_upper_case(tmp_string);
            add_to_table(table_line, tmp_string, &column);

            add_int_to_table(table_line, "%d", matches[m]->longest_contig, &column);
            
            pto = matches[m]->number_of_paths <= max_paths_to_show ? matches[m]->number_of_paths:max_paths_to_show;
            for (p=0; p<pto; p++) {             
                add_int_to_table(table_line, "%d", matches[m]->paths[p]->mid, &column);
            }           
            for (p=matches[m]->number_of_paths; p<max_paths_to_show; p++) {
                add_to_table(table_line, "", &column);
            }
                                                                
            make_flags_string(tmp_string, m);
            add_to_table(table_line, tmp_string, &column);
            
            for (p=0; p<3; p++) {
                if ((p < matches[m]->number_of_paths) && (matches[m]->statistics.coverage_av[p][0] > 0.0)) {
                    char temp[64];
                    
                    if (matches[m]->statistics.coverage_complete[p][0]) {
                        sprintf(temp, "%.2f", matches[m]->statistics.coverage_av[p][0]);
                    } else {
                        sprintf(temp, "(%.2f)", matches[m]->statistics.coverage_av[p][0]);
                    }                    
                    add_to_table(table_line, temp, &column);
                    //add_double_to_table(table_line, "%.2f", matches[m]->statistics.coverage_av[p][0], &column);
                } else {
                    add_to_table(table_line, "", &column);
                }
            }
            
            for (p=0; p<3; p++) {
                if ((p < matches[m]->number_of_paths) && (matches[m]->statistics.coverage_av[p][1] > 0.0)) {
                    char temp[64];
                    
                    if (matches[m]->statistics.coverage_complete[p][1]) {
                        sprintf(temp, "%.2f", matches[m]->statistics.coverage_av[p][1]);
                    } else {
                        sprintf(temp, "(%.2f)", matches[m]->statistics.coverage_av[p][1]);
                    }                    
                    add_to_table(table_line, temp, &column);
                    //add_double_to_table(table_line, "%.2f", matches[m]->statistics.coverage_av[p][1], &column);
                } else {
                    add_to_table(table_line, "", &column);
                }
            }
                        
            for (p=0; p<3; p++) {
                if (p < matches[m]->number_of_paths) {
                    add_double_to_table(table_line, "%.2f", matches[m]->statistics.coverage_pc[p][0], &column);
                } else {
                    add_to_table(table_line, "", &column);
                }
            }
            
            for (p=0; p<3; p++) {
                if (p < matches[m]->number_of_paths) {
                    add_double_to_table(table_line, "%.2f", matches[m]->statistics.coverage_pc[p][1], &column);
                } else {
                    add_to_table(table_line, "", &column);
                }
            }
            
            add_double_to_table(table_line, "%.2f", matches[m]->statistics.coverage_difference[0], &column);
            add_double_to_table(table_line, "%.2f", matches[m]->statistics.coverage_difference[1], &column);
            
            add_int_to_table(table_line, "%d", matches[m]->statistics.q_total, &column);
            
#ifdef EXTRA_STATS
            add_long_double_to_table(table_line, "%Le", matches[m]->statistics.cpnp, &column);
            add_long_double_to_table(table_line, "%.4Lf", matches[m]->statistics.p_ratio[0], &column);
            add_long_double_to_table(table_line, "%.4Lf", matches[m]->statistics.p_ratio[1], &column);
            add_long_double_to_table(table_line, "%.4Lf", matches[m]->statistics.p_value[0][0], &column);
            add_long_double_to_table(table_line, "%.4Lf", matches[m]->statistics.p_value[1][0], &column);
            add_long_double_to_table(table_line, "%.4Lf", matches[m]->statistics.p_value[0][1], &column);
            add_long_double_to_table(table_line, "%.4Lf", matches[m]->statistics.p_value[1][1], &column);
#endif
            
            strcat(table_line, "\n");
            if (table_fp) {
                fputs(table_line, table_fp);
            }

            n_output_tab++;
            
            if (csv_fp) {
                if ((!(matches[m]->flags & FLAG_IGNORE)) && (!(matches[m]->flags & FLAG_REPEAT))) {
                    sprintf(csv_line, "%d,%d,%d,%s,%d,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n",
                            rank_count,
                            matches[m]->match_number,
                            matches[m]->number_of_paths,
                            convert_type_for_csv(matches[m]->statistics.type, tmp_string),
                            matches[m]->longest_contig,
                            matches[m]->paths[0]->mid,
                            matches[m]->paths[1]->mid,
                            matches[m]->statistics.coverage_complete[0][0] ? matches[m]->statistics.coverage_av[0][0]:0.0,
                            matches[m]->statistics.coverage_complete[1][0] ? matches[m]->statistics.coverage_av[1][0]:0.0,
                            matches[m]->statistics.coverage_complete[0][1] ? matches[m]->statistics.coverage_av[0][1]:0.0,
                            matches[m]->statistics.coverage_complete[1][1] ? matches[m]->statistics.coverage_av[1][1]:0.0,
                            matches[m]->statistics.coverage_complete[0][0] ? matches[m]->statistics.coverage_pc[0][0]:0.0,
                            matches[m]->statistics.coverage_complete[1][0] ? matches[m]->statistics.coverage_pc[1][0]:0.0,
                            matches[m]->statistics.coverage_complete[0][1] ? matches[m]->statistics.coverage_pc[0][1]:0.0,
                            matches[m]->statistics.coverage_complete[1][1] ? matches[m]->statistics.coverage_pc[1][1]:0.0,
                            matches[m]->statistics.coverage_difference[0],
                            matches[m]->statistics.coverage_difference[1],
                            matches[m]->statistics.q_total);
                    fputs(csv_line, csv_fp);
                    n_output_csv++;
                }
            }
            
            rank_count++;                    
        }
    }
    if (table_fp) {
        fclose(table_fp);
    }
    
    if (csv_fp) {
        fclose(csv_fp);
    }
    
    if (table_fp) {
        log_and_screen_printf("%d matches output to table file.\n", n_output_tab);
    }
    
    if (csv_fp) {
        log_and_screen_printf("%d matches output to CSV file.\n", n_output_csv);
    }

    log_write_timestamp(0);
    log_printf(" ENDED Output rank table\n");
}

/*----------------------------------------------------------------------*
 * Function: get_formatted_contig                                       *
 * Purpose:  Make a contig string with lower case letters representing  *
 *           differences.                                               *
 * Params:   out -> char buffer to store string.                        *
 *           m = match number.                                          *
 *           p = path number.                                           *
 * Returns:  Pointer to string.                                         *
 *----------------------------------------------------------------------*/
char* get_formatted_contig(char* out, int m, int p)
{
    int all_same_length = 1;
    int all_same;
    int i,j;
    
    for (i=1; i<matches[m]->number_of_paths; i++) {
        if (matches[m]->paths[i]->mid != matches[m]->paths[0]->mid) {
            all_same_length = 0;
            break;
        }
    }
    
    if (all_same_length) {
        for (i=0; i<matches[m]->paths[p]->mid; i++) {
            all_same = 1;
            for (j=1; j<matches[m]->number_of_paths; j++) {
                if (matches[m]->paths[j]->contig_mid[i] != matches[m]->paths[0]->contig_mid[i]) {
                    all_same = 0;
                    break;
                }
            }
            if (all_same) {
                out[i] = matches[m]->paths[p]->contig_mid[i];
            } else {
                out[i] = tolower(matches[m]->paths[p]->contig_mid[i]);
            }
        }
        out[i] = 0;
    } else {
        make_lower_case_copy(out, matches[m]->paths[p]->contig_mid);
    }       
    return out;
}

/*----------------------------------------------------------------------*
 * Function: output_rank_contigs                                        *
 * Purpose:  Output ranked list of contigs                              *
 * Params:   None.                                                      *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void output_rank_contigs(void)
{
    FILE *fp;
    char header_line[1024];
    char *formatted_seq;
    char previous_type[64];
    char modified_type[64];
    int n_output = 0;
    int rank_count = 1;
    int m, p;
    
    log_newline();
    log_write_timestamp(0);
    log_printf(" START Output rank contigs\n");

    printf("\nWriting rank contig file...\n");
    previous_type[0] = 0;

    formatted_seq = malloc(max_line_length);
    if (!formatted_seq) {
        printf("Error: Insufficient memory to allocate space to store line.\n");
        exit(1);
    }       
    
    fp = fopen(rank_contig_filename, "w");
    for (m=0; m<number_of_matches; m++) {
        if (strcmp(previous_type, matches[m]->statistics.type) != 0) {
            strcpy(previous_type, matches[m]->statistics.type);
            rank_count = 1;
        }

        if (matches[m]->number_of_paths > 0) {
            if (strlen(matches[m]->statistics.type) > 2) {
                sprintf(modified_type, "%c_%s", matches[m]->statistics.type[0], matches[m]->statistics.type+2);
            } else {
                strcpy(modified_type, matches[m]->statistics.type);
            }
            
            for (p=0; p<matches[m]->number_of_paths; p++) {
                sprintf(header_line, ">Type_%s_Rank_%d_Match_%d_Path_%d\n", modified_type, rank_count, matches[m]->match_number, p);
                fputs(header_line, fp);
                fputs(matches[m]->paths[p]->contig_pre, fp);
                get_formatted_contig(formatted_seq, m, p);
                fputs(formatted_seq, fp);
                fputs(matches[m]->paths[p]->contig_post, fp);
                fputs("\n", fp);
            }

            rank_count++;
        }
        
        n_output++;
    }
    
    fclose(fp);
    log_and_screen_printf("%d matches output.\n", n_output);

    log_write_timestamp(0);
    log_printf(" ENDED Output rank contigs\n");
}

/*----------------------------------------------------------------------*
 * Function: look_for_quality_scores_in_fastq                           *
 * Purpose:  Look through all kmers in fastq files and see if they are  *
 *           in our hash table of first kmers in bubbles. If match is   *
 *           found, store the quality scores.                           *
 * Params:   input_filename -> pointer to filename of files             *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void look_for_quality_scores_in_fastq(char* input_filename)
{
    Sequence * seq;
    FILE * fp_fnames;
    FILE * fp_file;
    int count_file   = 0;
    char filename[1024];
    int max_read_length = 2000;
    long long entry_length;
    
    log_newline();
    log_write_timestamp(0);
    log_printf(" START Look for quality scores in FASTQ\n");
    
    printf("\nLooking at quality scores...\n");
    
    // Open file of file names
    fp_fnames= fopen(input_filename, "r");
    if (!fp_fnames) {
        printf("Can't open file of files\n");
        exit(1);
    }
    
    // Preallocate the memory used to read the sequences
    seq = malloc(sizeof(Sequence));
    if (seq == NULL) {
        fputs("Out of memory trying to allocate Sequence\n",stderr);
        exit(1);
    }
    // TODO: Always defaults to offset 33 (Sanger), but for completeness, should this be a parameter?
    alloc_sequence(seq, max_read_length, max_line_length, 33);
    
    // max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
    int max_windows = max_read_length/(kmer_size+1);
    
    // Number of possible kmers in a 'perfect' read
    int max_kmers   = max_read_length-kmer_size+1;
    
    // Preallocate the space of memory used to keep the sliding_windows.
    // NB: this space of memory is reused for every call 
    // -- with the view to avoid memory fragmentation
    // NB: this space needs to preallocate memory for orthogonal situations: 
    //    * a good read -> few windows, many kmers per window
    //    * a bad read  -> many windows, few kmers per window
    KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
    if (windows == NULL) {
        fputs("Out of memory trying to allocate a KmerArraySet",stderr);
        exit(1);
    }  
    
    // Allocate memory for the sliding windows 
    binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);
    
    // For each file 
    while (!feof(fp_fnames)) {
        long long count_bad_reads = 0;
        short colour = 0;
        fscanf(fp_fnames, "%s %hd\n", filename, &colour);
        
        log_printf("Opening file %s\n", filename);
        
        fp_file = fopen(filename, "r");
        if (fp_file == NULL) {
            printf("cannot open file:%s\n",filename);
            exit(1);
        }
        
        long long seq_length = 0;
        count_file++;
        
        while ((entry_length = read_sequence_from_fastq(fp_file, seq, max_read_length))) {
            int i, j;
            
            seq_length += (long long) entry_length;
            
            char quality_cut_off = 0;
            int nkmers = get_sliding_windows_from_sequence_with_quality(seq->seq, seq->qual, entry_length, quality_cut_off, kmer_size, windows, max_windows, max_kmers);
            
            if (nkmers == 0) {
                count_bad_reads++;
            } else {
                // For each window
                for(i=0;i<windows->nwindows;i++) {
                    KmerSlidingWindow * current_window = &(windows->window[i]);
                    
                    // For each kmer in window
                    for(j=0;j<current_window->nkmers;j++) { 
                        BinaryKmer tmp_kmer;
                        
                        Element *e = hash_table_find(element_get_key(&current_window->kmer[j], kmer_size, &tmp_kmer), hash_table);

                        if (e) {
                            element_add_quality_string(e, colour, current_window->quality_strings[j].quality);
                        }
                    }
                }
            }
        }
                
        printf("- File: %i Filename: %s Colour: %d\n", count_file, filename, colour);
    }

    log_write_timestamp(0);
    log_printf(" ENDED Look for quality scores in FASTQ\n");
}

/*----------------------------------------------------------------------*
 * Function: display_sequence_and_coverage                              *
 * Purpose:  Display sequence, with coverage aligned underneath.        *
 * Params:   m = match to show                                          *
 *           width = width of screen                                    *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void display_sequence_and_coverage(int m, int width)
{
    int c, i, j, pos, p;
    int start, end;
    char *line;
    FILE *fp = fopen(fasta_filename, "r");
    
    line = malloc(max_line_length);
    if (!line) {
        printf("Error: Insufficient memory to allocate space to store line.\n");
        exit(1);
    }    
    
    if (!fp) {
        printf("Error: can't open %s\n", fasta_filename);
        exit(1);
    }
    
    for (p=0; p<matches[m]->number_of_paths; p++) {
        fseek(fp, matches[m]->paths[p]->contig_file_offset, SEEK_SET);
        if (!fgets(line, max_line_length, fp)) {
            printf("Error: couldn't read line.\n");
        }
        
        pos = 0;
        start = 0;
        for (i=0; i<matches[m]->paths[p]->contig_length; i++) {
            printf(" %c  ", line[i]);
            pos+=4;
            if ((i == matches[m]->paths[p]->contig_length-1) || ((pos+4) > width)) {
                pos=0;
                printf("\n");
                end = i;
                for (c=0; c<number_of_colours; c++) {
                    for (j=start; j<=end; j++) {
                        printf("%3d ", matches[m]->paths[p]->colour_coverage[c]->coverage[j]);
                    }
                    printf("\n");
                }
                start = i+1;
                printf("\n");
            }
        }
        printf("\n");
    }
    
    fclose(fp);
}

/*----------------------------------------------------------------------*
 * Function: parse_expected_coverage_string                             *
 * Purpose:  Get values out of an expected coverage string in the       *
 *           options file.                                              *
 * Params:   s -> the coverage string, eg. "1,10,50,50"                 *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void parse_expected_coverage_string(char *s)
{
    char* p;
    int colour, tolerance, path;

    p = strtok(s, ",");
    if (p) {
        colour=atoi(p);
        if ((colour >=0) && (colour < MAX_COLOURS)) {
            p = strtok(NULL, ",");
            tolerance = atoi(p);
            if ((tolerance >=0) && (tolerance <=100)) {
                path = 0;
                while (p) {
                    p = strtok(NULL, ",");
                    if (p) {
                        expected_coverage_pc[path++][colour] = atoi(p);                                 
                    }
                };
                coverage_pc_tolerance[colour] = tolerance;
            } else {
                printf("Error: tolerance out of range in expected coverage file.\n");
            }
        } else {
            printf("Error: colour out of range in expected coverage file.\n");
        }
    }                        
}

/*----------------------------------------------------------------------*
 * Function: read_options                                               *
 * Purpose:  Read options file                                          *
 * Params:   None                                                       *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void read_options_file(void)
{
    FILE *fp;
    char line[1024];
    char* param;
    char* value;
    
    log_newline();
    log_write_timestamp(0);
    log_printf(" START Read options file\n");

    if (!options_filename) {
        return;
    }
    
    fp = fopen(options_filename, "r");
    if (!fp) {
        printf("Error: Can't open options file.\n");
        exit(1);
    }

    while (!feof(fp)) {
        while (fgets(line, 1024, fp)) {
            if (line[0] != '#') {
                make_upper_case(line);
                param = strtok(line, " ");
                value = strtok(NULL, "\"");
                if (strcmp(param, "EXPECTEDCOVERAGE") == 0) {
                    if (value) {
                        parse_expected_coverage_string(value);
                    } else {
                        log_and_screen_printf("Warning: Bad value for parameter ExpectedCoverage.\n");
                    }                    
                } else if (strcmp(param, "MINIMUMCONTIGSIZE") == 0) {
                    if (value) {
                        minimum_contig_size = atoi(value);
                    }
                } else {
                    log_and_screen_printf("Warning: Unknown parameter (%s) in options file.\n", param);
                }
            }
        }
    }
    
    log_write_timestamp(0);
    log_printf(" ENDED Read options file\n");
}



/*----------------------------------------------------------------------*
 * Function: output_specified_match_contigs                             *
 * Purpose:  Read a file of match numbers and output contigs for those  *
 *           matches.                                                   *
 * Params:   in -> input filename.                                      *
 *           out -> output filename.                                    *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void output_specified_match_contigs(char *in, char* out)
{
    FILE* fp_in;
    FILE* fp_out;
    int m, p;
    int n_output = 0;
    char header_line[1024];
    char formatted_seq[4096];

    log_newline();
    log_write_timestamp(0);
    log_printf(" START Output specified match contigs\n");

    fp_in = fopen(in, "r");
    if (!fp_in) {
        printf("Error: Can't open specified matches filename.\n");
        exit(1);
    }

    fp_out = fopen(out, "w");
    if (!fp_out) {
        printf("Error: Can't open specified matches output ilename.\n");
        exit(1);
    }
    
    while (!feof(fp_in))
    {
        fscanf(fp_in, "%d\n", &m);

        printf("Looking for %d from %d\n",m, number_of_matches);
        if ((m>0) && (m<=number_of_matches)) {
            sprintf(header_line, ">Match_%d_Position_%d_", matches[m]->match_number, matches[m]->paths[0]->pre+1);

            for (p=0; p<matches[m]->number_of_paths; p++) {
                if (p > 0) {
                    strcat(header_line, "_or_");
                }
                sprintf(header_line + strlen(header_line), "%c", matches[m]->paths[p]->contig_mid[0]);
            }
                                
            fputs(header_line, fp_out);
            fputs("\n", fp_out);
            fputs(matches[m]->paths[0]->contig_pre, fp_out);
            get_formatted_contig(formatted_seq, m, 0);
            formatted_seq[0] = 'N';
            fputs(formatted_seq, fp_out);
            fputs(matches[m]->paths[0]->contig_post, fp_out);
            fputs("\n", fp_out);

            n_output++;             
        }
    }
        
    fclose(fp_in);
    fclose(fp_out);
    
    printf("%d contigs output.\n", n_output);

    log_write_timestamp(0);
    log_printf(" ENDED Output specified match contigs\n");
}

/*----------------------------------------------------------------------*
 * Function: check_inverted_duplicate                                   *
 * Purpose:  Called when there is a suspicion that ma is a reverse      *
 *           complement of mb.                                          *
 * Params:   ma -> first match to check                                 *
 *           mb -> second match to check                                *
 * Returns:  1 if duplicate, 0 if not                                   *
 *----------------------------------------------------------------------*/
int check_inverted_duplicate(Match* ma, Match* mb)
{
    int n_paths = ma->number_of_paths;
    int p, q;
    char* seq_a;
    char* rev_a;
    char* seq_b;
    int n_matches = 0;
    int match_found = 0;
    
    for (p=0; p<n_paths; p++) {
        seq_a = malloc(ma->paths[p]->pre + ma->paths[p]->mid + ma->paths[p]->post + 1);
        rev_a = malloc(ma->paths[p]->pre + ma->paths[p]->mid + ma->paths[p]->post + 1);
        
        if (!seq_a || !rev_a) {
            printf("Error: Out of memory checking inverted duplicates.\n");
            exit(1);
        }
        
        strcpy(seq_a, ma->paths[p]->contig_pre);
        strcat(seq_a, ma->paths[p]->contig_mid);
        strcat(seq_a, ma->paths[p]->contig_post);
        make_reverse_compliment(seq_a, rev_a);
                
        match_found = 0;
        for (q=0; q<n_paths; q++) {
            seq_b = malloc(mb->paths[q]->pre + mb->paths[q]->mid + mb->paths[q]->post + 1);
            if (!seq_b) {
                printf("Error: Out of memory checking inverted duplicates.\n");
                exit(1);
            }
            
            strcpy(seq_b, mb->paths[q]->contig_pre);
            strcat(seq_b, mb->paths[q]->contig_mid);
            strcat(seq_b, mb->paths[q]->contig_post);
            
            if (strcmp(seq_b, rev_a) == 0) {
                match_found = 1;
            }
            
            free(seq_b);
            
            if (match_found) {
                break;
            }
        }
        
        free(seq_a);
        free(rev_a);
        
        if (match_found) {
            n_matches++;
        } else {
            break;
        }
    }
    
    return (n_matches == n_paths ? 1:0);
}

/*----------------------------------------------------------------------*
 * Function: deduplicate_matches                                        *
 * Purpose:  Where one match is the reverse compliment of another, mark *
 *           the match (and all paths) as repeats.                      *
 * Params:   None                                                       *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void deduplicate_matches(void)
{
    Match* key;
    int m;
        
    key = malloc(sizeof(Match));
    if (!key) {
        printf("Error: couldn't get memory for MatchPath structure.\n");
        exit(1);
    }                 
    key->paths[0] = malloc(sizeof(MatchPath));
    if (!key->paths[0]) {
        printf("Error: couldn't get memory for MatchPath structure.\n");
        exit(1);
    }
    key->number_of_paths = 1;
    key->paths[0]->first_contig_kmer = malloc(kmer_size+1);
    if (!key->paths[0]->first_contig_kmer) {
        printf("Error: couldn't get memory for MatchPath structure.\n");
        exit(1);
    }
    
    log_newline();
    log_write_timestamp(0);
    log_printf(" START Sorting for deduplicates\n");
    
    printf("\nSorting matches for deduplicates...\n");  
    qsort(matches, number_of_matches, sizeof(Match*), match_alpha_cmp);
    
    log_write_timestamp(0);
    log_printf(" ENDED Sorting for deduplicates\n");

    log_newline();
    log_write_timestamp(0);
    log_printf(" START Deduplicating\n");

    for (m=0; m<number_of_matches; m++)
    {
        if (!(matches[m]->flags & FLAG_REPEAT)) {
            // Get last kmer, reverse it, then that is what we are searching for
            make_reverse_compliment(matches[m]->paths[0]->last_contig_kmer, key->paths[0]->first_contig_kmer);
                                        
            // Try and find it
            Match** mp = bsearch(&key, matches, number_of_matches, sizeof(Match*), match_alpha_cmp);
            if (mp) {
                Match* ma = *mp;
                if (ma->match_number != matches[m]->match_number) {
                    if (ma->number_of_paths == matches[m]->number_of_paths) {
                        if (check_inverted_duplicate(matches[m], ma) == 1) {
                            log_printf("Warning: Match %d is a reverse compliment copy of %d. Marked %d as a repeat.\n", matches[m]->match_number, ma->match_number, ma->match_number);
                            set_flag(&(ma->flags), FLAG_REPEAT);
                            n_deduplicates++;
                        }
                    }
                }
            }
        }
    }
    
    free(key->paths[0]->first_contig_kmer);
    free(key->paths[0]);
    free(key);
    
    log_printf("%d duplicates detected and marked.\n", n_deduplicates);
    log_write_timestamp(0);
    log_printf(" ENDED Deduplicating\n");    
}

/*----------------------------------------------------------------------*
 * Function: parse_string                                               *
 * Purpose:  Return a string parameter (command line parsing)           *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 *           i -> pointer to argument number counter                    *
 * Returns:  a string and updates i                                     *
 *----------------------------------------------------------------------*/
char* parse_string(int argc, char* argv[], int* i)
{
    char* token = 0;
    
    if (strlen(argv[*i]) > 2)
        token = argv[*i] + 2;
    else if (*i < (argc-1)) {
        *i = *i + 1;
        token = argv[*i];
    }
    
    if ((token) && (token[0] == '-'))
        token = 0;
    
    return token;
}

/*----------------------------------------------------------------------*
 * Function: parse_int                                                  *
 * Purpose:  Return an integer parameter (command line parsing)         *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 *           i -> pointer to argument number counter                    *
 * Returns:  a number and updates i                                     *
 *----------------------------------------------------------------------*/
int parse_int(int argc, char* argv[], int* i)
{
    int v = -1;
    char* token = parse_string(argc, argv, i);
    
    if (token)
        v = atoi(token);
    
    return v;
}

/*----------------------------------------------------------------------*
 * Function: parse_command_line_args                                    *
 * Purpose:  Deal with command line arguments                           *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void parse_command_line_args(int argc, char* argv[])
{
    int i = 1;

    if (argc < 4)
    {
        printf("\nbubbleparse\n\n");
        printf("Syntax: bubbleparse [-f filename] [-i filename] [-k int] [options]\n");
        printf("where [-f filename] specifies the base filename of the files to parse\n");
        printf("      [-i filename] specifies the name of a file of files containing the original read files\n");
        printf("      [-k int] specifies kmer size\n");
        printf("      [-t filename] specifies the filename for a plain text ranked table of SNPs\n");
        printf("      [-c filename] specifies the filename for a CSV ranked table of SNPs\n");
        printf("      [-r filename] specifies the filename for a ranked contig FASTA file\n");
        printf("      [-l filename] specifies the filename for a rank log file\n");
        //printf("      [-m filename] to cancel the full parse, but to output match contigs contained in file specified\n");  
        printf("      [-d filename] to output a log file\n");
        printf("      [-o filename] specifies an options file\n");
        printf("      [-p int] specifies the FASTQ quality score offset (default 64)\n");
        printf("      [-x] to remove reverse compliment duplicates due to X nodes\n");
        //printf("      [-y int] specifies the maximum line length in kb (default 8192).\n");
        printf("\n");
        exit(1);
    }

    while(i < argc)
    {
        char* parameter = argv[i];
        if (parameter[0] == '-')
        {
            switch (parameter[1])
            {
                case 'c':
                    rank_csv_filename = parse_string(argc, argv, &i);
                    if (!rank_csv_filename) {
                        printf("Error: Invalid filename for -c parameter.\n");
                        exit(1);
                    }
                    break;
                case 'd':
                    debug_log_filename = parse_string(argc, argv, &i);
                    if (!debug_log_filename) {
                        printf("Error: Invalid filename for -d parameter.\n");
                        exit(1);
                    }
                case 'e':
                    expected_filename = parse_string(argc, argv, &i);
                    if (!expected_filename) {
                        printf("Error: Invalid filename for -e parameter.\n");
                        exit(1);
                    }
                    break;
                case 'f':
                    base_filename = parse_string(argc, argv, &i);
                    if (!base_filename) {
                        printf("Error: Invalid filename for -f parameter.\n");
                        exit(1);
                    }
                    break;
                case 'i':
                    file_of_filenames = parse_string(argc, argv, &i);
                    use_quality_scores = 1;
                    if (!file_of_filenames) {
                        printf("Error: Invalid filename for -i parameter.\n");
                        exit(1);
                    }
                    break;
                case 'k':
                    kmer_size = parse_int(argc, argv, &i);
                    if (kmer_size >= (NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32)) {
                        printf("Error: This version compiled only for kmers up to %i.\n", (NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32)-1);
                        exit(1);
                    }
                    break;
                case 'l':
                    rank_log_filename = parse_string(argc, argv, &i);
                    if (!rank_log_filename) {
                        printf("Error: Invalid filename for -l parameter.\n");
                        exit(1);
                    }
                    break;
                case 'm':
                    match_selection_filename = parse_string(argc, argv, &i);
                    if (!match_selection_filename) {
                        printf("Error: Must specify a filename for -m parameter.\n");
                        exit(1);
                    }
                    break;
                case 'o':
                    options_filename = parse_string(argc, argv, &i);
                    if (!options_filename) {
                        printf("Error: Invalid filename for -o parameter.\n");
                        exit(1);
                    }
                    break;
                case 'p':
                    quality_offset = parse_int(argc, argv, &i);
                    if ((quality_offset != 64) && (quality_offset != 33) && (quality_offset != 59)) {
                        printf("Warning: Unrecognised FASTQ quality offset... running anyway...\n");
                    }
                    break;
                case 'r':
                    rank_contig_filename = parse_string(argc, argv, &i);
                    if (!rank_contig_filename) {
                        printf("Error: Invalid filename for -r parameter.\n");
                        exit(1);
                    }
                    break;
                case 't':
                    rank_table_filename = parse_string(argc, argv, &i);
                    if (!rank_table_filename) {
                        printf("Error: Invalid filename for -t parameter.\n");
                        exit(1);
                    }
                    break;
                case 'x':
                    deduplicate = TRUE;
                    break;
                case 'y':
                    max_line_length = parse_int(argc, argv, &i) * 1024;
                    printf("Max line length set to %d bytes.\n", max_line_length);
                    break;
                default:
                    printf("Error: Invalid parameter %c\n", parameter[1]);
                    exit(1);
            }
        }
        
        i++;
    }
    
    if ((kmer_size < 1) || (kmer_size > 255)) {
        printf("Error: Invalid kmer size.\n");
        exit(1);
    }
    
    //if (!file_of_filenames) {
    //    printf("Error: You must specify an input file of files.\n");
    //    exit(1);
    //}
    
    if (!base_filename) {
        printf("Error: You must specify a .gv output file.\n");
        exit(1);
    }   
    
    if (!rank_log_filename && !match_selection_filename && !rank_table_filename && !rank_csv_filename && !rank_contig_filename) {
        printf("Error: You must specify at least one output file.\n");
        exit(1);
    }
}

/*----------------------------------------------------------------------*
 * Debug/test functions                                                 *
 *----------------------------------------------------------------------*/
#ifdef DEBUG_MAKE_TEST_DATA
int randint(int max)
{
    double d = (double)rand() / (double)RAND_MAX;   
    return (int)((double)max*d);
}

void make_random_quality_string(char* s)
{
    int i;
    for (i=0; i<kmer_size; i++) {
        s[i] = quality_offset + randint(40) + 20;
    }
    s[kmer_size] = 0;
}

void make_test_quality_data(void)
{   
    void make_quality_strings(dBNode * node) {
        int c, i;
        char qs[kmer_size+1];
        for (c=0; c<number_of_colours; c++) {
            QualityStringArray *qa = &node->quality_string_arrays[c];
            
            assert(qa != 0);            
            assert(qa->number_of_strings >= 0);
            assert(qa->limit >= qa->number_of_strings);
            
            for (i=qa->number_of_strings; i<qa->limit; i++) {
                make_random_quality_string(qs);
                element_add_quality_string(node, c, qs);
            }
        }
    }
    
    printf("\n\n\n*** Warning ***\n Debug mode: Creating test quality data...\n\n\n");
    hash_table_traverse(&make_quality_strings, hash_table);
}
#endif

void print_parameter_values(void)
{
    int match_size = sizeof(Match)+sizeof(RankingStats)+((sizeof(MatchPath)+kmer_size+1)*MAX_PATHS_PER_MATCH);
    
    printf("\nbubbleparse\n\n");
    printf("                Max k: %i\n", (NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32)-1);
    printf("     Input fasta file: %s\n", fasta_filename);
    printf("  Input coverage file: %s\n", coverage_filename);
    if (use_quality_scores) {
        printf("  Input file of files: %s\n", file_of_filenames);
    }
    if (match_selection_filename) {
        printf(" Match selection file: %s\n", match_selection_filename);
        printf("Selected matches file: %s\n", selected_matches_filename);
    }
    if (rank_table_filename) {
        printf("      Rank table file: %s\n", rank_table_filename);
    }
    if (rank_csv_filename) {
        printf("        Rank CSV file: %s\n", rank_csv_filename);
    }
    if (rank_contig_filename) {
        printf("     Rank contig file: %s\n", rank_contig_filename);
    }
    if (rank_log_filename) {
        printf("        Rank log file: %s\n", rank_log_filename);
    }    
    printf("            Kmer size: %d\n", kmer_size);
    printf("           Match size: %i\n", match_size);
    printf("        c0 expected %%: %.2f/%.2f\n", expected_coverage_pc[0][0], expected_coverage_pc[1][0]);
    printf("         c0 tolerance: %.2f\n", coverage_pc_tolerance[0]);
    printf("        c1 expected %%: %.2f/%.2f\n", expected_coverage_pc[0][1], expected_coverage_pc[1][1]);
    printf("         c1 tolerance: %.2f\n", coverage_pc_tolerance[1]);
    printf("  Minimum contig size: %d\n", minimum_contig_size); 
    printf(" Quality score offset: %i", quality_offset);

    if (quality_offset == 33) {
        printf(" (Sanger format)\n");
    } else if (quality_offset == 64) {
        printf(" (Solexa/Illumina format)\n");
    } else {
        printf(" (Unknown format)\n");
    }
}

void handler(int sig) {
    void *array[10];
    size_t size;
    
    // get void*'s for all entries on the stack
    size = backtrace(array, 10);
    
    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, 2);
    log_printf("\nError: signal %d:\n", sig);
    exit(1);
}

/*----------------------------------------------------------------------*
 * Function: main                                                       *
 * Purpose:  Entry point to program.                                    *
 * Params:   argc = number of arguments                                 *
 *           argv -> array of arguments                                 *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
int main (int argc, char * argv[])
{
    char time_string[64];
    int i;
    int n=16;
    int b=100;

    // Install signal handler for catching segmentation faults
    signal(SIGSEGV, handler);    
    
    // Parse command line args
    parse_command_line_args(argc, argv);
    
    // Debugging log
    if (debug_log_filename) {
        log_start(debug_log_filename);
    }
    
    sprintf(fasta_filename, "%s.fasta", base_filename);
    sprintf(coverage_filename, "%s.coverage", base_filename);
    sprintf(selected_matches_filename, "%s.selected", base_filename);
    
    matches = malloc(max_matches * sizeof(Match*));
    for (i=0; i<max_matches; i++) matches[i] = 0;
        
    expected_coverage_pc[0][0] = 50;
    expected_coverage_pc[1][0] = 50;
    expected_coverage_pc[0][1] = 100;
    expected_coverage_pc[1][1] = 0; 
    coverage_pc_tolerance[0] = 10;
    coverage_pc_tolerance[1] = 10;
    
    read_options_file();    
    print_parameter_values();

    make_time_string(time_string);
    printf("\nStarting at %s...\n\n", time_string);
    
    // Create hash table
    printf("Creating hash table (entry size %libytes, total size %iMb)\n", sizeof(Element), (int)(pow(2, n)*b*sizeof(Element)/1024/1024));
    for (i=0; i<max_matches; i++) {
        matches[i] = 0;
    }
    hash_table = hash_table_new(n, b, 10, kmer_size); 
        
    read_fasta_file(fasta_filename);
    read_coverage_file(coverage_filename);

    if (match_selection_filename) {
        output_specified_match_contigs(match_selection_filename, selected_matches_filename);
    } else {        
        generate_overall_stats();
        
        if (use_quality_scores) {
            look_for_quality_scores_in_fastq(file_of_filenames);
        }
        
        #ifdef DEBUG_MAKE_TEST_DATA
        make_test_quality_data();
        #endif
        
        score_matches();
        
        if (deduplicate) {
            deduplicate_matches();
        }
        
        rank_matches();
        
        if (rank_table_filename || rank_csv_filename) {
            output_rank_table();
        }
        
        if (rank_contig_filename) {
            output_rank_contigs();
        }
    }
        
    make_time_string(time_string);
    printf("\nFinished at %s.\n", time_string);
    log_printf("\nFinished.\n");

    return 0;
}
