#include <inttypes.h>
#include "kstring.h"

/*****************
 * Strand macros *
 *****************/
#define GTF_IGNORE_STRAND     0
#define GTF_SAME_STRAND       1
#define GTF_OPPOSITE_STRAND   2
#define GTF_EXACT_SAME_STRAND 3

/***********************
 * Overlap type macros *
 ***********************/
#define GTF_MATCH_ANY     0
#define GTF_MATCH_EXACT   1
#define GTF_MATCH_CONTAIN 2
#define GTF_MATCH_WITHIN  3
#define GTF_MATCH_START   4
#define GTF_MATCH_END     5

typedef struct {
    int32_t key;
    int32_t val;
} Attribute;

/*! @typedef
 @abstract Structure for a single GTF line
 @field	 chrom         Index into the chrom hash table
 @field  source        Index into the source hash table
 @field  feature       Index into the feature hash table
 @field  start         0-based starting position
 @field  end           1-based end position
 @field  score         The score field. A value of DBL_MAX indicates a "."
 @field  strand        0: '+'; 1: '-'; 3: '.'
 @field  frame         0: '0'; 1: '1'; 2: '2'; 3: '.'
 @field  gene_id       Index into the gene_id hash table
 @field  transcript_id Index into the transcript_id hash table
 @discussion Positions are 0-based half open ([start, end)), like BED files.
*/

typedef struct GTFentry {
    int32_t chrom;
    int32_t source;
    int32_t feature;
    uint32_t start;
    uint32_t end;
    double score;
    uint8_t strand:4, frame:4;
    int32_t gene_id;
    int32_t transcript_id;
    uint32_t labelIdx;
    int nAttributes;
    Attribute **attrib;
    struct GTFentry *left, *right;
} GTFentry;

typedef struct {
    kstring_t chrom;
    kstring_t feature;
    kstring_t source;
    uint32_t start;
    uint32_t end;
    double score;
    uint8_t strand:4, frame: 4;
    kstring_t gene;
    kstring_t transcript;
    int nAttributes;
    Attribute **attrib;
} GTFline;

typedef struct GTFnode {
    uint32_t center;
    GTFentry *starts, *ends;
    struct GTFnode *left, *right;
} GTFnode;

typedef struct {
    int32_t chrom;
    uint32_t n_entries;
    void **tree;
} GTFchrom;

typedef struct hashTableElement {
    int32_t val;
    struct hashTableElement *next;
} hashTableElement;

typedef struct {
    uint64_t l, m;
    hashTableElement **elements;
    char **str;
} hashTable;

typedef struct {
    int32_t n_targets, m;
    int balanced;
    hashTable *htChroms;
    hashTable *htSources;
    hashTable *htFeatures;
    hashTable *htAttributes;
    GTFchrom **chroms;
} GTFtree;

typedef struct {
    int32_t l, m;
    GTFentry **overlaps;
    GTFtree *tree;
} overlapSet;

typedef struct {
    int32_t l, m;
    overlapSet **os;
} overlapSetList;

typedef struct {
    int32_t l, m;
    int32_t *IDs;
    uint32_t *cnts;
    hashTable *ht;
} uniqueSet;

//A function that can be applied to all entries in a GTF/BED/etc. file as it's
//being processed. The pointer as input is currently a GTFline *. The return
//value is 0 (ignore entry) or 1 (keep entry).
typedef int (*FILTER_FUNC)(void*);
typedef int (*FILTER_ENTRY_FUNC)(GTFtree *, GTFentry *);

//A function used to compare to GTFentry items to see if the intersect in some
//way (e.g., due to sharing a gene_id). This is used to intersect overlapsets.
typedef int (*COMPARE_FUNC)(GTFentry *, GTFentry *);

//gtf.c
GTFtree * initGTFtree(void);
void destroyGTFtree(GTFtree *t);
void sortGTF(GTFtree *o);
void printGTFtree(GTFtree *t);
void printGTFvineStart(GTFentry *e, const char *chrom, const char *str);
void printGTFvineStartR(GTFentry *e, const char *chrom, const char *str);
int addGTFentry(GTFtree *t, char *chrom, uint32_t start, uint32_t end, uint8_t strand, char *transcriptID, uint32_t labelIDX, double score);
int addEnrichmententry(GTFtree *t, char *chrom, uint32_t start, uint32_t end, uint8_t strand, double score, char *feature);
int hasOverlaps(GTFtree *t, uint32_t *minOverlap);

//hashTable.c
hashTable *initHT(uint64_t size);
void destroyHTelement(hashTableElement *e);
void destroyHT(hashTable *ht);
int32_t addHTelement(hashTable *ht, char *s);
uint64_t hashString(char *s);
int strExistsHT(hashTable *ht, char *s);
int32_t str2valHT(hashTable *ht, char *s);
char *val2strHT(hashTable *ht, int32_t val);
int hasAttribute(GTFtree *t, GTFentry *e, char *str);
char *getAttribute(GTFtree *t, GTFentry *e, char *str); //NULL if the attribute isn't there

//findOverlaps.c
//overlapSet functions
overlapSet *os_init(GTFtree *t);
void os_reset(overlapSet *os);
void os_destroy(overlapSet *os);
overlapSet *os_grow(overlapSet *os);
void os_exclude(overlapSet *os, int i);
void os_requireAttributes(overlapSet *os, char **keys, char **vals, int len);
void os_requireSource(overlapSet *os, char *val);
void os_requireFeature(overlapSet *os, char *val);
overlapSet *os_intersect(overlapSet *os1, overlapSet *os2, COMPARE_FUNC f);
//overlapSetList functions
overlapSetList *osl_init(void);
void osl_reset(overlapSetList *osl);
void osl_destroy(overlapSetList *osl);
void osl_push(overlapSetList *osl, overlapSet *os);
void osl_grow(overlapSetList *osl);
overlapSet *osl_intersect(overlapSetList *osl, COMPARE_FUNC f);
overlapSet *osl_union(overlapSetList *osl);
//uniqueSet functions
void us_destroy(uniqueSet *us);
uint32_t us_cnt(uniqueSet *us, int32_t i);
char *us_val(uniqueSet *us, int32_t i);
//Driver functions
overlapSet * findOverlaps(overlapSet *os, GTFtree *t, char *chrom, uint32_t start, uint32_t end, int strand, int matchType, int strandType, int keepOS, FILTER_ENTRY_FUNC ffunc);
