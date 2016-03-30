#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "gtf.h"

/*******************************************************************************
*
* Comparison functions
*
* These are according to Allen's Interval Algebra
*
*******************************************************************************/
static inline int rangeAny(uint32_t start, uint32_t end, GTFentry *e) {
    if(end <= e->start) return -1;
    if(start >= e->end) return 1;
    return 0;
}

static inline int rangeContains(uint32_t start, uint32_t end, GTFentry *e) {
    if(e->start >= start && e->end <= end) return 0;
    if(e->end < end) return -1;
    return 1;
}

static inline int rangeWithin(uint32_t start, uint32_t end, GTFentry *e) {
    if(start >= e->start && end <= e->end) return 0;
    if(start < e->start) return -1;
    return 1;
}

static inline int rangeExact(uint32_t start, uint32_t end, GTFentry *e) {
    if(start == e->start && end == e->end) return 0;
    if(start < e->start) return -1;
    if(end < e->end) return -1;
    return 1;
}

static inline int rangeStart(uint32_t start, uint32_t end, GTFentry *e) {
    if(start == e->start) return 0;
    if(start < e->start) return -1;
    return 1;
}

static inline int rangeEnd(uint32_t start, uint32_t end, GTFentry *e) {
    if(end == e->end) return 0;
    if(end < e->end) return -1;
    return 1;
}

static inline int exactSameStrand(int strand, GTFentry *e) {
    return strand == e->strand;
}

static inline int sameStrand(int strand, GTFentry *e) {
    if(strand == 3 || e->strand == 3) return 1;
    if(strand == e->strand) return 1;
    return 0;
}

static inline int oppositeStrand(int strand, GTFentry *e) {
    if(strand == 3 || e->strand == 3) return 1;
    if(strand != e->strand) return 1;
    return 0;
}

/*******************************************************************************
*
* OverlapSet functions
*
*******************************************************************************/
overlapSet *os_init(GTFtree *t) {
    overlapSet *os = calloc(1, sizeof(overlapSet));
    assert(os);
    os->tree = t;
    return os;
}

void os_reset(overlapSet *os) {
    int i;
    for(i=0; i<os->l; i++) os->overlaps[i] = NULL;
    os->l = 0;
}

void os_destroy(overlapSet *os) {
    if(os->overlaps) free(os->overlaps);
    free(os);
}

overlapSet *os_grow(overlapSet *os) {
    int i;
    os->m++;
    kroundup32(os->m);
    os->overlaps = realloc(os->overlaps, os->m * sizeof(GTFentry*));
    assert(os->overlaps);
    for(i=os->l; i<os->m; i++) os->overlaps[i] = NULL;

    return os;
}

static void os_push(overlapSet *os, GTFentry *e) {
    if(os->l+1 >= os->m) os = os_grow(os);
    os->overlaps[os->l++] = e;
}

overlapSet *os_dup(overlapSet *os) {
    int i;
    overlapSet *os2 = os_init(os->tree);
    for(i=0; i<os->l; i++) os_push(os2, os->overlaps[i]);
    return os2;
}

void os_exclude(overlapSet *os, int i) {
    int j;
    for(j=i; j<os->l-1; j++) os->overlaps[j] = os->overlaps[j+1];
    os->overlaps[--os->l] = NULL;
}

static int os_sortFunc(const void *a, const void *b) {
    GTFentry *pa = *(GTFentry**) a;
    GTFentry *pb = *(GTFentry**) b;

    if(pa->start < pb->start) return -1;
    if(pb->start < pa->start) return 1;
    if(pa->end < pb->end) return -1;
    if(pb->end < pa->end) return 1;
    return 0;
}

static void os_sort(overlapSet *os) {
    qsort((void *) os->overlaps, os->l, sizeof(GTFentry**), os_sortFunc);
}

//Non-existant keys/values will be ignored
void os_requireAttributes(overlapSet *os, char **key, char **val, int len) {
    int i, j, k, filter;
    int32_t keyHash, valHash;

    for(i=0; i<len; i++) {
        if(!os->l) break;

        keyHash = str2valHT(os->tree->htAttributes, key[i]);
        valHash = str2valHT(os->tree->htAttributes, val[i]);
        assert(keyHash>=0);
        assert(valHash>=0);
        for(j=0; j<os->l; j++) {
            filter = 1;
            for(k=0; k<os->overlaps[j]->nAttributes; k++) {
                if(os->overlaps[j]->attrib[k]->key == keyHash) {
                    if(os->overlaps[j]->attrib[k]->val == valHash) {
                        filter = 0;
                        break;
                    }
                }
            }
            if(filter) {
                os_exclude(os, j);
                j--; //os_exclude shifts everything
            }
        }
    }
}

//This is an inefficient implementation. It would be faster to sort according
//to COMPARE_FUNC and then do an O(n) merge.
overlapSet *os_intersect(overlapSet *os1, overlapSet *os2, COMPARE_FUNC f) {
    overlapSet *os = os_init(os1->tree);
    int i, j;

    for(i=0; i<os1->l; i++) {
        for(j=0; j<os2->l; j++) {
            if(f(os1->overlaps[i],os2->overlaps[j]) == 0) {
                os_push(os, os1->overlaps[i]);
                os_exclude(os2, j);
                break;
            }
        }
    }

    return os;
}
/*******************************************************************************
*
* OverlapSetList functions
*
*******************************************************************************/
overlapSetList *osl_init() {
    overlapSetList *osl = calloc(1, sizeof(overlapSetList));
    assert(osl);
    return osl;
}

void osl_reset(overlapSetList *osl) {
    int i;
    for(i=0; i<osl->l; i++) os_destroy(osl->os[i]);
    osl->l = 0;
}

void osl_destroy(overlapSetList *osl) {
    osl_reset(osl);
    if(osl->os) free(osl->os);
    free(osl);
}

void osl_grow(overlapSetList *osl) {
    int i;
    osl->m++;
    kroundup32(osl->m);
    osl->os = realloc(osl->os, osl->m * sizeof(overlapSet*));
    assert(osl->os);
    for(i=osl->l; i<osl->m; i++) osl->os[i] = NULL;
}

void osl_push(overlapSetList *osl, overlapSet *os) {
    if(osl->l+1 >= osl->m) osl_grow(osl);
    osl->os[osl->l++] = os;
}

//The output needs to be destroyed
overlapSet *osl_intersect(overlapSetList *osl, COMPARE_FUNC f) {
    int i;
    if(!osl->l) return NULL;

    overlapSet *osTmp, *os = os_dup(osl->os[0]);
    for(i=1; i<osl->l; i++) {
        osTmp = os_intersect(os, osl->os[i], f);
        os_destroy(os);
        os = osTmp;
        if(os->l == 0) break;
    }
    return os;
}

//Returns 1 if the node is in the overlapSet, otherwise 0.
int os_contains(overlapSet *os, GTFentry *e) {
    int i;
    for(i=0; i<os->l; i++) {
        if(os->overlaps[i] == e) return 1;
    }
    return 0;
}

//This could be made much more efficient
overlapSet *osl_union(overlapSetList *osl) {
    int i, j;
    if(!osl->l) NULL;

    if(!osl->os) return NULL;
    if(!osl->os[0]) return NULL;
    overlapSet *os = os_dup(osl->os[0]);
    for(i=1; i<osl->l; i++) {
        for(j=0; j<osl->os[i]->l; j++) {
            if(!os_contains(os, osl->os[i]->overlaps[j])) {
                os_push(os, osl->os[i]->overlaps[j]);
            }
        }
    }
    return os;
}

/*******************************************************************************
*
* uniqueSet functions
*
*******************************************************************************/
static uniqueSet *us_init(hashTable *ht) {
    uniqueSet *us = calloc(1, sizeof(uniqueSet));
    assert(us);
    us->ht = ht;
    return us;
}

void us_destroy(uniqueSet *us) {
    if(!us) return;
    if(us->IDs) {
        free(us->IDs);
        free(us->cnts);
    }
    free(us);
}

static uniqueSet *us_grow(uniqueSet *us) {
    int i;
    us->m++;
    kroundup32(us->m);
    us->IDs = realloc(us->IDs, us->m * sizeof(int32_t));
    assert(us->IDs);
    us->cnts = realloc(us->cnts, us->m * sizeof(uint32_t));
    assert(us->cnts);
    for(i=us->l; i<us->m; i++) {
        us->IDs[i] = -1;
        us->cnts[i] = 0;
    }

    return us;
}

static void us_push(uniqueSet *us, int32_t ID) {
    if(us->l+1 >= us->m) us = us_grow(us);
    us->IDs[us->l] = ID;
    us->cnts[us->l++] = 1;
}

static void us_inc(uniqueSet *us) {
    assert(us->l<=us->m);
    us->cnts[us->l-1]++;
}

uint32_t us_cnt(uniqueSet *us, int32_t i) {
    assert(i<us->l);
    return us->cnts[i];
}

char *us_val(uniqueSet *us, int32_t i) {
    if(i>=us->l) return NULL;
    return val2strHT(us->ht, us->IDs[i]);
}

/*******************************************************************************
*
* Overlap set count/unique functions
*
*******************************************************************************/
static int int32_t_cmp(const void *a, const void *b) {
    int32_t ia = *((int32_t*) a);
    int32_t ib = *((int32_t*) b);
    return ia-ib;
}

int32_t cntAttributes(overlapSet *os, char *attributeName) {
    int32_t IDs[os->l], i, j, key, last, n = 0;
    if(!strExistsHT(os->tree->htAttributes, attributeName)) return n;

    key = str2valHT(os->tree->htAttributes, attributeName);
    for(i=0; i<os->l; i++) {
        IDs[i] = -1;
        for(j=0; j<os->overlaps[i]->nAttributes; j++) {
            if(os->overlaps[i]->attrib[j]->key == key) {
                IDs[i] = os->overlaps[i]->attrib[j]->val;
                break;
            }
        }
    }
    qsort((void*) IDs, os->l, sizeof(int32_t), int32_t_cmp);

    last = IDs[0];
    n = (last >= 0) ? 1 : 0;
    for(i = 1; i<os->l; i++) {
        if(IDs[i] != last) {
            n++;
            last = IDs[i];
        }
    }
    return n;
}

uniqueSet *uniqueAttributes(overlapSet *os, char *attributeName) {
    if(!os) return NULL;
    if(os->l == 0) return NULL;
    int32_t IDs[os->l], i, j, key, last;
    if(!strExistsHT(os->tree->htAttributes, attributeName)) return NULL;
    uniqueSet *us = us_init(os->tree->htAttributes);

    key = str2valHT(os->tree->htAttributes, attributeName);
    for(i=0; i<os->l; i++) {
        IDs[i] = -1;
        for(j=0; j<os->overlaps[i]->nAttributes; j++) {
            if(os->overlaps[i]->attrib[j]->key == key) {
                IDs[i] = os->overlaps[i]->attrib[j]->val;
                break;
            }
        }
    }
    qsort((void*) IDs, os->l, sizeof(int32_t), int32_t_cmp);

    last = -1;
    for(i=0; i<os->l; i++) {
        if(IDs[i] != last || last < 0) {
            us_push(us, IDs[i]);
            last = IDs[i];
        } else {
            us_inc(us);
        }
    }

    if(us->l) return us;
    us_destroy(us);
    return NULL;
}

/*******************************************************************************
*
* Node iterator functions
*
*******************************************************************************/
//bit 1: go left, bit 2: go right (a value of 3 is then "do both")
static int centerDirection(uint32_t start, uint32_t end, GTFnode *n) {
    if(n->center >= start && n->center < end) return 3;
    if(n->center < start) return 2;
    return 1;
}

static int matchingStrand(GTFentry *e, int strand, int strandType) {
    if(strandType == GTF_IGNORE_STRAND) return 1;

    if(strandType == GTF_SAME_STRAND) {
        return sameStrand(strand, e);
    } else if(strandType == GTF_OPPOSITE_STRAND) {
        return oppositeStrand(strand, e);
    } else if(strandType == GTF_EXACT_SAME_STRAND) {
        return exactSameStrand(strand, e);
    }

    fprintf(stderr, "[matchingStrand] Unknown strand type %i. Assuming a match.\n", strandType);
    return 1;
}

static void filterStrand(overlapSet *os, int strand, int strandType) {
    int i;

    if(strandType == GTF_IGNORE_STRAND) return;

    for(i=os->l-1; i>=0; i--) {
        if(strandType == GTF_SAME_STRAND) {
            if(!sameStrand(strand, os->overlaps[i])) os_exclude(os, i);
        } else if(strandType == GTF_OPPOSITE_STRAND) {
            if(!oppositeStrand(strand, os->overlaps[i])) os_exclude(os, i);
        } else if(strandType == GTF_EXACT_SAME_STRAND) {
            if(!exactSameStrand(strand, os->overlaps[i])) os_exclude(os, i);
        }
    }
}

static void pushOverlaps(overlapSet *os, GTFtree *t, GTFentry *e, uint32_t start, uint32_t end, int comparisonType, int direction, FILTER_ENTRY_FUNC ffunc) {
    int dir;
    int keep = 1;
    if(!e) return;

    if(ffunc) keep = ffunc(t, e);

    switch(comparisonType) {
    case GTF_MATCH_EXACT :
        if((dir = rangeExact(start, end, e)) == 0) {
            if(keep) os_push(os, e);
        }
        break;
    case GTF_MATCH_WITHIN :
        if((dir = rangeAny(start, end, e)) == 0) {
            if(keep) if(rangeWithin(start, end ,e) == 0) os_push(os, e);
        }
        break;
    case GTF_MATCH_CONTAIN :
        if((dir = rangeAny(start, end, e)) == 0) {
            if(keep) if(rangeContains(start, end, e) == 0) os_push(os, e);
        }
        break;
    case GTF_MATCH_START :
        if((dir = rangeStart(start, end, e)) == 0) {
            if(keep) os_push(os, e);
        }
        break;
    case GTF_MATCH_END :
        if((dir = rangeEnd(start, end, e)) == 0) {
            if(keep) os_push(os, e);
        }
        break;
    default :
        if((dir = rangeAny(start, end, e)) == 0) {
            if(keep) os_push(os, e);
        }
        break;
    }

    if(direction) {
        if(dir > 0) return;
        pushOverlaps(os, t, e->right, start, end, comparisonType, direction, ffunc);
    } else {
        if(dir < 0) return;
        pushOverlaps(os, t, e->left, start, end, comparisonType, direction, ffunc);
    }
}

static int32_t countOverlapsEntry(GTFtree *t, GTFentry *e, uint32_t start, uint32_t end, int strand, int matchType, int strandType, int direction, int32_t max, FILTER_ENTRY_FUNC ffunc) {
    int dir;
    int32_t cnt = 0;
    if(!e) return cnt;
    
    switch(matchType) {
    case GTF_MATCH_EXACT :
        if((dir = rangeExact(start, end, e)) == 0) {
            cnt = 1;
        }
        break;
    case GTF_MATCH_WITHIN :
        if((dir = rangeAny(start, end, e)) == 0) {
            if(rangeWithin(start, end, e) == 0) cnt = 1;
        }
        break;
    case GTF_MATCH_CONTAIN :
        if((dir = rangeAny(start, end, e)) == 0) {
            if(rangeContains(start, end, e) == 0) cnt = 1;
        }
        break;
    case GTF_MATCH_START :
        if((dir = rangeStart(start, end, e)) == 0) {
            cnt = 1;
        }
        break;
    case GTF_MATCH_END :
        if((dir = rangeEnd(start, end, e)) == 0) {
            cnt = 1;
        }
        break;
    default :
        if((dir = rangeAny(start, end, e)) == 0) {
            cnt = 1;
        }
        break;
    }

    if(cnt) {
        if(!matchingStrand(e, strand, strandType)) cnt = 0;
    }

    if(cnt && ffunc) {
        if(!ffunc(t, e)) cnt = 0;
    }

    if(max && cnt >= max) return max;

    if(direction) {
        if(dir > 0) return cnt;
        return cnt + countOverlapsEntry(t, e->right, start, end, strand, matchType, strandType, direction, max, ffunc);
    } else {
        if(dir < 0) return cnt;
        return cnt + countOverlapsEntry(t, e->left, start, end, strand, matchType, strandType, direction, max, ffunc);
    }
}

static void pushOverlapsNode(overlapSet *os, GTFtree *t, GTFnode *n, uint32_t start, uint32_t end, int matchType, FILTER_ENTRY_FUNC ffunc) {
    int dir;
    if(!n) return;
    dir = centerDirection(start, end, n);

    if(dir&1) {
        pushOverlaps(os, t, n->starts, start, end, matchType, 1, ffunc);
        pushOverlapsNode(os, t, n->left, start, end, matchType, ffunc);
    } 
    if(dir&2) {
        if(dir!=3) pushOverlaps(os, t, n->ends, start, end, matchType, 0, ffunc);
        pushOverlapsNode(os, t, n->right, start, end, matchType, ffunc);
    }
}

static int32_t countOverlapsNode(GTFtree *t, GTFnode *n, uint32_t start, uint32_t end, int strand, int matchType, int strandType, int32_t max, FILTER_ENTRY_FUNC ffunc) {
    int32_t cnt = 0;
    int dir;
    if(!n) return cnt;
    dir = centerDirection(start, end, n);

    if(dir&1) {
        cnt += countOverlapsEntry(t, n->starts, start, end, strand, matchType, strandType, 1, max, ffunc);
        if(max && cnt >= max) return max;
        cnt += countOverlapsNode(t, n->left, start, end, strand, matchType, strandType, max, ffunc);
        if(max && cnt >= max) return max;
    } 
    if(dir&2) {
        if(dir!=3) cnt += countOverlapsEntry(t, n->starts, start, end, strand, matchType, strandType, 0, max, ffunc);
        if(max && cnt >= max) return max;
        cnt += countOverlapsNode(t, n->right, start, end, strand, matchType, strandType, max, ffunc);
        if(max && cnt >= max) return max;
    }
    return cnt;
}

/*******************************************************************************
*
* Driver functions for end use.
*
*******************************************************************************/
overlapSet * findOverlaps(overlapSet *os, GTFtree *t, char *chrom, uint32_t start, uint32_t end, int strand, int matchType, int strandType, int keepOS, FILTER_ENTRY_FUNC ffunc) {
    int32_t tid = str2valHT(t->htChroms, chrom);
    overlapSet *out = os;

    if(out && !keepOS) os_reset(out);
    else if(!out) out = os_init(t);

    if(tid<0) return out;
    if(!t->balanced) {
        fprintf(stderr, "[findOverlaps] The tree has not been balanced! No overlaps will be returned.\n");
        return out;
    }

    pushOverlapsNode(out, t, (GTFnode*) t->chroms[tid]->tree, start, end, matchType, ffunc);
    if(out->l) filterStrand(out, strand, strandType);
    if(out->l) os_sort(out);

    return out;
}

int32_t countOverlaps(GTFtree *t, char *chrom, uint32_t start, uint32_t end, int strand, int matchType, int strandType, FILTER_ENTRY_FUNC ffunc) {
    int32_t tid = str2valHT(t->htChroms, chrom);
    if(tid<0) return 0;

    if(!t->balanced) {
        fprintf(stderr, "[countOverlaps] The tree has not been balanced! No overlaps will be returned.\n");
        return 0;
    }

    return countOverlapsNode(t, (GTFnode*) t->chroms[tid]->tree, start, end, strand, matchType, strandType, 0, ffunc);
}

int overlapsAny(GTFtree *t, char *chrom, uint32_t start, uint32_t end, int strand, int matchType, int strandType, FILTER_ENTRY_FUNC ffunc) {
    int32_t tid = str2valHT(t->htChroms, chrom);
    if(tid<0) return 0;

    if(!t->balanced) {
        fprintf(stderr, "[overlapsAny] The tree has not been balanced! No overlaps will be returned.\n");
        return 0;
    }

    return countOverlapsNode(t, (GTFnode*) t->chroms[tid]->tree, start, end, strand, matchType, strandType, 1, ffunc);
}
