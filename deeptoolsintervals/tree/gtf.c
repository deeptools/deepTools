#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include "gtf.h"

//Nodes for the interval tree
typedef struct {
    int32_t center;
    int32_t l;
    GTFentry *start;
    GTFentry *end;
    struct treeNode *left;
    struct treeNode *right;
} treeNode;

//The sizes shouldn't be preset...
GTFtree * initGTFtree() {
    GTFtree *t = calloc(1, sizeof(GTFtree));
    assert(t);

    //Initialize the hash tables
    t->htChroms = initHT(128);
    t->htSources = initHT(128);
    t->htFeatures = initHT(128);
    t->htAttributes = initHT(128);

    return t;
}

void destroyGTFentry(GTFentry *e) {
    int32_t i;
    if(!e) return;
    if(e->right) destroyGTFentry(e->right);
    for(i=0; i<e->nAttributes; i++) {
        if(e->attrib[i]) free(e->attrib[i]);
    }
    if(e->attrib) free(e->attrib);
    free(e);
}

void destroyGTFnode(GTFnode *n) {
    if(n->left) destroyGTFnode(n->left);
    if(n->starts) destroyGTFentry(n->starts);
    if(n->right) destroyGTFnode(n->right);
    free(n);
}

void destroyGTFchrom(GTFchrom *c, int balanced) {
    if(balanced) destroyGTFnode((GTFnode*) c->tree);
    else destroyGTFentry((GTFentry*) c->tree);
    free(c);
}

//This need to handle htTargets and htGenes still
void destroyGTFtree(GTFtree *t) {
    uint32_t i;
    for(i=0; i<t->n_targets; i++) {
        destroyGTFchrom(t->chroms[i], t->balanced);
    }

    destroyHT(t->htChroms);
    destroyHT(t->htSources);
    destroyHT(t->htFeatures);
    destroyHT(t->htAttributes);

    free(t->chroms);
    free(t);
}

void addChrom(GTFtree *t) {
    int i;

    t->n_targets++;
    //Grow if needed
    if(t->n_targets >= t->m) {
        t->m++;
        kroundup32(t->m);
        t->chroms = realloc(t->chroms, t->m * sizeof(GTFchrom *));
        assert(t->chroms);
        for(i=t->n_targets-1; i<t->m; i++) t->chroms[i] = NULL;
    }

    assert(!t->chroms[t->n_targets-1]); //We shouldn't be adding over anything...
    t->chroms[t->n_targets-1] = calloc(1,sizeof(GTFchrom));
    assert(t->chroms[t->n_targets-1]);
}


//Returns NULL on error
static Attribute *makeAttribute(GTFtree *t, char *value) {
    int32_t idx;
    Attribute *a = malloc(sizeof(Attribute));
    if(!a) return NULL;

    if(!strExistsHT(t->htAttributes, "transcript_id")) {
        idx = addHTelement(t->htAttributes, "transcript_id");
    } else {
        idx = str2valHT(t->htAttributes, "transcript_id");
    }
    a->key = idx;
    if(!strExistsHT(t->htAttributes, value)) {
        idx = addHTelement(t->htAttributes, value);
    } else {
        idx = str2valHT(t->htAttributes, value);
    }
    a->val = idx;

    return a;
}

/* This currently hard-codes the following:
    feature
    source
    frame
    all attributes (most are skipped)

    returns 1 on error
*/
int addGTFentry(GTFtree *t, char *chrom, uint32_t start, uint32_t end, uint8_t strand, char *transcriptID, uint32_t labelIDX, double score) {
    int32_t IDchrom, IDfeature, IDsource;
    char feature[] = "transcript", source[] = "deepTools";
    uint8_t frame = 3;
    GTFentry *e = NULL;
    Attribute *a = NULL;
    Attribute **attribs = calloc(1, sizeof(Attribute *));
    if(!attribs) return 1;

    //Get the chromosome ID
    if(!strExistsHT(t->htChroms, chrom)) {
        addChrom(t);
        IDchrom = addHTelement(t->htChroms, chrom);
    } else {
        IDchrom = str2valHT(t->htChroms, chrom);
    }

    //Handle the hard-coded stuff, which in case they're ever requested
    if(!strExistsHT(t->htSources, source)) {
        IDsource = addHTelement(t->htSources, source);
    } else {
        IDsource = str2valHT(t->htSources, source);
    }
    if(!strExistsHT(t->htFeatures, feature)) {
        IDfeature = addHTelement(t->htFeatures, feature);
    } else {
        IDfeature = str2valHT(t->htFeatures, feature);
    }

    //Create the attribute
    a = makeAttribute(t, transcriptID);
    if(!a) goto error;
    attribs[0] = a;

    //Initialize the entry
    e = malloc(sizeof(GTFentry));
    if(!e) goto error;
    e->right = NULL;

    e->chrom = IDchrom;
    e->feature = IDfeature;
    e->source = IDsource;
    e->start = start;
    e->end = end;
    e->strand = strand;
    e->frame = frame;
    e->score = score;
    e->nAttributes = 1;
    e->attrib = attribs;
    e->labelIdx = labelIDX;

    if(t->chroms[IDchrom]->tree) {
        e->left = ((GTFentry*) t->chroms[IDchrom]->tree)->left;
        e->left->right = e;
        ((GTFentry*) t->chroms[IDchrom]->tree)->left = e;
    } else {
        t->chroms[IDchrom]->tree = (void *) e;
        e->left = e;
    }
    t->chroms[IDchrom]->n_entries++;

    return 0;

error:
    if(attribs) free(attribs);
    if(a) free(a);
    if(e) free(e);
    return 1;
}

/* This currently hard-codes the following:
    name
    source
    frame
    all attributes

    returns 1 on error
*/
int addEnrichmententry(GTFtree *t, char *chrom, uint32_t start, uint32_t end, uint8_t strand, double score, char *feature) {
    int32_t IDchrom, IDfeature, IDsource;
    char source[] = "deepTools";
    uint8_t frame = 3;
    GTFentry *e = NULL;
    //Attribute **attribs = calloc(1, sizeof(Attribute *));
    //if(!attribs) return 1;

    //Get the chromosome ID
    if(!strExistsHT(t->htChroms, chrom)) {
        addChrom(t);
        IDchrom = addHTelement(t->htChroms, chrom);
    } else {
        IDchrom = str2valHT(t->htChroms, chrom);
    }

    //Handle the hard-coded stuff, in case they're ever requested
    if(!strExistsHT(t->htSources, source)) {
        IDsource = addHTelement(t->htSources, source);
    } else {
        IDsource = str2valHT(t->htSources, source);
    }

    if(!strExistsHT(t->htFeatures, feature)) {
        IDfeature = addHTelement(t->htFeatures, feature);
    } else {
        IDfeature = str2valHT(t->htFeatures, feature);
    }

    //Initialize the entry
    e = malloc(sizeof(GTFentry));
    if(!e) goto error;
    e->right = NULL;

    e->chrom = IDchrom;
    e->feature = IDfeature;
    e->source = IDsource;
    e->start = start;
    e->end = end;
    e->strand = strand;
    e->frame = frame;
    e->score = score;
    e->nAttributes = 0;
    e->attrib = NULL;

    if(t->chroms[IDchrom]->tree) {
        e->left = ((GTFentry*) t->chroms[IDchrom]->tree)->left;
        e->left->right = e;
        ((GTFentry*) t->chroms[IDchrom]->tree)->left = e;
    } else {
        t->chroms[IDchrom]->tree = (void *) e;
        e->left = e;
    }
    t->chroms[IDchrom]->n_entries++;

    return 0;

error:
    //if(attribs) free(attribs);
    //if(a) free(a);
    if(e) free(e);
    return 1;
}

/*******************************************************************************
*
* Sorting functions
*
*******************************************************************************/
GTFentry *getMiddleR(GTFentry *e, uint32_t pos) {
    uint32_t i;
    GTFentry *tmp, *o = e;

    if(!o->right) return o;
    for(i=1; i<pos; i++) {
        assert(o->right);
        o = o->right;
    }
    tmp = o;
    assert(o->right);
    o = o->right;
    tmp->right = NULL;
    return o;
}
GTFentry *getMiddleL(GTFentry *e, uint32_t pos) {
    uint32_t i;
    GTFentry *tmp, *o = e;

    if(!o->left) {
        return o;
    }
    for(i=1; i<pos; i++) {
        assert(o->left);
        o = o->left;
    }
    tmp = o;
    assert(o->left);
    o = o->left;
    tmp->left = NULL;
    return o;
}

int cmpRangesStart(GTFentry *a, GTFentry *b) {
    if(!b && !a) return 0;
    if(!b) return -1;
    if(!a) return 1;
    if(a->start < b->start) return -1;
    if(b->start < a->start) return 1;
    if(b->end < a->end) return 1;
    return -1;
}
int cmpRangesEnd(GTFentry *a, GTFentry *b) {
    if(!b && !a) return 0;
    if(!a) return 1;
    if(!b) return -1;
    if(a->end > b->end) return -1;
    if(b->end > a->end) return 1;
    if(a->start > b->start) return -1;
    return 1;
}

GTFentry *mergeSortStart(GTFentry *a, GTFentry *b) {
    GTFentry *o = a, *last;
    int i = cmpRangesStart(a,b);

    if(i<0) {
        o = a;
        a = a->right;
    } else if(i>0) {
        o = b;
        b = b->right;
    } else{
         return NULL;
    }
    last = o;
    last->right = NULL;

    while((i=cmpRangesStart(a,b))) {
        if(i>0) {
            last->right= b;
            last = b;
            b = b->right;
        } else {
            last->right= a;
            last = a;
            a = a->right;
        }
    }
    last->right = NULL;
    return o;
}

GTFentry *mergeSortEnd(GTFentry *a, GTFentry *b) {
    GTFentry *o = a, *last;
    int i = cmpRangesEnd(a,b);

    if(i<0) {
        o = a;
        a = a->left;
    } else if(i>0) {
        o = b;
        b = b->left;
    } else {
        return NULL;
    }
    last = o;
    last->left = NULL;

    while((i=cmpRangesEnd(a,b))) {
        if(i<0) {
            assert(a != last);
            last->left = a;
            last = a;
            a = a->left;
        } else {
            assert(b != last);
            last->left = b;
            last = b;
            b = b->left;
        }
    }
    last->left = NULL;
    return o;
}

GTFentry *sortTreeStart(GTFentry *e, uint32_t l) {
    if(l==1) return e;
    uint32_t half = l/2;
    GTFentry *middle = getMiddleR(e, half);

    return mergeSortStart(sortTreeStart(e,half), sortTreeStart(middle,half+(l&1)));
}

GTFentry *sortTreeEnd(GTFentry *e, uint32_t l) {
    if(l==1) {
        e->left = NULL; //The list is circular, so...
        return e;
    }
    uint32_t half = l/2;
    assert(e->left);
    assert(e != e->left);
    GTFentry *middle = getMiddleL(e, half);
    assert(e != middle);
    assert(e != e->left);

    return mergeSortEnd(sortTreeEnd(e,half), sortTreeEnd(middle,half+(l&1)));
}

/*******************************************************************************
*
* Functions for interval tree construction
*
*******************************************************************************/

//Note the returned object is the rightmost interval sorted by end position
GTFentry *sortChrom(GTFchrom *c) {
    GTFentry *e = ((GTFentry *)c->tree)->left;
    ((GTFentry*) c->tree)->left = NULL;
    c->tree = (void *) sortTreeStart((GTFentry *) c->tree, c->n_entries);
    e = sortTreeEnd(e, c->n_entries);
    return e;
}

uint32_t getCenter(GTFentry *ends) {
    GTFentry *slow = ends;
    GTFentry *fast = ends;

    while(fast->left && fast->left->left) {
        slow = slow->left;
        fast = fast->left->left;
    }

    return slow->end-1;
}

GTFentry *getMembers(GTFentry **members, GTFentry **rStarts, GTFentry *starts, uint32_t pos) {
    GTFentry *tmp, *newStarts = NULL;
    GTFentry *last = NULL, *lastMember = NULL;

    *members = NULL, *rStarts = NULL;

    while(starts && starts->start <= pos) {
        if(starts->end > pos) {
            tmp = starts->right;
            if(!*members) {
                lastMember = starts;
                *members = starts;
            } else {
                lastMember->right = starts;
                lastMember = starts;
            }
            starts->right = NULL;
            starts = tmp;
        } else {
            if(!newStarts) {
                newStarts = starts;
                last = starts;
            } else {
                last->right = starts;
                last = starts;
            }
            starts = starts->right;
        }
    }
    *rStarts = starts;
    if(lastMember) lastMember->right = NULL;
    if(last) last->right = NULL;
    assert(*members);
    return newStarts;
}
GTFentry *getRMembers(GTFentry **members, GTFentry **lEnds, GTFentry *ends, uint32_t pos) {
    GTFentry *tmp, *newEnds = NULL;
    GTFentry *last = NULL, *lastMember = NULL;

    *members = NULL, *lEnds = NULL;

    while(ends && ends->end > pos) {
        tmp = ends->left;
        if(ends->start <= pos) {
            if(!*members) {
                *members = ends;
                lastMember = ends;
            } else {
                lastMember->left = ends;
                lastMember = ends;
            }
        } else {
            if(!newEnds) {
                newEnds = ends;
                last = ends;
            } else {
                last->left = ends;
                last = ends;
            }
        }
        ends->left = NULL;
        ends = tmp;
    }
    *lEnds = ends;
    assert(*members);
    lastMember->left = NULL;
    if(newEnds) last->left = NULL;
    return newEnds;
}

GTFnode *makeIntervalTree(GTFentry *starts, GTFentry *ends) {
    uint32_t center = getCenter(ends);//, nMembers;
    GTFentry *rStarts = NULL; //getRStarts(starts, center);
    GTFentry *lEnds = NULL; //getLEnds(ends, center);
    GTFentry *memberStarts = NULL, *memberEnds = NULL;
    GTFnode *out = calloc(1, sizeof(GTFnode));
    assert(out);

    starts = getMembers(&memberStarts, &rStarts, starts, center);
    ends = getRMembers(&memberEnds, &lEnds, ends, center);

    out->center = center;
    out->starts = memberStarts;
    out->ends = memberEnds;
    if(lEnds && starts) {
        out->left = makeIntervalTree(starts, lEnds);
    } else {
        out->left = NULL;
    }
    if(rStarts && ends) {
        out->right = makeIntervalTree(rStarts, ends);
    } else {
        out->right = NULL;
    }

    return out;
}

void sortGTF(GTFtree *t) {
    int32_t i;
    GTFentry *ends;

    for(i=0; i<t->n_targets; i++) {
        ends = sortChrom(t->chroms[i]);
        t->chroms[i]->tree = (void*) makeIntervalTree((GTFentry*) t->chroms[i]->tree, ends);
    }
    t->balanced = 1;
}

int nodeHasOverlaps(GTFnode *node, int firstNode, uint32_t *lpos, uint32_t *minDistance) {
    int rv = 0;
    GTFentry *e = node->starts;

    // Go down the left
    if(node->left) {
        rv = nodeHasOverlaps(node->left, firstNode, lpos, minDistance);
        if(rv) return rv;
    } else if(firstNode) {
        //This only has to be specially set on the left-most node
        *lpos = e->end;
        *minDistance = e->start;
        e = e->right;
    }

    // Test this node
    while(e) {
        if(e->start < *lpos) {
            *minDistance = 0;
            return 1;
        }
        if(e->start - *lpos < *minDistance) *minDistance = e->start - *lpos;
        *lpos = e->end;
        e = e->right;
    }

    if(node->right) return nodeHasOverlaps(node->right, 0, lpos, minDistance);
    return rv;
}

int hasOverlapsChrom(GTFchrom *chrom, uint32_t *minDistance) {
    uint32_t lpos;
    if(chrom->n_entries < 2) return 0;
    return nodeHasOverlaps((GTFnode*) chrom->tree, 1, &lpos, minDistance);
}

// Given a GTF tree, returning 1 if ANY of the entries overlap with each other, 0 otherwise
// minDistance is updated to return the minimum distance between intervals. This will be 0 if there are overlaps.
int hasOverlaps(GTFtree *t, uint32_t *minDistance) {
    int32_t i;
    int rv = 0;
    *minDistance = (uint32_t) -1;

    for(i=0; i<t->n_targets; i++) {
        rv = hasOverlapsChrom(t->chroms[i], minDistance);
        if(rv) return rv;
    }
    return rv;
}

/*******************************************************************************
*
* Misc. functions
*
*******************************************************************************/
void printBalancedGTF(GTFnode *n, const char *chrom) {
    kstring_t ks, ks2;
    ks.s = NULL; ks.l = ks.m = 0;
    ks2.s = NULL; ks2.l = ks2.m = 0;
    kputs(chrom, &ks);
    kputc(':', &ks);
    kputuw(n->center, &ks);
    if(n->left) {
        kputs(chrom, &ks2);
        kputc(':', &ks2);
        kputuw(n->left->center, &ks2);
        printf("\t\"%s\" -> \"%s\";\n", ks.s, ks2.s);
        printBalancedGTF(n->left, chrom);
    }

    printf("\t\"%s:%"PRIu32"\" [shape=box];\n", chrom, n->center);
    GTFentry *e = n->starts;
    if(e) printGTFvineStart(e, chrom, ks.s);
    if(n->ends) printGTFvineStartR(n->ends, chrom, ks.s);

    if(n->right) {
        ks2.l = 0;
        kputs(chrom, &ks2);
        kputc(':', &ks2);
        kputuw(n->right->center, &ks2);
        printf("\t\"%s\" -> \"%s\";\n", ks.s, ks2.s);
        printBalancedGTF(n->right, chrom);
    }
    free(ks.s);
    if(ks2.s) free(ks2.s);
}

void printGTFvineR(GTFentry *e, const char* chrom) {
    if(e->left == e) return;
    if(!e->left) return;
    printf("\t\"%s:%"PRIu32"-%"PRIu32"\" -> \"%s:%"PRIu32"-%"PRIu32"\" [color=red];\n", chrom, e->start, e->end, chrom, e->left->start, e->left->end);
    printGTFvineR(e->left, chrom);
}
void printGTFvineStartR(GTFentry *e, const char *chrom, const char *str) {
    printf("\t\"%s\" -> \"%s:%"PRIu32"-%"PRIu32"\" [color=red];\n", str, chrom, e->start, e->end);
    if(e->left) printGTFvineR(e, chrom);
}

void printGTFvine(GTFentry *e, const char* chrom) {
    if(!e->right) return;
    printf("\t\"%s:%"PRIu32"-%"PRIu32"\" -> \"%s:%"PRIu32"-%"PRIu32"\";\n", chrom, e->start, e->end, chrom, e->right->start, e->right->end);
    printGTFvine(e->right, chrom);
}
void printGTFvineStart(GTFentry *e, const char *chrom, const char *str) {
    printf("\t\"%s\" -> \"%s:%"PRIu32"-%"PRIu32"\";\n", str, chrom, e->start, e->end);
    if(e->right) printGTFvine(e, chrom);
}

void printGTFtree(GTFtree *t) {
    int32_t i;
    const char *chromName;

    if(t->balanced) printf("digraph balancedTree {\n");
    else printf("digraph unbalancedTree {\n");

    for(i=0; i<t->n_targets; i++) {
        chromName = val2strHT(t->htChroms, i);
        if(t->balanced) {
            printBalancedGTF((GTFnode*) t->chroms[i]->tree, chromName);
        } else {
            printGTFvineStart((GTFentry*) t->chroms[i]->tree, chromName, chromName);
        }
    }
    printf("}\n");
}
