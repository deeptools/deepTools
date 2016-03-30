#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "murmur3.h"
#include "gtf.h"

uint64_t hashString(char *s) {
    int len = strlen(s);
    uint64_t hash_val[2];
    uint32_t seed = 0xAAAAAAAA;
#if UINTPTR_MAX == 0xffffffff
    MurmurHash3_x86_128((void *) s, len, seed, (void *) &hash_val);
#else
    MurmurHash3_x64_128((void *) s, len, seed, (void *) &hash_val);
#endif
    return hash_val[0];
}

hashTable *initHT(uint64_t size) {
    hashTable *ht = calloc(1, sizeof(hashTable));
    assert(ht);

    ht->elements = calloc(size, sizeof(hashTableElement*));
    assert(ht->elements);
    ht->str = calloc(size, sizeof(char*));
    assert(ht->str);

    ht->m = size;
    return ht;
}

void insertHTelement(hashTable *ht, hashTableElement *e, uint64_t hash) {
    uint64_t i = hash%ht->m;
    hashTableElement *curr = ht->elements[i];
    if(!curr) ht->elements[i] = e;
    else {
        while(curr->next) curr = curr->next;
        curr->next = e;
    }
}

static void rehashElement(hashTable *ht, hashTableElement *e) {
    hashTableElement *next = e->next;
    if(!e) return;
    uint64_t hash = hashString(ht->str[e->val]);
    e->next = NULL;
    insertHTelement(ht, e, hash);
    if(next) rehashElement(ht, next);
}

static void rehashHT(hashTable *ht) {
    int32_t i;
    hashTableElement *e;
    for(i=0; i<ht->l; i++) {
        if(ht->elements[i]) {
            e = ht->elements[i];
            ht->elements[i] = NULL;
            rehashElement(ht, e);
        }
    }
}

static void growHT(hashTable *ht) {
    int i;
    ht->m = ht->l+1;
    kroundup32(ht->m);
    ht->str = realloc(ht->str, ht->m*sizeof(char*));
    assert(ht->str);
    ht->elements = realloc(ht->elements, ht->m*sizeof(hashTableElement*));

    for(i=ht->l; i<ht->m; i++) {
        ht->str[i] = NULL;
        ht->elements[i] = NULL;
    }
    rehashHT(ht);
}

//don't do this if the element's already in the table!
int32_t addHTelement(hashTable *ht, char *s) {
    if(!s) return -1;
    uint64_t hash = hashString(s);
    int32_t val = ht->l++;
    if(ht->l >= ht->m) growHT(ht);
    ht->str[val] = strdup(s);

    hashTableElement *e = calloc(1, sizeof(hashTableElement));
    assert(e);
    e->val = val;
    insertHTelement(ht, e, hash);
    return val;
}

void destroyHTelement(hashTableElement *e) {
    hashTableElement *next = e->next;
    free(e);
    if(next) destroyHTelement(next);
}

void destroyHT(hashTable *ht) {
    int i;

    for(i=0; i<ht->l; i++) free(ht->str[i]);

    for(i=0; i<ht->m; i++) {
        if(ht->elements[i]) destroyHTelement(ht->elements[i]);
    }
    free(ht->elements);
    free(ht->str);
    free(ht);
}

int strExistsHT(hashTable *ht, char *s) {
    if(!s) return 0;
    uint64_t h = hashString(s);
    hashTableElement *curr = ht->elements[h%ht->m];
    while(curr) {
        if(strcmp(ht->str[curr->val], s) == 0) return 1;
        curr = curr->next;
    }
    return 0;
}

//Returns -1 if not present
int32_t str2valHT(hashTable *ht, char *s) {
    if(!s) return -1;
    uint64_t h = hashString(s);
    hashTableElement *curr = ht->elements[h%ht->m];
    while(curr) {
        if(strcmp(ht->str[curr->val], s) == 0) return curr->val;
        curr = curr->next;
    }
    return -1;
}

//Returns NULL on error
char *val2strHT(hashTable *ht, int32_t val) {
    if(val<0) return NULL;
    if(val>=ht->l) return NULL;
    return ht->str[val];
}

int hasAttribute(GTFtree *t, GTFentry *e, char *str) {
    int32_t i, key = str2valHT(t->htAttributes, str);

    for(i=0; i<e->nAttributes; i++) {
        if(e->attrib[i]->key == key) return 1;
    }
    return 0;
}

//Returns NULL if the entry lacks the attribute
char *getAttribute(GTFtree *t, GTFentry *e, char *str) {
    int32_t i, key = str2valHT(t->htAttributes, str);

    for(i=0; i<e->nAttributes; i++) {
        if(e->attrib[i]->key == key) return val2strHT(t->htAttributes, e->attrib[i]->val);
    }
    return NULL;
}
