#include <Python.h>
#include <assert.h>
#include <inttypes.h>
#include "tree.h"
#include <assert.h>
#include <float.h>

static void pyGTFDealloc(pyGTFtree_t *self) {
    if(self->t) destroyGTFtree(self->t);
    PyObject_DEL(self);
}

#if PY_MAJOR_VERSION >= 3
//Return 1 iff obj is a ready unicode type
int PyString_Check(PyObject *obj) {
    if(PyUnicode_Check(obj)) {
        return PyUnicode_READY(obj)+1;
    }
    return 0;
}

char *PyString_AsString(PyObject *obj) {
    return PyBytes_AsString(PyUnicode_AsASCIIString(obj));
}

PyObject *PyString_FromString(char *s) {
    return PyUnicode_FromString(s);
}
#endif

//Will return 1 for long or int types currently
int isNumeric(PyObject *obj) {
#if PY_MAJOR_VERSION < 3
    if(PyInt_Check(obj)) return 1;
#endif
    return PyLong_Check(obj);
}

//On error, throws a runtime error, so use PyErr_Occurred() after this
uint32_t Numeric2Uint(PyObject *obj) {
    long l;
#if PY_MAJOR_VERSION < 3
    if(PyInt_Check(obj)) {
        return (uint32_t) PyInt_AsLong(obj);
    }
#endif
    l = PyLong_AsLong(obj);
    //Check bounds
    if(l > 0xFFFFFFFF) {
        PyErr_SetString(PyExc_RuntimeError, "Length out of bounds for a bigWig file!");
        return (uint32_t) -1;
    }
    return (uint32_t) l;
}

static PyObject *pyGTFinit(PyObject *self, PyObject *args) {
    GTFtree *t = NULL;
    pyGTFtree_t *pt;

    t = initGTFtree();
    if(!t) return NULL;

    pt = PyObject_New(pyGTFtree_t, &pyGTFtree);
    if(!pt) goto error;

    pt->t = t;
    return (PyObject*) pt;

error:
    if(t) destroyGTFtree(t);
    PyErr_SetString(PyExc_RuntimeError, "Received an error during tree initialization!");
    return NULL;
}

static PyObject *pyAddEntry(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    char *chrom = NULL, *name = NULL, *sscore = NULL;
    uint32_t start, end, labelIdx;
    double score;
    uint8_t strand;
    unsigned long lstrand, lstart, lend, llabelIdx;

    if(!(PyArg_ParseTuple(args, "skkskks", &chrom, &lstart, &lend, &name, &lstrand, &llabelIdx, &sscore))) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEntry received an invalid or missing argument!");
        return NULL;
    }

    //Convert all of the longs
    if(lstart >= (uint32_t) -1 || lend >= (uint32_t) -1 || lend <= lstart) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEntry received invalid bounds!");
        return NULL;
    }
    start = (uint32_t) lstart;
    end = (uint32_t) lend;
    if(lstrand != 0 && lstrand != 1 && lstrand != 3) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEntry received an invalid strand!");
        return NULL;
    }
    strand = (uint8_t) lstrand;
    if(llabelIdx >= (uint32_t) -1) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEntry received an invalid label idx (too large)!");
        return NULL;
    }
    labelIdx = (uint32_t) llabelIdx;

    //Handle the score
    if(strcmp(sscore, ".") == 0) {
        score = DBL_MAX;
    } else {
        score = strtod(sscore, NULL);
    }

    //Actually add the entry
    if(addGTFentry(t, chrom, start, end, strand, name, labelIdx, score)) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEntry received an error while inserting an entry!");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *pyAddEnrichmentEntry(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    char *chrom = NULL, *sscore = NULL, *feature = NULL;
    uint32_t start, end;
    double score;
    uint8_t strand;
    unsigned long lstrand, lstart, lend;

    if(!(PyArg_ParseTuple(args, "skkkss", &chrom, &lstart, &lend, &lstrand, &sscore, &feature))) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEnrichmentEntry received an invalid or missing argument!");
        return NULL;
    }

    //Convert all of the longs
    if(lstart >= (uint32_t) -1 || lend >= (uint32_t) -1 || lend <= lstart) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEnrichmentEntry received invalid bounds!");
        return NULL;
    }
    start = (uint32_t) lstart;
    end = (uint32_t) lend;
    if(lstrand != 0 && lstrand != 1 && lstrand != 3) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEnrichmentEntry received an invalid strand!");
        return NULL;
    }
    strand = (uint8_t) lstrand;

    //Handle the score
    if(strcmp(sscore, ".") == 0) {
        score = DBL_MAX;
    } else {
        score = strtod(sscore, NULL);
    }

    //Actually add the entry
    if(addEnrichmententry(t, chrom, start, end, strand, score, feature)) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEnrichmentEntry received an error while inserting an entry!");
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *pyVine2Tree(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    sortGTF(t);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *pyPrintGTFtree(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    printGTFtree(t);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *pyCountEntries(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    uint32_t nEntries = 0;
    unsigned long lnEntries = 0;
    int32_t i;
    PyObject *out = NULL;

    for(i=0; i<t->n_targets; i++) {
        nEntries += t->chroms[i]->n_entries;
    }
    lnEntries = (unsigned long) nEntries;
    out = PyLong_FromUnsignedLong(lnEntries);

    return out;
}

static PyObject *pyIsTree(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    if(t->balanced) Py_RETURN_TRUE;
    Py_RETURN_FALSE;
}

static PyObject *pyHasOverlaps(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    int rv;
    uint32_t minDistance = (uint32_t) -1;
    unsigned long long ominDistance;
    PyObject *otuple = NULL, *oval = NULL;

    rv = hasOverlaps(t, &minDistance);
    ominDistance = minDistance; // ominDistance should have at least as much space as minDistance
    otuple = PyTuple_New(2);
    if(!otuple) {
        PyErr_SetString(PyExc_RuntimeError, "Could not allocate space for a tuple!\n");
        return NULL;
    }
    oval = PyLong_FromUnsignedLongLong(ominDistance);
    if(!oval) {
        PyErr_SetString(PyExc_RuntimeError, "Could not allocate space for a single integer!\n");
        return NULL;
    }

    if(rv) {
        Py_INCREF(Py_True);
        PyTuple_SET_ITEM(otuple, 0, Py_True);
    } else {
        Py_INCREF(Py_False);
        PyTuple_SET_ITEM(otuple, 0, Py_False);
    }
    PyTuple_SetItem(otuple, 1, oval);

    return otuple;
}

static PyObject *pyFindOverlaps(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    char *chrom = NULL, *name = NULL, *transcript_id = NULL, strandChar;
    int32_t i;
    uint32_t start, end;
    int strand = 3, strandType = 0, matchType = 0;
    unsigned long lstrand, lstart, lend, lmatchType, lstrandType, llabelIdx;
    overlapSet *os = NULL;
    PyObject *olist = NULL, *otuple = NULL, *includeStrand = Py_False, *oscore = NULL;

    if(!(PyArg_ParseTuple(args, "skkkkksO", &chrom, &lstart, &lend, &lstrand, &lmatchType, &lstrandType, &transcript_id, &includeStrand))) {
        PyErr_SetString(PyExc_RuntimeError, "pyFindOverlaps received an invalid or missing argument!");
        return NULL;
    }

    //I'm assuming that this is never called outside of the module
    strandType = (int) lstrandType;
    strand = (int) lstrand;
    matchType = (int) matchType;
    start = (uint32_t) lstart;
    end = (uint32_t) lend;

    os = findOverlaps(NULL, t, chrom, start, end, strand, matchType, strandType, 0, NULL);

    // Did we receive an error?
    if(!os) {
        PyErr_SetString(PyExc_RuntimeError, "findOverlaps returned NULL!");
        return NULL;
    }

    // Convert the overlapSet to a list of tuples
    olist = PyList_New(os->l);
    if(!olist) goto error;
    for(i=0; i<os->l; i++) {
        // Make the tuple
        if(includeStrand == Py_True) {
            otuple = PyTuple_New(6);
        } else {
            otuple = PyTuple_New(5);
        }
        if(!otuple) goto error;
        lstart = (unsigned long) os->overlaps[i]->start;
        lend = (unsigned long) os->overlaps[i]->end;
        name = getAttribute(t, os->overlaps[i], transcript_id);
        llabelIdx = (unsigned long) os->overlaps[i]->labelIdx;
        strandChar = '.';
        if(os->overlaps[i]->strand == 0) {
            strandChar = '+';
        } else if(os->overlaps[i]->strand == 1) {
            strandChar = '-';
        }
        if (os->overlaps[i]->score == DBL_MAX) {
            oscore = Py_BuildValue("s", ".");
        } else {
            oscore = Py_BuildValue("d", os->overlaps[i]->score);
        }
        if(!oscore) goto error;
        if(includeStrand == Py_True) {
            otuple = Py_BuildValue("(kkskcO)", lstart, lend, name, llabelIdx, strandChar, oscore);
        } else {
            otuple = Py_BuildValue("(kkskO)", lstart, lend, name, llabelIdx, oscore);
        }
        if(!otuple) goto error;

        // Add the tuple
        if(PyList_SetItem(olist, i, otuple)) goto error;
        otuple = NULL;
    }
    os_destroy(os);

    return olist;

error:
    if(otuple) Py_DECREF(otuple);
    if(olist) Py_DECREF(olist);
    PyErr_SetString(PyExc_RuntimeError, "findOverlaps received an error!");
    return NULL;
}

static PyObject *pyFindOverlappingFeatures(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    char *chrom = NULL;
    int32_t i;
    uint32_t start, end;
    int strand = 3, strandType = 0, matchType = 0;
    unsigned long lstrand, lstart, lend, lmatchType, lstrandType;
    overlapSet *os = NULL;
    PyObject *olist = NULL, *ostring = NULL;

    if(!(PyArg_ParseTuple(args, "skkkkk", &chrom, &lstart, &lend, &lstrand, &lmatchType, &lstrandType))) {
        PyErr_SetString(PyExc_RuntimeError, "pyFindOverlaps received an invalid or missing argument!");
        return NULL;
    }

    //I'm assuming that this is never called outside of the module
    strandType = (int) lstrandType;
    strand = (int) lstrand;
    matchType = (int) matchType;
    start = (uint32_t) lstart;
    end = (uint32_t) lend;

    os = findOverlaps(NULL, t, chrom, start, end, strand, matchType, strandType, 0, NULL);

    // Did we receive an error?
    if(!os) {
        PyErr_SetString(PyExc_RuntimeError, "findOverlaps returned NULL!");
        return NULL;
    }

    if(!os->l) {
        os_destroy(os);
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Convert the overlapSet to a list of tuples
    olist = PyList_New(os->l);
    if(!olist) goto error;
    for(i=0; i<os->l; i++) {
        //Make the python string
        ostring = PyString_FromString(val2strHT(t->htFeatures, os->overlaps[i]->feature));
        if(!ostring) goto error;

        // Add the item
        if(PyList_SetItem(olist, i, ostring)) goto error;
        ostring = NULL;
    }
    os_destroy(os);

    return olist;

error:
    if(ostring) Py_DECREF(ostring);
    if(olist) Py_DECREF(olist);
    PyErr_SetString(PyExc_RuntimeError, "findOverlappingFeatures received an error!");
    return NULL;
}

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_tree(void) {
    PyObject *res;
    errno = 0;

    if(PyType_Ready(&pyGTFtree) < 0) return NULL;
    res = PyModule_Create(&treemodule);
    if(!res) return NULL;

    Py_INCREF(&pyGTFtree);
    PyModule_AddObject(res, "pyGTFtree", (PyObject *) &pyGTFtree);

    return res;
}
#else
//Python2 initialization
PyMODINIT_FUNC inittree(void) {
    errno = 0; //Sometimes libpython2.7.so is missing some links...
    if(PyType_Ready(&pyGTFtree) < 0) return;
    Py_InitModule3("tree", treeMethods, "A module for handling GTF files for deepTools");
}
#endif
