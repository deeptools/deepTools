#include <Python.h>
#include <assert.h>
#include <inttypes.h>
#include "tree.h"
#include <assert.h>

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

//I don't know what happens if PyBytes_AsString(NULL) is used...
char *PyString_AsString(PyObject *obj) {
    return PyBytes_AsString(PyUnicode_AsASCIIString(obj));
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
    char *chrom = NULL, *name = NULL;
    uint32_t start, end, labelIdx;
    uint8_t strand;
    unsigned long lstrand, lstart, lend, llabelIdx;

    if(!(PyArg_ParseTuple(args, "skkskk", &chrom, &lstart, &lend, &name, &lstrand, &llabelIdx))) {
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

    //Actually add the entry
    if(addGTFentry(t, chrom, start, end, strand, name, labelIdx)) {
        PyErr_SetString(PyExc_RuntimeError, "pyAddEntry received an error while inserting an entry!");
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

static PyObject *pyFindOverlaps(pyGTFtree_t *self, PyObject *args) {
    GTFtree *t = self->t;
    char *chrom = NULL, *name = NULL, *transcript_id = NULL;
    int32_t i;
    uint32_t start, end;
    int strand = 3, strandType = 0, matchType = 0;
    unsigned long lstrand, lstart, lend, lmatchType, lstrandType, llabelIdx;
    overlapSet *os = NULL;
    PyObject *olist = NULL, *otuple = NULL;

    if(!(PyArg_ParseTuple(args, "skkkkks", &chrom, &lstart, &lend, &lstrand, &lmatchType, &lstrandType, &transcript_id))) {
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
        // Make the tuple
        otuple = PyTuple_New(4);
        if(!otuple) goto error;
        lstart = (unsigned long) os->overlaps[i]->start;
        lend = (unsigned long) os->overlaps[i]->end;
        name = getAttribute(t, os->overlaps[i], transcript_id);
        llabelIdx = (unsigned long) os->overlaps[i]->labelIdx;
        otuple = Py_BuildValue("(kksk)", lstart, lend, name, llabelIdx);
        if(!otuple) goto error;

        // Add the tuple
        if(PyList_SetItem(olist, i, otuple)) goto error;
        otuple = NULL;
    }

    return olist;

error:
    if(otuple) Py_DECREF(otuple);
    if(olist) Py_DECREF(olist);
    PyErr_SetString(PyExc_RuntimeError, "findOverlaps received an error!");
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
