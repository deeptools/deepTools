#include <Python.h>
#include <structmember.h>
#include "gtf.h"

typedef struct {
    PyObject_HEAD
    GTFtree *t;
} pyGTFtree_t;

/* 
  Remove all asserts and ensure that the new return values are honoured.

  Profile the code with a test to eliminate unneeded cruft.
*/

static PyObject *pyGTFinit(PyObject *self, PyObject *args);
static PyObject *pyAddEntry(pyGTFtree_t *self, PyObject *args);
static PyObject *pyAddEnrichmentEntry(pyGTFtree_t *self, PyObject *args);
static PyObject *pyVine2Tree(pyGTFtree_t *self, PyObject *args);
static PyObject *pyPrintGTFtree(pyGTFtree_t *self, PyObject *args);
static PyObject *pyCountEntries(pyGTFtree_t *self, PyObject *args);
static PyObject *pyFindOverlaps(pyGTFtree_t *self, PyObject *args);
static PyObject *pyFindOverlappingFeatures(pyGTFtree_t *self, PyObject *args);
static PyObject *pyIsTree(pyGTFtree_t *self, PyObject *args);
static PyObject *pyHasOverlaps(pyGTFtree_t *self, PyObject *args);
static void pyGTFDealloc(pyGTFtree_t *self);

static PyMethodDef treeMethods[] = {
    {"initTree", (PyCFunction) pyGTFinit, METH_VARARGS,
"Initialize the tree\n"},
    {"addEntry", (PyCFunction) pyAddEntry, METH_VARARGS,
"Some documentation for pyAddEntry\n"},
    {"addEnrichmentEntry", (PyCFunction) pyAddEnrichmentEntry, METH_VARARGS,
"Some documentation for pyAddEnrichmentEntry\n"},
    {"finish", (PyCFunction) pyVine2Tree, METH_VARARGS,
"This must be called after ALL entries from ALL files have been added.\n"},
    {"printGTFtree", (PyCFunction) pyPrintGTFtree, METH_VARARGS,
"Prints a text representation in dot format.\n"},
    {"countEntries", (PyCFunction) pyCountEntries, METH_VARARGS,
"Count the number of entries in a GTFtree\n"},
    {"isTree", (PyCFunction) pyIsTree, METH_VARARGS,
"Return True if the object is a tree\n"},
    {"hasOverlaps", (PyCFunction) pyHasOverlaps, METH_VARARGS,
"Returns a tuple with the first value True if ANY of the entries in the tree overlap (ignoring strand) and False otherwise. The second value in the tuple is the minimum distance between intervals (0 on overlap).\n"},
    {"findOverlaps", (PyCFunction) pyFindOverlaps, METH_VARARGS,
"Find overlapping intervals\n"},
    {"findOverlappingFeatures", (PyCFunction) pyFindOverlappingFeatures, METH_VARARGS,
"Find overlapping intervals, returning a list of features\n"},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
struct treemodule_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct treemodule_state*)PyModule_GetState(m))

static PyModuleDef treemodule = {
    PyModuleDef_HEAD_INIT,
    "tree",
    "A python module creating/accessing GTF-based interval trees with associated meta-data",
    -1,
    treeMethods,
    NULL, NULL, NULL, NULL
};
#endif


//Should set tp_dealloc, tp_print, tp_repr, tp_str, tp_members
static PyTypeObject pyGTFtree = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,              /*ob_size*/
#endif
    "pyGTFtree",          /*tp_name*/
    sizeof(pyGTFtree_t),           /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)pyGTFDealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash*/
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    PyObject_GenericGetAttr,   /*tp_getattro*/
    PyObject_GenericSetAttr,   /*tp_setattro*/
    0,                         /*tp_as_buffer*/
#if PY_MAJOR_VERSION >= 3
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
#else
    Py_TPFLAGS_HAVE_CLASS,     /*tp_flags*/
#endif
    "GTF tree",                /*tp_doc*/
    0,                         /*tp_traverse*/
    0,                         /*tp_clear*/
    0,                         /*tp_richcompare*/
    0,                         /*tp_weaklistoffset*/
    0,                         /*tp_iter*/
    0,                         /*tp_iternext*/
    treeMethods,               /*tp_methods*/
    0,                         /*tp_members*/
    0,                         /*tp_getset*/
    0,                         /*tp_base*/
    0,                         /*tp_dict*/
    0,                         /*tp_descr_get*/
    0,                         /*tp_descr_set*/
    0,                         /*tp_dictoffset*/
    0,                         /*tp_init*/
    0,                         /*tp_alloc*/
    0,                         /*tp_new*/
    0,0,0,0,0,0
};
