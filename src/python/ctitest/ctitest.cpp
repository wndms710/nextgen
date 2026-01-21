#include <Python.h>

int factorial_(int n){
  if(n<=1)
    return 1;
  else
   return n * factorial_(n-1);
}

static PyObject* factorial(PyObject* self, PyObject* args){
  int n;
  if (!PyArg_ParseTuple(args,"i",&n))
    return NULL;
  int result = factorial_(n);
  return Py_BuildValue("i",result);
}

static PyMethodDef mainMethods[] = {
  {"factorial",factorial,METH_VARARGS,"Calculate the factorial of int n"},
  {NULL,NULL,0,NULL}
};

static PyModuleDef ctitest = {
  PyModuleDef_HEAD_INIT,
  "ctitest","Factorial Calculation",
  -1,
 mainMethods
};

PyMODINIT_FUNC PyInit_ctitest(void){
  return PyModule_Create(&ctitest);
}
