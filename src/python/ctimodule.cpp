#ifndef WITHOUT_PYTHON 
#include <Python.h>
#endif

// for PYEXC use PyExc_TypeError, PyExc_ValueError, or PyExc_SystemError, 
#define PYERR(PYEXC,MESSAGE) { std::stringstream ss; ss << MESSAGE; PyErr_SetString(PYEXC,ss.str().c_str()); return NULL; }

#include <iostream>
#include "Chemtable.hpp"

AbstractChemtable3D* chemtable = NULL;
double a=1;

int factorial_(int n) {
  a += 1.0;
  std::cout << "a increased to " << a << std::endl;
  if(n<=1)
    return 1;
  else
    return n * factorial_(n-1);
}

double pow3_(double d) {
  return d*a;
}

void chemtable_load_(const char* filename) { 
  if (chemtable != NULL) {
    delete chemtable;
    chemtable = NULL;
  }
  initChemtable(chemtable,filename);
}

double chemtable_lookup_(const char* name,double Z,double Zvar,double C) { 
  std::cout << "chemtable_lookup_: " << name << " Z,Zvar,C: " << Z << " " << Zvar << " " << C << std::endl; 
  assert(chemtable != NULL);
  double result;
  chemtable->lookup(&result,name,&Z,&Zvar,&C,1);
  return result;
}

double chemtable_pressure_() { 
  assert(chemtable != NULL);
  return chemtable->pressure;
}


#ifndef WITHOUT_PYTHON 
static PyObject* factorial(PyObject* self, PyObject* args){
  int n;
  if (!PyArg_ParseTuple(args,"i",&n))
    return NULL;
  std::cout << "JUST GOT n: " << n << std::endl;
  int result = factorial_(n);
  return Py_BuildValue("i",result);
}

static PyObject* pow3(PyObject* self, PyObject* args){
  double n;
  if (!PyArg_ParseTuple(args,"d",&n))
    return NULL;
  std::cout << "JUST GOT n: " << n << std::endl;
  double result = pow3_(n);
  return Py_BuildValue("d",result);
}

static PyObject* chemtable_load(PyObject* self, PyObject* args){
  const char* filename; 
  if (!PyArg_ParseTuple(args, "s", &filename))
    return NULL;
  try {
    chemtable_load_(filename);
  }
  catch(int ierr) {
    PYERR(PyExc_ValueError,"Cannot load chemtable \"" << filename << "\"");
  }
  int result = 1;
  return Py_BuildValue("i",result);
}

static PyObject* chemtable_lookup(PyObject* self, PyObject* args){
  const char* name; 
  double Z,Zvar,C;
  if (!PyArg_ParseTuple(args, "sddd", &name, &Z, &Zvar, &C))
    return NULL;
  if (chemtable == NULL) {
    PYERR(PyExc_SystemError,"No chemtable loaded. Use chemtable_load(...");
  }
  double result;
  try {
    result = chemtable_lookup_(name,Z,Zvar,C);
  }
  catch(int ierr) {
    PYERR(PyExc_ValueError,"No chemtable data available for \"" << name << "\"");
  }
  return Py_BuildValue("d",result);
}

static PyObject* chemtable_pressure(PyObject* self, PyObject* args){
  const char* name; 
  double Z,Zvar,C;
  if (!PyArg_ParseTuple(args,""))
    return NULL;
  if (chemtable == NULL) {
    PYERR(PyExc_SystemError,"No chemtable loaded. Use chemtable_load(...");
  }
  double result = chemtable_pressure_();
  return Py_BuildValue("d",result);
}

static PyMethodDef mainMethods[] = {
  {"factorial",factorial,METH_VARARGS,"Calculate the factorial of n"},
  {"pow3",pow3,METH_VARARGS,"Calculate the pow n**3"},
  {"chemtable_load",chemtable_load,METH_VARARGS,"loads a particular chemtable"},
  {"chemtable_lookup",chemtable_lookup,METH_VARARGS,"lookup a var in the chemtable"},
  {"chemtable_pressure",chemtable_pressure,METH_VARARGS,"return chemtable pressure"},
  {NULL,NULL,0,NULL}
};

static PyModuleDef cti = {
  PyModuleDef_HEAD_INIT,
  "cti","Cascade Technologies Inc python interface",
  -1,
 mainMethods
};

PyMODINIT_FUNC PyInit_cti(void){
  return PyModule_Create(&cti);
}
#endif

#ifdef WITHOUT_PYTHON
int main(int argc, char* argv[]) {
  
  if (argc != 2) {
    std::cout << "Usage:\n   ./a.out chemtable-name" << std::endl;
    return -1;
  }
  
  chemtable_load_(argv[1]);
  cout << "\nrho lookup: " << chemtable_lookup_("rho",0.0,0.0,0.0) << endl;
  
  return 0;
}
#endif
