#include "hhgrid.h"

#include <stdlib.h> 
#include <assert.h>
#include <unistd.h>


int python_check_errors()
{
  if(!PyErr_Occurred()) return 0;
  PyErr_Print();
  return 1;
}


PyObject* grid_initialize(const char* grid_name)
{
    if( access(grid_name, F_OK) == -1 || access(grid_name, R_OK) == -1 )
      {
	printf("ERROR: Failed to read grid at %s\n", grid_name);
	return NULL;
      }

    PyObject* pModule = PyImport_ImportModule("creategrid");
    if (python_check_errors()) return NULL;

    PyObject* pClass = PyObject_GetAttrString(pModule, "CreateGrid");
    if (python_check_errors()) return NULL;

    PyObject* pGridName = PyString_FromString(grid_name);
    if (python_check_errors()) return NULL;

    PyObject* pGridNameTuple = PyTuple_Pack(1,pGridName);
    if (python_check_errors()) return NULL;

    PyObject* pInstance = PyInstance_New(pClass, pGridNameTuple, NULL);
    if (python_check_errors()) return NULL;
    
    assert(pModule        != NULL);
    assert(pClass         != NULL);
    assert(pGridName      != NULL);
    assert(pGridNameTuple != NULL);
    assert(pInstance      != NULL);

    Py_DECREF(pModule);
    Py_DECREF(pClass);
    Py_DECREF(pGridName);
    Py_DECREF(pGridNameTuple);

    return pInstance;
};


double grid_virt(PyObject* grid, double s, double t)
{
    PyObject* pResult = PyObject_CallMethod(grid, "GetAmplitude", "(ff)", s, t);
    assert(pResult != NULL);
    double result = PyFloat_AsDouble(pResult);
    Py_DECREF(pResult);
    return result;
};
