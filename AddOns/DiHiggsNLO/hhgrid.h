#ifndef _hhgrid_H_
#define _hhgrid_H_

#include <Python.h>

PyObject* grid_initialize(const char* grid_name);
double grid_virt(PyObject* grid, double s, double t);
int python_check_errors(); 

#endif
