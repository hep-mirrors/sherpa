%module Sherpa
%include "Exception.i"
%include "Initialization_Handler.i"
%include "Matrix_Element_Handler.i"
%include "Cluster_Amplitude.i"
%include "Process_Base.i"
%include "Message.i"
%include "Flavour.i"
%include "Vec4.i"
%include "iostream.i"
%{
#include <SHERPA/Main/Sherpa.H>
  %}

%catches (ATOOLS::Exception) SHERPA::Sherpa::InitializeTheRun(int, char**);

// A typemap is required in order to be able to pass
// the python arguments to SHERPA::Sherpa::InitializeTheRun(int, char**)
%typemap(in) char ** {
  // Check if is a list
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **) malloc((size+1)*sizeof(char *));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o))
	$1[i] = PyString_AsString(PyList_GetItem($input,i));
      else {
	PyErr_SetString(PyExc_TypeError,"Sherpa execution argument list must contain strings");
	free($1);
	return NULL;
      }
    }
    $1[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,"Sherpa execution argument is not a list");
    return NULL;
  }
 }

// This cleans up the char ** array for which memory was allocated
%typemap(freearg) char ** {
  free((char *) $1);
 }

namespace SHERPA {

  class Initialization_Handler;

  class Sherpa : public ATOOLS::Terminator_Object {
    
  public:
    Sherpa();
    ~Sherpa();
    bool InitializeTheRun(int, char**);
    bool SummarizeRun();
    bool GenerateOneEvent();
    bool InitializeTheEventHandler();
    inline Initialization_Handler * GetInitHandler() const;
    
  };
}
