%module Sherpa

%include "std_vector.i"
namespace std {
  %template(vectori) vector<int>;
  %template(vectord) vector<double>;
};

%include "Terminator_Objects.i"
%include "Exception.i"
%include "Flavour.i"
%include "Vec4.i"
%include "Particle.i"
%include "Blob.i"
%include "Blob_List.i"
%include "MEProcess.i"
%include "Random.i"
%include "Model_Base.i"
%include "Rambo.i"

%{
#include <SHERPA/Main/Sherpa.H>
  // Inclusion of these headers is required here because
  // of the static pointers to RNG and model that are
  // made available:
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Model_Base.H"

#include <cstring>
  %}

%catches (const ATOOLS::Exception&) SHERPA::Sherpa::InitializeTheRun();

// A typemap is required in order to be able to pass
// the python arguments to SHERPA::Sherpa::Sherpa()
%typemap(in) char ** {
  // Check if is a list
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **) malloc((size+1)*sizeof(char *));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($input,i);
      Py_ssize_t str_size = 0;
      const char* str = PyUnicode_AsUTF8AndSize(o, &str_size);
      $1[i] = (char *)malloc((str_size + 1) * sizeof(char));
      std::strcpy($1[i], str);
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

  class Sherpa : public ATOOLS::Terminator_Object {
    
  public:
    Sherpa(int, char**);
    ~Sherpa();
    bool InitializeTheRun();
    bool SummarizeRun();
    bool GenerateOneEvent();
    bool InitializeTheEventHandler();
    long int NumberOfEvents() const;
    const ATOOLS::Blob_List &GetBlobList() const;
    double GetMEWeight(const ATOOLS::Cluster_Amplitude &ampl,const int mode=0) const;
    double TotalXS();
    double TotalErr();
    
  };
}

// Make the global pointers
// to RNG and model availeble
namespace ATOOLS {
  extern Random* ran;
}
namespace MODEL {
  Model_Base* s_model;
}
