%{
#include <vector>
#include <string>
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Exception.H"
#include "AddOns/Python/MEProcess.H"
%}

%catches (ATOOLS::Exception) MEProcess::Initialize();

namespace SHERPA{
  class Sherpa;
}
namespace ATOOLS{
  class Cluster_Amplitude;
  class ColorID;
}
namespace PHASIC{
  class Process_Base;
}

class MEProcess{

public:

  MEProcess(SHERPA::Sherpa* Generator);
  ~MEProcess();
  void AddInFlav(const int &id);
  void AddOutFlav(const int &id);
  void AddInFlav(const int &id, const int &col1, const int &col2);
  void AddOutFlav(const int &id, const int &col1, const int &col2);
  double GenerateColorPoint();
  void Initialize();

  void SetMomentum(int, double, double, double, double);
  void SetMomenta(ATOOLS::Vec4D_Vector&);

  // Get the momenta that were previously set
  ATOOLS::Vec4D_Vector GetMomenta();

  double TestPoint(const double& sqrts);
  double MatrixElement();
  double CSMatrixElement();
  double MEProcess::GetFlux();
  inline ATOOLS::Cluster_Amplitude* GetAmp()
  {return m_amp;}

  %extend {
    PyObject* SetMomenta(PyObject* list_list) {
      if (!PyList_Check(list_list)){
	PyErr_SetString(PyExc_TypeError,"Argument of SetMomenta must be a list of lists");
	return NULL;
      }
      ATOOLS::Vec4D_Vector vec4_vec;
      for (int i(0); i<PySequence_Length(list_list); i++)
	{
	  PyObject* momentum = PySequence_GetItem(list_list,i);
	  if(!PyList_Check(momentum)){
	    PyErr_SetString(PyExc_TypeError,"Argument of SetMomenta must be a list of lists");
	    return NULL;
	  }
	  if(PySequence_Length(momentum)!=4){
	    PyErr_SetString(PyExc_TypeError,"Momenta must have four components");
	    return NULL;
	  }
	  vec4_vec.push_back(ATOOLS::Vec4D( PyFloat_AsDouble(PySequence_GetItem(momentum,0)),
					    PyFloat_AsDouble(PySequence_GetItem(momentum,1)),
					    PyFloat_AsDouble(PySequence_GetItem(momentum,2)),
					    PyFloat_AsDouble(PySequence_GetItem(momentum,3)) ));
	  $self->SetMomenta(vec4_vec);
	}
      return PyInt_FromLong(1);
    };
  }

};

