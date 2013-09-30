%{
#include <vector>
#include <string>
#include "ATOOLS/Phys/Flavour.H"
#include "AddOns/Python/MEProcess.H"
%}

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
  bool HasColorIntegrator();
  void SetColors();
  void Initialize();

  void SetMomentum(int, double, double, double, double);

  double MatrixElement();
  double CSMatrixElement();
  inline ATOOLS::Cluster_Amplitude* GetAmp()
  {return m_amp;}

};

