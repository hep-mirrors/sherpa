#ifndef MEProcess_H
#define MEProcess_H

#include <vector>
#include <string>

#include "ATOOLS/Phys/Flavour.H"

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
  
private:
  
  ATOOLS::Cluster_Amplitude* m_amp;
  MEProcess();
  SHERPA::Sherpa* m_gen;
  std::string m_name;
  PHASIC::Process_Base *m_proc;
  std::vector<std::vector<int> > m_colcombinations;
  std::vector<int> m_gluinds;
  std::vector<int> m_quainds;
  std::vector<int> m_quabarinds;
  std::vector<int> m_inpdgs;
  std::vector<int> m_outpdgs;
  std::vector<unsigned int> m_mom_inds;

  int m_ncolinds;
  unsigned int m_nin;
  unsigned int m_nout;

  void SetMomentumIndices(const std::vector<int> &pdgs);

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

  void SetMomenta(const std::vector<double*> &p);
  void SetMomentum(int, double, double, double, double);

  double MatrixElement();
  double CSMatrixElement();
  inline ATOOLS::Cluster_Amplitude* GetAmp()
  {return m_amp;}
};

#endif
