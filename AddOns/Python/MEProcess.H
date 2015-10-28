#ifndef MEProcess_H
#define MEProcess_H

#include <vector>
#include <string>

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

namespace SHERPA{
  class Sherpa;
}
namespace ATOOLS{
  class Cluster_Amplitude;
  class ColorID;
  class Mass_Selector;
}
namespace PHASIC{
  class Process_Base;
  class Rambo;
}

class MEProcess{
  
private:

  std::string m_name;
  ATOOLS::Cluster_Amplitude* p_amp;
  //ATOOLS::Mass_Selector      m_msel;
  SHERPA::Sherpa*            p_gen;
  PHASIC::Process_Base *     p_proc;
  PHASIC::Rambo*             p_rambo;

  size_t                         m_ncolinds;
  std::vector<std::vector<int> > m_colcombinations;
  std::vector<int>               m_gluinds, m_quainds, m_quabarinds;
  std::vector<int>               m_inpdgs, m_outpdgs;
  std::vector<size_t>            m_mom_inds;

  // Store flavs here separately in the same order
  // in which they were added (for performance reasons)
  ATOOLS::Flavour_Vector         m_flavs; 

  size_t m_npsp,m_nin,m_nout;

  void SetMomentumIndices(const std::vector<int> &pdgs);
  PHASIC::Process_Base* FindProcess();

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

  size_t NumberOfPoints();

  void SetMomenta(size_t n);
  void SetMomenta(const std::vector<double*> &p);
  void SetMomenta(const ATOOLS::Vec4D_Vector &p);
  void SetMomentum(const size_t &id, const double &E , const double &px,
                                     const double &py, const double &pz);
  void SetMomentum(const size_t &id, const ATOOLS::Vec4D &p);
  void SetColor   (const size_t &id, const ATOOLS::ColorID &col);

  ATOOLS::Vec4D_Vector TestPoint(const double& sqrts);
  double MatrixElement();
  double CSMatrixElement();

  double GetFlux();
  std::string GeneratorName();

  ATOOLS::Flavour GetFlav(size_t i);
  inline ATOOLS::Flavour GetInFlav(size_t i)  { return GetFlav(i); }
  inline ATOOLS::Flavour GetOutFlav(size_t i) { return GetFlav(i+m_nin); }

  inline ATOOLS::Cluster_Amplitude* GetAmp() { return p_amp; }
  inline std::string Name() { return m_name; }
};

#endif
