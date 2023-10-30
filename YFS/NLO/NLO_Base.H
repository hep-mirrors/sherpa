#ifndef YFS_NLO_BASE_H
#define YFS_NLO_BASE_H

#include "YFS/Main/YFS_Base.H"
#include "YFS/NLO/Virtual.H"
#include "YFS/NLO/Real.H"
#include "YFS/NLO/RealReal.H"
#include "YFS/NLO/RealVirtual.H"
#include "YFS/Main/Griffin_Interface.H"
#include "YFS/Main/YFS_Form_Factor.H"
#include "YFS/Main/Define_Dipoles.H"
#include "ATOOLS/Math/Histogram_2D.H"
#include "ATOOLS/Math/Histogram.H"

#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"




namespace YFS {
typedef std::map<std::string, ATOOLS::Histogram_2D *> histograms2d;
typedef std::map<std::string, ATOOLS::Histogram *> histograms1d;

class NLO_Base : YFS_Base {
public:

  NLO_Base();
  ~NLO_Base();

  YFS::Virtual     *p_virt;
  YFS::Real        *p_real;
  YFS::RealReal    *p_realreal;
  YFS::RealVirtual *p_realvirt;
  YFS::Define_Dipoles *p_nlodipoles;
  Griffin::Griffin_Interface *p_griffin;
  YFS::YFS_Form_Factor *p_yfsFormFact;
  PHASIC::Single_Process *p_born;
  YFS::Define_Dipoles *p_dipoles, *p_global_dipoles;
  Vec4D_Vector m_plab, m_reallab, m_FSRPhotons, m_ISRPhotons, m_bornMomenta, m_eikmom;
  Vec4D m_nvec;
  ATOOLS::Flavour_Vector m_flavs;
  histograms2d m_histograms2d;
  histograms1d m_histograms1d;
  ATOOLS::Momenta_Stretcher m_stretcher;

  double m_evts, m_recola_evts;
  int m_is_isr;

  void Init(Flavour_Vector &flav, Vec4D_Vector &plab, Vec4D_Vector &born);
  void InitializeVirtual(const PHASIC::Process_Info& pi);
  void InitializeReal(const PHASIC::Process_Info& pi);
  void InitializeRealReal(const PHASIC::Process_Info& pi);
  void InitializeRealVirtual(const PHASIC::Process_Info& pi);
  double CalculateVirtual();
  double CalculateReal();
  double CalculateReal(Vec4D k, int submode=0);
  double CalculateRealVirtual();
  double CalculateRealReal();
  void MakeHardMomenta();
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k);
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k, ATOOLS::Vec4D &kk);
  void MapInitial(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k);
  void MapFinal(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k);
  void CheckMassReg();
  void CheckMassRegRV(const ATOOLS::Vec4D_Vector &p);
  double CalculateNLO();
  bool CheckPhotonForReal(const Vec4D &k);
  bool CheckPhotonCollinear(const Vec4D_Vector &p, const Vec4D &k);
  void CheckMasses(Vec4D_Vector &p);
  bool CheckMomentumConservation(Vec4D_Vector p);

  inline void SetBorn(double b) { m_born = b;}
  inline void SetMomentum(Vec4D_Vector &p) { m_plab = p;}
  inline void SetBornProc(PHASIC::Single_Process *p) {p_born = p;}
  inline void SetResonanceMode(int &mode) {m_resonance_mode = mode;}

};


}
#endif