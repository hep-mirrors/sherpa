#ifndef YFS_NLO_BASE_H
#define YFS_NLO_BASE_H

#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/FSR.H"
#include "YFS/NLO/Virtual.H"
#include "YFS/NLO/Real.H"
#include "YFS/NLO/RealReal.H"
#include "YFS/NLO/RealVirtual.H"
#include "YFS/NLO/VirtualVirtual.H"
#include "YFS/Main/YFS_Form_Factor.H"
#include "YFS/Main/Define_Dipoles.H"
#include "ATOOLS/Math/Histogram_2D.H"
#include "ATOOLS/Math/Histogram.H"

#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "PHASIC++/Channels/Multi_Channel.H"




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
  YFS::VirtualVirtual *p_vv;
  YFS::Define_Dipoles *p_nlodipoles;
  YFS::YFS_Form_Factor *p_yfsFormFact;
  YFS::FSR  *p_fsr;
  PHASIC::Single_Process *p_born;
  YFS::Define_Dipoles *p_dipoles;
  Vec4D_Vector m_plab, m_reallab, m_FSRPhotons, m_ISRPhotons, m_bornMomenta, m_eikmom, m_borngamma;
  Vec4D m_nvec, m_realsum, m_ksum;
  ATOOLS::Flavour_Vector m_flavs;
  histograms2d m_histograms2d;
  histograms1d m_histograms1d;
  ATOOLS::Momenta_Stretcher m_stretcher;
  PHASIC::Multi_Channel *p_fsrchannel;
  PHASIC::Cut_Data *p_cuts;

  double m_evts, m_recola_evts, m_ranTheta, m_ranPhi, m_real,m_ravg, m_eikeex, m_subloc;
  int m_is_isr;
  bool m_failcut;

  void Init(Flavour_Vector &flav, Vec4D_Vector &plab, Vec4D_Vector &born);
  void InitializeVirtual(const PHASIC::Process_Info& pi);
  void InitializeReal(const PHASIC::Process_Info& pi);
  void InitializeRealReal(const PHASIC::Process_Info& pi);
  void InitializeRealVirtual(const PHASIC::Process_Info& pi);
  void InitializeVV(const PHASIC::Process_Info& pi);
  double CalculateVirtual();
  double CalculateReal();
  double CalculateRealReal();
  double CalculateReal(Vec4D k, int submode=0);
  double CalculateRealReal(Vec4D k1, Vec4D k2, int fsrmode1=0, int fsrmode2=0);
  double CalculateRealVirtual();
  double CalculateRealVirtual(Vec4D k, int submode=0);
  double CalculateVV();
  void MakeHardMomenta();
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k);
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k1, ATOOLS::Vec4D &k2);
  void MapInitial(ATOOLS::Vec4D_Vector &p);
  void CheckMassReg();
  double CalculateNLO();
  bool CheckPhotonForReal(const Vec4D &k);
  bool CheckPhotonForReal(const Vec4D &k, const Vec4D_Vector &p);
  void CheckMasses(Vec4D_Vector &p, int realmode=0);
  bool CheckMomentumConservation(Vec4D_Vector p);
  void CheckRealSub(Vec4D k, int mode);
  void CheckRealVirtualSub(Vec4D k);
  void CheckRealRealSub(Vec4D k1, Vec4D k2, int fsr1, int fsr2);
  void RandomRotate(Vec4D &k);
  void RescaleMasses(Vec4D_Vector &p, std::vector<double> masses);
  inline void SetBorn(double b) { m_born = b;}
  inline void SetMomentum(Vec4D_Vector &p) { m_plab = p;}
  inline void SetBornProc(PHASIC::Single_Process *p) {p_born = p;}
  inline void SetFSR(YFS::FSR *p) {p_fsr = p;}
  inline bool HasVirtual() {return m_looptool;}
  inline bool HasReal() {return m_realtool;}
  inline bool HasRR() {return m_rrtool;}
};


}
#endif