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




namespace YFS {

class NLO_Base : YFS_Base {
public:

  NLO_Base();
  ~NLO_Base();

  YFS::Virtual     *p_virt;
  YFS::Real        *p_real;
  YFS::RealReal    *p_realreal;
  YFS::RealVirtual *p_realvirt;
  Griffin::Griffin_Interface *p_griffin;
  YFS::YFS_Form_Factor *p_yfsFormFact;
  YFS::Define_Dipoles *p_dipoles;
  Vec4D_Vector m_plab, m_reallab, m_FSRPhotons, m_ISRPhotons, m_bornMomenta;
  ATOOLS::Flavour_Vector m_flavs;

  void Init(Flavour_Vector &flav, Vec4D_Vector &plab, Vec4D_Vector &born);
  void InitializeVirtual(const PHASIC::Process_Info& pi);
  void InitializeReal(const PHASIC::Process_Info& pi);
  void InitializeRealReal(const PHASIC::Process_Info& pi);
  void InitializeRealVirtual(const PHASIC::Process_Info& pi);
  double CalculateVirtual();
  double CalculateReal();
  double CalculateReal(Vec4D &k);
  double CalculateRealVirtual();
  double CalculateRealReal();
  double CollinearReal(Vec4D &k, Vec4D_Vector &p);

  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k);
  double CalculateNLO();
  bool CheckPhotonForReal(const Vec4D &k);

  inline void SetBorn(double b) { m_born = b;}

};


}
#endif