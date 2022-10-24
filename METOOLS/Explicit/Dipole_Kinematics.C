#include "METOOLS/Explicit/Dipole_Kinematics.H"

#include "METOOLS/Explicit/Vertex.H"
#include "PDF/Main/NLOMC_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE METOOLS::Vertex_Key
#define OBJECT_TYPE    METOOLS::Dipole_Kinematics
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

using namespace METOOLS;
using namespace ATOOLS;


Dipole_Kinematics::Dipole_Kinematics
(Dipole_Info *const info,Current *const i,Current *const j,
 Current *const k,Current *const ijt,Current *const kt):
  p_i(i), p_j(j), p_k(k), p_ijt(ijt), p_kt(kt),
  m_type(0), m_swap(0), m_trig(1),
  p_info(info),
  m_mi2(0.0), m_mj2(0.0), m_mij2(0.0), m_mk2(0.0),
  m_mi(0.0), m_mj(0.0), m_mij(0.0), m_mk(0.0),
  m_ym(0.0), m_yp(1.0), m_f(0.0), m_a(0.0),
  p_subevt(NULL), p_nlomc(NULL)
{
  if (p_i) m_mi2=sqr(m_mi=p_i->Flav().Mass());
  if (p_j) m_mj2=sqr(m_mj=p_j->Flav().Mass());
  m_phase[1]=m_phase[0]=0.0;
  m_res[2]=m_res[1]=m_res[0]=0.0;
  if ((p_i && p_i->Direction()==0) ||
      (p_j && p_j->Direction()==0) ||
      p_k->Direction()==0)
    THROW(fatal_error,"Missing current information");
  if (p_k->Direction()>0) m_type|=2;
  if ((p_i && p_i->Direction()>0) ||
      (p_j && p_j->Direction()>0)) {
    if (p_i && p_j && p_j->Direction()>0) {
      std::swap<Current*>(p_i,p_j);
      m_swap=1;
    }
    m_type|=1;
  }
  m_mij2 = sqr(m_mij=p_ijt->Flav().Mass());
  m_mk2 = sqr(m_mk=p_k->Flav().Mass());
}

Dipole_Kinematics::~Dipole_Kinematics()
{
}


void Dipole_Kinematics::SetNLOMC(PDF::NLOMC_Base *const mc)
{
  p_nlomc=mc;
  p_info->SetSubType(p_nlomc->SubtractionType());
  if (p_nlomc->SubtractionType()==1) p_info->SetKappa(1.0);
}


std::ostream &METOOLS::operator<<
  (std::ostream &str,const Dipole_Kinematics &k)
{
  return str<<*k.JI()<<","<<*k.JJ()<<"<->"<<*k.JK()<<" "<<k.Type();
}


// TODO: do this right
#include "METOOLS/Explicit/CS_Dipole_Kinematics.H"

DECLARE_GETTER(Dipole_Kinematics,"Dipole_Kinematics",Dipole_Kinematics,Vertex_Key);

Dipole_Kinematics * ATOOLS::Getter<Dipole_Kinematics,Vertex_Key,Dipole_Kinematics>::operator()(const Parameter_Type & info) const {
  
  return new CS_Dipole_Kinematics(info.p_dinfo,info.m_j[0],info.m_j[1],info.p_k,info.p_c,info.p_kt);
}

void Getter<Dipole_Kinematics,Vertex_Key,Dipole_Kinematics>::PrintInfo(std::ostream &str,const size_t width) const {
  str<<"Dipole Kinematics\n";
}
