#include "ATOOLS/Math/Tensor.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/UFO/UFO_Model.H"
#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Org/Data_Reader.H"

#define CF 1.33333333333333333
#define TR 0.5

using namespace ATOOLS;
using namespace EXTRAXS;
using namespace PHASIC;
using namespace std;

//------------------------------------------------------------------------------
// gq -> ll q
//------------------------------------------------------------------------------

namespace EXTRAXS {

  class XS_gqllq_CSS_approx : public ME2_Base {
  private:
    EXTRAXS::ME2_Base * p_bornme;
    bool   m_anti;
    double m_alphasdef;

  public:

    XS_gqllq_CSS_approx(const Process_Info& pi,const Flavour_Vector& fl);
    ~XS_gqllq_CSS_approx();

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    double LOME2(const Vec4D&, const Vec4D&, const Vec4D&, const Vec4D&,
                 const Vec4D&, int);
  };
}

XS_gqllq_CSS_approx::XS_gqllq_CSS_approx
(const Process_Info& pi, const Flavour_Vector& fl) : ME2_Base(pi, fl)
{
  Process_Info pico(pi);
  m_anti=false;
  pico.m_fi.m_ps.erase(pico.m_fi.m_ps.end()-1);
  pico.m_ii.m_ps[0].m_fl=pico.m_ii.m_ps[1].m_fl.Bar();
  if (pico.m_ii.m_ps[0].m_fl.IsAnti()) {
    std::swap(pico.m_ii.m_ps[0].m_fl,pico.m_ii.m_ps[1].m_fl);
    m_anti=true;
  }
  pico.m_fi.m_nloqcdtype=nlo_type::lo;
  p_bornme = dynamic_cast<ME2_Base*>(PHASIC::Tree_ME2_Base::GetME2(pico));
  if (!p_bornme) THROW(fatal_error,"no born me found.");
  m_alphasdef = MODEL::as->Default();
  PRINT_INFO("initialised XS_gqllq_CSS_approx2");
}

XS_gqllq_CSS_approx::~XS_gqllq_CSS_approx()
{
  if (p_bornme) delete p_bornme;
}

double XS_gqllq_CSS_approx::operator()
(const ATOOLS::Vec4D_Vector& p)
{
  double res(0);
  // (Gqb) q -> ll
  if (m_anti)
    res+=LOME2(p[0],p[4],p[1],p[2],p[3],1);
  // (Gq) qb -> ll
  else
    res+=LOME2(p[0],p[4],p[1],p[2],p[3],0);
  return res;
}

double XS_gqllq_CSS_approx::LOME2(const Vec4D& pi, const Vec4D& pj,
                                  const Vec4D& pk, const Vec4D& k1,
                                  const Vec4D& k2, int ij)
{
  DEBUG_FUNC("");
//  msg_Debugging()<<"pi: "<<pi<<std::endl;
//  msg_Debugging()<<"pj: "<<pj<<std::endl;
//  msg_Debugging()<<"pk: "<<pk<<std::endl;
//  msg_Debugging()<<"k1: "<<k1<<std::endl;
//  msg_Debugging()<<"k2: "<<k2<<std::endl;
  double pipj(pi*pj), pipk(pi*pk), pjpk(pj*pk);
  double xiab=(pipk-pipj-pjpk)/pipk;
  double vi=pipj/pipk;

  Vec4D pijt=xiab*pi;
  Vec4D pkt=pk;

  double kp(sqrt(2.*pipk*vi*(1.-xiab-vi)));
  Vec4D kperp(pj-(1.-xiab-vi)/xiab*pijt-vi*pkt);
  double phi(0.);
  if      ((kperp[1]>=0. && kperp[2]>=0.) || (kperp[1]>=0. && kperp[2]<0.))
    phi=acos(kperp[2]/kp);
  else if ((kperp[1]<0. && kperp[2]<0.) || (kperp[1]<0. && kperp[2]>=0.))
    phi=-acos(kperp[2]/kp)+2.*M_PI;
  else THROW(fatal_error,"Could not determine phi.");

  Vec4D K(pi-pj+pk), Kt(pijt+pkt);
  ATOOLS::Lorentz_Ten2D Lambda = MetricTensor()
                                 - 2./(K+Kt).Abs2()*BuildTensor(Kt+K,Kt+K)
                                 + 2./Kt.Abs2()*BuildTensor(Kt,K);

  Vec4D k1t = Contraction(Lambda,2,k1);
  Vec4D k2t = Contraction(Lambda,2,k2);

  msg_Debugging()<<"pijt: "<<pijt<<std::endl;
  msg_Debugging()<<"pkt:  "<<pkt<<std::endl;
  msg_Debugging()<<"k1t:  "<<k1t<<std::endl;
  msg_Debugging()<<"k2t:  "<<k2t<<std::endl;

  Vec4D_Vector moms(4);
  moms[ij]=pijt;
  moms[1-ij]=pkt;
  moms[2]=k1t;
  moms[3]=k2t;
  double born(-48.*(*p_bornme)(moms));
  // SF = 8*pi*TR/(2pipj*x) * (1-2x(1-x))
  double split(8.*M_PI/((pi+pj).Abs2()*xiab)*(1.-2.*xiab*(1.-xiab)));
  msg_Debugging()<<8.*M_PI*m_alphasdef<<std::endl;
  msg_Debugging()<<"M2 = "<<born<<" ,  SF = "<<split*m_alphasdef<<std::endl;
  return -born*split*m_alphasdef;
}

DECLARE_TREEME2_GETTER(XS_gqllq_CSS_approx,"XS_gqllq_CSS_approx")
Tree_ME2_Base *ATOOLS::Getter
<Tree_ME2_Base,Process_Info,XS_gqllq_CSS_approx>::
operator()(const Process_Info &pi) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  Data_Reader read(" ",";","!","=");
  if (read.GetValue<int>("EXTRAXS_CSS_APPROX_ME",0)==0) return NULL;
  if (pi.m_fi.NLOType()!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=5) return NULL;
  if (fl[1].IsQuark()  && fl[4]==fl[1] &&
      fl[0].IsGluon()  &&
      fl[2].IsLepton() && fl[3]==fl[2].Bar()) {
    if (pi.m_maxcpl[0]==1 && pi.m_maxcpl[1]==2 &&
	pi.m_mincpl[0]==1 && pi.m_mincpl[1]==2) {
      return new XS_gqllq_CSS_approx(pi,fl);
    }
  }
  return NULL;
}


//------------------------------------------------------------------------------
// qq -> ll g
//------------------------------------------------------------------------------

namespace EXTRAXS {

  class XS_qqllg_CSS_approx : public ME2_Base {
  private:
    EXTRAXS::ME2_Base * p_bornme;
    double m_alphasdef;

  public:

    XS_qqllg_CSS_approx(const Process_Info& pi,const Flavour_Vector& fl);
    ~XS_qqllg_CSS_approx();

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    double LOME2(const Vec4D&, const Vec4D&, const Vec4D&, const Vec4D&,
                 const Vec4D&, int);
  };
}

XS_qqllg_CSS_approx::XS_qqllg_CSS_approx
(const Process_Info& pi, const Flavour_Vector& fl) : ME2_Base(pi, fl)
{
  Process_Info pico(pi);
  pico.m_fi.m_ps.erase(pico.m_fi.m_ps.end()-1);
  pico.m_fi.m_nloqcdtype=nlo_type::lo;
  p_bornme = dynamic_cast<ME2_Base*>(PHASIC::Tree_ME2_Base::GetME2(pico));
  m_alphasdef = MODEL::as->Default();
  PRINT_INFO("initialised XS_qqllg_CSS_approx2");
}

XS_qqllg_CSS_approx::~XS_qqllg_CSS_approx()
{
  if (p_bornme) delete p_bornme;
}

double XS_qqllg_CSS_approx::operator()
(const ATOOLS::Vec4D_Vector& p)
{
  double res(0);
  // (qG) qb -> ll
  res+=LOME2(p[0],p[4],p[1],p[2],p[3],0);
  // q (qbG) -> ll
  res+=LOME2(p[1],p[4],p[0],p[2],p[3],1);
  return res;
}

double XS_qqllg_CSS_approx::LOME2(const Vec4D& pi, const Vec4D& pj,
                                  const Vec4D& pk, const Vec4D& k1,
                                  const Vec4D& k2, int ij)
{
  DEBUG_FUNC("");
//  msg_Debugging()<<"pi: "<<pi<<std::endl;
//  msg_Debugging()<<"pj: "<<pj<<std::endl;
//  msg_Debugging()<<"pk: "<<pk<<std::endl;
//  msg_Debugging()<<"k1: "<<k1<<std::endl;
//  msg_Debugging()<<"k2: "<<k2<<std::endl;
  double pipj(pi*pj), pipk(pi*pk), pjpk(pj*pk);
  double xiab=(pipk-pipj-pjpk)/pipk;
  double vi=pipj/pipk;

  Vec4D pijt=xiab*pi;
  Vec4D pkt=pk;

  double kp(sqrt(2.*pipk*vi*(1.-xiab-vi)));
  Vec4D kperp(pj-(1.-xiab-vi)/xiab*pijt-vi*pkt);
  double phi(0.);
  if      ((kperp[1]>=0. && kperp[2]>=0.) || (kperp[1]>=0. && kperp[2]<0.))
    phi=acos(kperp[2]/kp);
  else if ((kperp[1]<0. && kperp[2]<0.) || (kperp[1]<0. && kperp[2]>=0.))
    phi=-acos(kperp[2]/kp)+2.*M_PI;
  else THROW(fatal_error,"Could not determine phi.");

  Vec4D K(pi-pj+pk), Kt(pijt+pkt);
  ATOOLS::Lorentz_Ten2D Lambda = MetricTensor()
                                 - 2./(K+Kt).Abs2()*BuildTensor(Kt+K,Kt+K)
                                 + 2./Kt.Abs2()*BuildTensor(Kt,K);

  Vec4D k1t = Contraction(Lambda,2,k1);
  Vec4D k2t = Contraction(Lambda,2,k2);

  msg_Debugging()<<"pijt: "<<pijt<<std::endl;
  msg_Debugging()<<"pkt:  "<<pkt<<std::endl;
  msg_Debugging()<<"k1t:  "<<k1t<<std::endl;
  msg_Debugging()<<"k2t:  "<<k2t<<std::endl;

  Vec4D_Vector moms(4);
  moms[ij]=pijt;
  moms[1-ij]=pkt;
  moms[2]=k1t;
  moms[3]=k2t;
  double born(-48.*(*p_bornme)(moms));
  // SF = 8*pi*CF/(2pipj*x) * (2/(1-x)-(1+x))
  double split(8.*M_PI/((pi+pj).Abs2()*xiab)*(2./(1.-xiab)-(1.+xiab)));
  msg_Debugging()<<8.*M_PI*m_alphasdef<<std::endl;
  msg_Debugging()<<"M2 = "<<born<<" ,  SF = "<<split*m_alphasdef<<std::endl;
  return -born*split*m_alphasdef;
}

DECLARE_TREEME2_GETTER(XS_qqllg_CSS_approx,"XS_qqllg_CSS_approx")
Tree_ME2_Base *ATOOLS::Getter
<Tree_ME2_Base,Process_Info,XS_qqllg_CSS_approx>::
operator()(const Process_Info &pi) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  Data_Reader read(" ",";","!","=");
  if (read.GetValue<int>("EXTRAXS_CSS_APPROX_ME",0)==0) return NULL;
  if (pi.m_fi.NLOType()!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=5) return NULL;
  if (fl[0].IsQuark()  && fl[1]==fl[0].Bar() &&
      fl[4].IsGluon()  &&
      fl[2].IsLepton() && fl[3]==fl[2].Bar()) {
    if (pi.m_maxcpl[0]==1 && pi.m_maxcpl[1]==2 &&
	pi.m_mincpl[0]==1 && pi.m_mincpl[1]==2) {
      return new XS_qqllg_CSS_approx(pi,fl);
    }
  }
  return NULL;
}



