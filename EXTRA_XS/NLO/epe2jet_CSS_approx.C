#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/UFO/UFO_Model.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/External_ME_Args.H"

#define CF 1.33333333333333333
#define TR 0.5

using namespace ATOOLS;
using namespace EXTRAXS;
using namespace PHASIC;
using namespace std;

//------------------------------------------------------------------------------
// e g -> e q qbar
//------------------------------------------------------------------------------

namespace EXTRAXS {

  class XS_egeqq_CSS_approx : public ME2_Base {
  private:
    EXTRAXS::ME2_Base * p_bornme1, * p_bornme2;
    double m_alphasdef;

  public:

    XS_egeqq_CSS_approx(const External_ME_Args& args);
    ~XS_egeqq_CSS_approx();

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    double LOME2(const Vec4D&, const Vec4D&, const Vec4D&, const Vec4D&,
                 const Vec4D&, int);
  };
}

XS_egeqq_CSS_approx::XS_egeqq_CSS_approx
(const External_ME_Args& args) : ME2_Base(args)
{
  PRINT_INFO("initialising XS_egeqq_CSS_approx2");

  Flavour_Vector inflavs1  = args.m_inflavs;
  Flavour_Vector inflavs2  = args.m_inflavs;
  Flavour_Vector outflavs1 = args.m_outflavs;
  Flavour_Vector outflavs2 = args.m_outflavs;

  inflavs1[1] = args.m_outflavs[1];
  inflavs2[1] = args.m_outflavs[2];

  outflavs1.erase(outflavs1.end()-1);
  outflavs2.erase(outflavs2.end()-2);

  p_bornme1 = dynamic_cast<ME2_Base*>
    (PHASIC::Tree_ME2_Base::GetME2(External_ME_Args(inflavs1,
						    outflavs1,
						    {0,2})));

  p_bornme2 = dynamic_cast<ME2_Base*>
    (PHASIC::Tree_ME2_Base::GetME2(External_ME_Args(inflavs2,
						    outflavs2,
						    {0,2})));

  if (!p_bornme1 || !p_bornme2) THROW(fatal_error,"no born me found.");
  m_alphasdef = MODEL::as->Default();
  PRINT_INFO("initialised XS_egeqq_CSS_approx2");
}

XS_egeqq_CSS_approx::~XS_egeqq_CSS_approx()
{
  if (p_bornme1) delete p_bornme1;
  if (p_bornme2) delete p_bornme2;
}

double XS_egeqq_CSS_approx::operator()
(const ATOOLS::Vec4D_Vector& p)
{
  double res(0);
  // e q -> e q
  res+=LOME2(p[1],p[3],p[4],p[0],p[2],1);
  // e qb -> e qb
  res+=LOME2(p[1],p[4],p[3],p[0],p[2],2);
  return res;
}

double XS_egeqq_CSS_approx::LOME2(const Vec4D& pi, const Vec4D& pj,
                                  const Vec4D& pk, const Vec4D& k1,
                                  const Vec4D& k2, int bornterm)
{
  DEBUG_FUNC("");
//  msg_Debugging()<<"pi: "<<pi<<std::endl;
//  msg_Debugging()<<"pj: "<<pj<<std::endl;
//  msg_Debugging()<<"pk: "<<pk<<std::endl;
//  msg_Debugging()<<"k1: "<<k1<<std::endl;
//  msg_Debugging()<<"k2: "<<k2<<std::endl;

  double pipj(pi*pj), pipk(pi*pk), pjpk(pj*pk);
  double xika=(pipj+pipk-pjpk)/(pipj+pipk);

  Vec4D pijt=xika*pi;
  Vec4D pkt=pk+pj-(1.-xika)*pi;

//  msg_Debugging()<<"pijt: "<<pijt<<std::endl;
//  msg_Debugging()<<"pkt:  "<<pkt<<std::endl;
//  msg_Debugging()<<"k1t:  "<<k1<<std::endl;
//  msg_Debugging()<<"k2t:  "<<k2<<std::endl;

  Vec4D_Vector moms(4);
  moms[1]=pijt;
  moms[3]=pkt;
  moms[0]=k1;
  moms[2]=k2;
  msg_Debugging()<<"(0): "<<moms[0]<<std::endl;
  msg_Debugging()<<"(1): "<<moms[1]<<std::endl;
  msg_Debugging()<<"(2): "<<moms[2]<<std::endl;
  msg_Debugging()<<"(3): "<<moms[3]<<std::endl;
  double born(0.);
  if (bornterm==1) born=(*p_bornme1)(moms);
  else born=(*p_bornme2)(moms);
  // SF = 8*pi*TR/(2pipj*x) * (1-2x(1-x))
  double split(8.*M_PI/((pi+pj).Abs2()*xika)*.5*(1.-2.*xika*(1.-xika)));
  msg_Debugging()<<8.*M_PI*m_alphasdef<<std::endl;
  msg_Debugging()<<"M2 = "<<born<<" ,  SF = "<<split*m_alphasdef<<std::endl;
  return born*split*m_alphasdef*CouplingFactor(1,0);
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_egeqq_CSS_approx,"XS_egeqq_CSS_approx")
Tree_ME2_Base *ATOOLS::Getter
<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::XS_egeqq_CSS_approx>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  Settings& s = Settings::GetMainSettings();
  if (!s["EXTRAXS_CSS_APPROX_ME"].Get<bool>()) return NULL;
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=5) return NULL;
  if (fl[0].IsLepton() && fl[2]==fl[0] &&
      fl[1].IsGluon()  &&
      fl[3].IsQuark()  && fl[4]==fl[3].Bar()) {
    if (args.m_orders[0]==1 && args.m_orders[1]==2) {
      return new XS_egeqq_CSS_approx(args);
    }
  }
  return NULL;
}


//------------------------------------------------------------------------------
// e q -> e g q
//------------------------------------------------------------------------------

namespace EXTRAXS {

  class XS_eqegq_CSS_approx : public ME2_Base {
  private:
    EXTRAXS::ME2_Base * p_bornme;
    double m_alphasdef;

  public:

    XS_eqegq_CSS_approx(const External_ME_Args& args);
    ~XS_eqegq_CSS_approx();

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    double LOME2FI(const Vec4D&, const Vec4D&, const Vec4D&, const Vec4D&,
                   const Vec4D&, int);
    double LOME2IF(const Vec4D&, const Vec4D&, const Vec4D&, const Vec4D&,
                   const Vec4D&, int);
  };
}

XS_eqegq_CSS_approx::XS_eqegq_CSS_approx
(const External_ME_Args& args) : ME2_Base(args)
{
  PRINT_INFO("initialising XS_eqegq_CSS_approx2");
  Flavour_Vector outflavs = args.m_outflavs;
  outflavs.erase(outflavs.end()-2);
  External_ME_Args bargs(args.m_inflavs, outflavs,{0,2});
  p_bornme = dynamic_cast<ME2_Base*>(PHASIC::Tree_ME2_Base::GetME2(bargs));
  if (!p_bornme) THROW(fatal_error,"no born me found.");
  m_alphasdef = MODEL::as->Default();
  PRINT_INFO("initialised XS_eqegq_CSS_approx2");
}

XS_eqegq_CSS_approx::~XS_eqegq_CSS_approx()
{
  if (p_bornme) delete p_bornme;
}

double XS_eqegq_CSS_approx::operator()
(const ATOOLS::Vec4D_Vector& p)
{
  double res(0);
  // e q -> e q(qG)
  res+=LOME2IF(p[1],p[3],p[4],p[0],p[2],3);
  // e q(qG) -> e q
  res+=LOME2FI(p[4],p[3],p[1],p[0],p[2],1);
  return res;
}

double XS_eqegq_CSS_approx::LOME2FI(const Vec4D& pi, const Vec4D& pj,
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
  double xija=(pjpk+pipk-pipj)/(pipk+pjpk);
  double zi=pjpk/(pipk+pjpk);

  Vec4D pijt=pi+pj-(1.-xija)*pk;
  Vec4D pkt=xija*pk;

//  msg_Debugging()<<"pijt: "<<pijt<<std::endl;
//  msg_Debugging()<<"pkt:  "<<pkt<<std::endl;
//  msg_Debugging()<<"k1t:  "<<k1<<std::endl;
//  msg_Debugging()<<"k2t:  "<<k2<<std::endl;

  Vec4D_Vector moms(4);
  moms[3]=pijt;
  moms[1]=pkt;
  moms[0]=k1;
  moms[2]=k2;
  msg_Debugging()<<"(0): "<<moms[0]<<std::endl;
  msg_Debugging()<<"(1): "<<moms[1]<<std::endl;
  msg_Debugging()<<"(2): "<<moms[2]<<std::endl;
  msg_Debugging()<<"(3): "<<moms[3]<<std::endl;
  double born((*p_bornme)(moms));
  // SF = 8*pi*CF/(2pipj*x) * (2/(2-z-x)-(1+x))
  double split(8.*M_PI/((pi+pj).Abs2()*xija)*4./3.*(2./(2.-(1.-zi)-xija)-(1.+(1.-zi))));
  msg_Debugging()<<8.*M_PI*m_alphasdef<<std::endl;
  msg_Debugging()<<"M2 = "<<born<<" ,  SF = "<<split*m_alphasdef<<std::endl;
  return born*split*m_alphasdef*CouplingFactor(1,0);
}

double XS_eqegq_CSS_approx::LOME2IF(const Vec4D& pi, const Vec4D& pj,
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
  double xika=(pipj+pipk-pjpk)/(pipj+pipk);
  double ui=pipj/(pipj+pipk);

  Vec4D pijt=xika*pi;
  Vec4D pkt=pk+pj-(1.-xika)*pi;

//  msg_Debugging()<<"pijt: "<<pijt<<std::endl;
//  msg_Debugging()<<"pkt:  "<<pkt<<std::endl;
//  msg_Debugging()<<"k1t:  "<<k1<<std::endl;
//  msg_Debugging()<<"k2t:  "<<k2<<std::endl;

  Vec4D_Vector moms(4);
  moms[1]=pijt;
  moms[3]=pkt;
  moms[0]=k1;
  moms[2]=k2;
  msg_Debugging()<<"(0): "<<moms[0]<<std::endl;
  msg_Debugging()<<"(1): "<<moms[1]<<std::endl;
  msg_Debugging()<<"(2): "<<moms[2]<<std::endl;
  msg_Debugging()<<"(3): "<<moms[3]<<std::endl;
  double born((*p_bornme)(moms));
  // SF = 8*pi*CF/(2pipj*x) * (2/(1-x+u)-(1+x))
  double split(8.*M_PI/((pi+pj).Abs2()*xika)*4./3.*(2./(1.-xika+ui)-(1.+xika)));
  msg_Debugging()<<8.*M_PI*m_alphasdef<<std::endl;
  msg_Debugging()<<"M2 = "<<born<<" ,  SF = "<<split*m_alphasdef<<std::endl;
  return born*split*m_alphasdef*CouplingFactor(1,0);
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_eqegq_CSS_approx,"XS_eqegq_CSS_approx")
Tree_ME2_Base *ATOOLS::Getter
<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::XS_eqegq_CSS_approx>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  Settings& s = Settings::GetMainSettings();
  if (!s["EXTRAXS_CSS_APPROX_ME"].Get<bool>()) return NULL;
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=5) return NULL;
  if (fl[0].IsLepton() && fl[2]==fl[0] &&
      fl[3].IsGluon()  &&
      fl[1].IsQuark()  && fl[4]==fl[1]) {
    if (args.m_orders[0]==1 && args.m_orders[1]==2) {
      return new XS_eqegq_CSS_approx(args);
    }
  }
  return NULL;
}




