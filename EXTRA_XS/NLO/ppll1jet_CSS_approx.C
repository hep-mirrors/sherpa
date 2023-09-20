#include "ATOOLS/Math/Tensor.H"
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
// gq -> ll q
//------------------------------------------------------------------------------

namespace EXTRAXS {

  class XS_gqllq_CSS_approx : public ME2_Base {
  private:
    EXTRAXS::ME2_Base * p_bornme;
    bool   m_swap;
    double m_alphasdef;

  public:

    XS_gqllq_CSS_approx(const External_ME_Args& args,
			const int swap);
    ~XS_gqllq_CSS_approx();

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    double LOME2(const Vec4D&, const Vec4D&, const Vec4D&, const Vec4D&,
                 const Vec4D&, int);
    void FillCombinations(std::set<std::pair<size_t,size_t> > &combs,
			  std::map<size_t,ATOOLS::Flavour_Vector> &fls);
  };
}

XS_gqllq_CSS_approx::XS_gqllq_CSS_approx
(const External_ME_Args& args,const int swap) :
  ME2_Base(args), m_swap(swap)
{
  Flavour_Vector inflavs  = args.m_inflavs;
  Flavour_Vector outflavs = args.m_outflavs;
  outflavs.erase(outflavs.end()-1);
  inflavs[m_swap] = inflavs[1-m_swap].Bar();

  External_ME_Args bargs(inflavs, outflavs, {0,2});
  p_bornme = dynamic_cast<ME2_Base*>(PHASIC::Tree_ME2_Base::GetME2(bargs));
  if (!p_bornme) THROW(fatal_error,"no born me found.");
  m_alphasdef = MODEL::as->Default();
  m_oqcd=1;
  m_oew =2;
  PRINT_INFO("initialised XS_gqllq_CSS_approx2");
}

XS_gqllq_CSS_approx::~XS_gqllq_CSS_approx()
{
  if (p_bornme) delete p_bornme;
}

void XS_gqllq_CSS_approx::FillCombinations
(std::set<std::pair<size_t,size_t> > &combs,
 std::map<size_t,ATOOLS::Flavour_Vector> &fls)
{
  combs.insert(std::pair<size_t,size_t>(1<<m_swap,1<<4));
  fls[(1<<m_swap)|(1<<4)].push_back(Flavour(m_flavs[4]));
  combs.insert(std::pair<size_t,size_t>(1<<2,1<<3));
  fls[(1<<2)|(1<<3)].push_back(Flavour(kf_Z));
  fls[(1<<2)|(1<<3)].push_back(Flavour(kf_photon));
  CompleteCombinations(combs,fls);
}

double XS_gqllq_CSS_approx::operator()
(const ATOOLS::Vec4D_Vector& p)
{
  double res(0);
  res+=LOME2(p[m_swap],p[4],p[1-m_swap],p[2],p[3],m_swap);
  return res;
}

double XS_gqllq_CSS_approx::LOME2(const Vec4D& pi, const Vec4D& pj,
                                  const Vec4D& pk, const Vec4D& k1,
                                  const Vec4D& k2, int ij)
{
  DEBUG_FUNC("");
  double pipj(pi*pj), pipk(pi*pk), pjpk(pj*pk);
  double xiab=(pipk-pipj-pjpk)/pipk;

  Vec4D pijt=xiab*pi;
  Vec4D pkt=pk;

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
  double born((*p_bornme)(moms));
  // SF = 8*pi*TR/(2pipj*x) * (1-2x(1-x))
  double split(8.*M_PI/((pi+pj).Abs2()*xiab)*.5*(1.-2.*xiab*(1.-xiab)));
  msg_Debugging()<<8.*M_PI*m_alphasdef<<std::endl;
  msg_Debugging()<<"M2 = "<<born<<" ,  SF = "<<split*m_alphasdef<<std::endl;
  return born*split*m_alphasdef*CouplingFactor(1,0);
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_gqllq_CSS_approx,"XS_gqllq_CSS_approx")
Tree_ME2_Base *ATOOLS::Getter
<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::XS_gqllq_CSS_approx>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  Settings& s = Settings::GetMainSettings();
  if (!s["EXTRAXS_CSS_APPROX_ME"].Get<bool>()) return NULL;
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=5) return NULL;
  if (fl[1].IsQuark()  && fl[4]==fl[1] &&
      fl[0].IsGluon()  &&
      fl[2].IsLepton() && fl[3]==fl[2].Bar()) {
    if (args.m_orders[0]==1 && args.m_orders[1]==2) {
      return new XS_gqllq_CSS_approx(args,0);
    }
  }
  if (fl[0].IsQuark()  && fl[4]==fl[0] &&
      fl[1].IsGluon()  &&
      fl[2].IsLepton() && fl[3]==fl[2].Bar()) {
    if (args.m_orders[0]==1 && args.m_orders[1]==2) {
      return new XS_gqllq_CSS_approx(args,1);
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

    XS_qqllg_CSS_approx(const External_ME_Args& args);
    ~XS_qqllg_CSS_approx();

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    double LOME2(const Vec4D&, const Vec4D&, const Vec4D&, const Vec4D&,
                 const Vec4D&, int);
    void FillCombinations(std::set<std::pair<size_t,size_t> > &combs,
			  std::map<size_t,ATOOLS::Flavour_Vector> &fls);
  };
}

XS_qqllg_CSS_approx::XS_qqllg_CSS_approx
(const External_ME_Args& args ) : ME2_Base(args)
{
  Flavour_Vector inflavs  = args.m_inflavs;
  Flavour_Vector outflavs = args.m_outflavs;
  outflavs.erase(outflavs.end()-1);
  External_ME_Args bargs(inflavs,outflavs,{0,2});
  p_bornme = dynamic_cast<ME2_Base*>(PHASIC::Tree_ME2_Base::GetME2(bargs));
  if (p_bornme) THROW(fatal_error,"no born me found.");
  m_alphasdef = MODEL::as->Default();
  m_oqcd=1;
  m_oew =2;
  PRINT_INFO("initialised XS_qqllg_CSS_approx2");
}

XS_qqllg_CSS_approx::~XS_qqllg_CSS_approx()
{
  if (p_bornme) delete p_bornme;
}

void XS_qqllg_CSS_approx::FillCombinations
(std::set<std::pair<size_t,size_t> > &combs,
 std::map<size_t,ATOOLS::Flavour_Vector> &fls)
{
  combs.insert(std::pair<size_t,size_t>(1<<0,1<<4));
  fls[(1<<0)|(1<<4)].push_back(Flavour(m_flavs[0].Bar()));
  combs.insert(std::pair<size_t,size_t>(1<<1,1<<4));
  fls[(1<<1)|(1<<4)].push_back(Flavour(m_flavs[1].Bar()));
  combs.insert(std::pair<size_t,size_t>(1<<2,1<<3));
  fls[(1<<2)|(1<<3)].push_back(Flavour(kf_Z));
  fls[(1<<2)|(1<<3)].push_back(Flavour(kf_photon));
  CompleteCombinations(combs,fls);
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
  double pipj(pi*pj), pipk(pi*pk), pjpk(pj*pk);
  double xiab=(pipk-pipj-pjpk)/pipk;

  Vec4D pijt=xiab*pi;
  Vec4D pkt=pk;

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
  double born((*p_bornme)(moms));
  // SF = 8*pi*CF/(2pipj*x) * (2/(1-x)-(1+x))
  double split(8.*M_PI/((pi+pj).Abs2()*xiab)*4./3.*(2./(1.-xiab)-(1.+xiab)));
  msg_Debugging()<<8.*M_PI*m_alphasdef<<std::endl;
  msg_Debugging()<<"M2 = "<<born<<" ,  SF = "<<split*m_alphasdef<<std::endl;
  return born*split*m_alphasdef*CouplingFactor(1,0);
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_qqllg_CSS_approx,"XS_qqllg_CSS_approx")
Tree_ME2_Base *ATOOLS::Getter
<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::XS_qqllg_CSS_approx>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  Settings& s = Settings::GetMainSettings();
  if (!s["EXTRAXS_CSS_APPROX_ME"].Get<bool>()) return NULL;
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=5) return NULL;
  if (fl[0].IsQuark()  && fl[1]==fl[0].Bar() &&
      fl[4].IsGluon()  &&
      fl[2].IsLepton() && fl[3]==fl[2].Bar()) {
    if (args.m_orders[0]==1 && args.m_orders[1]==2) {
      return new XS_qqllg_CSS_approx(args);
    }
  }
  return NULL;
}



