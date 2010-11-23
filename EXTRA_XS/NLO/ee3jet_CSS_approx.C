#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Org/Data_Reader.H"

#define CF 1.33333333333333333

using namespace ATOOLS;
using namespace EXTRAXS;
using namespace PHASIC;
using namespace std;

namespace EXTRAXS {

  class XS_ee3jet_CSS_approx : public ME2_Base {
  private:
    EXTRAXS::ME2_Base * p_bornme;
    double m_alphasdef;

  public:

    XS_ee3jet_CSS_approx(const Process_Info& pi,const Flavour_Vector& fl);
    ~XS_ee3jet_CSS_approx();

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    double LOME2(const Vec4D&, const Vec4D&, const Vec4D&, const Vec4D&,
		 const Vec4D&, int);
  };
}

XS_ee3jet_CSS_approx::XS_ee3jet_CSS_approx
(const Process_Info& pi, const Flavour_Vector& fl) : ME2_Base(pi, fl)
{
  PRINT_INFO("initialised XS_ee3jet_CSS_approx");
  Process_Info pico(pi);
  pico.m_fi.m_ps.erase(pico.m_fi.m_ps.begin());
  pico.m_fi.m_nloqcdtype=nlo_type::lo;
  PRINT_INFO(pico);
  p_bornme = dynamic_cast<ME2_Base*>(PHASIC::Tree_ME2_Base::GetME2(pico));
  m_alphasdef = (*MODEL::as)(rpa.gen.CplScale());
  PRINT_INFO("initialised XS_ee3jet_CSS_approx2");
}

XS_ee3jet_CSS_approx::~XS_ee3jet_CSS_approx()
{
  if (p_bornme) delete p_bornme;
}

double XS_ee3jet_CSS_approx::operator()
(const ATOOLS::Vec4D_Vector& p)
{
  double res(0);
  res+=LOME2(p[0],p[1],p[3],p[2],p[4],2);
  res+=LOME2(p[0],p[1],p[4],p[2],p[3],3);
  return res;
}

double XS_ee3jet_CSS_approx::LOME2(const Vec4D& p0, const Vec4D& p1,
				   const Vec4D& pi, const Vec4D& pj, 
				   const Vec4D& pk, int ij)
{
  Vec4D Q(pi+pj+pk);
  double Q2=Q*Q, mij2=0., mk2=0.;
  double lrat(Q2/(Q2-(pi+pj).Abs2()));
  Vec4D pkt(lrat*(pk-(Q*pk/Q2)*Q)+(Q2+mk2-mij2)/(2.0*Q2)*Q);
  Vec4D pijt(Q-pkt);
  Vec4D_Vector moms(4);
  moms[0]=p0;
  moms[1]=p1;
  moms[ij]=pijt;
  moms[3-ij+2]=pkt;
  double born((*p_bornme)(moms));
  double zi=(pi*pk)/(pi*pk+pj*pk), yijk=(pi*pj)/(pi*pj+pi*pk+pj*pk);
  double split(8.*M_PI*CF/((pi+pj).Abs2())*(2./(1.-zi+yijk*zi)-(1.+zi)));
  return born*split*m_alphasdef;
}

DECLARE_TREEME2_GETTER(XS_ee3jet_CSS_approx_Getter,
                   "XS_ee3jet_CSS_approx")
Tree_ME2_Base *XS_ee3jet_CSS_approx_Getter::operator()
(const Process_Info &pi) const
{
  Data_Reader read(" ",";","!","=");
  if (read.GetValue<int>("EXTRAXS_CSS_APPROX_ME",0)==0) return NULL;
  if (pi.m_fi.NLOType()!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=5) return NULL;
  if (fl[0].IsLepton() && fl[1]==fl[0].Bar() &&
      fl[2].IsGluon()  &&
      fl[3].IsQuark()  && fl[3]==fl[4].Bar()) {
    if ((pi.m_oqcd==1 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
      return new XS_ee3jet_CSS_approx(pi,fl);
    }
  }
  return NULL;
}


