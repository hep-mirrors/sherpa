#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Phys/Flow.H"
#include "MODEL/Main/Model_Base.H"

#include "EXTRA_XS/Main/ME2_Base.H"

#define PropID(i,j) ((1<<i)|(1<<j))

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;

namespace EXTRAXS {
  class XS_PP_ffbar : public ME2_Base {
  private:

    int     m_r;
    double  m_cpl, m_m2;

  public:

    XS_PP_ffbar(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };
}

XS_PP_ffbar::XS_PP_ffbar(const Process_Info& pi, const Flavour_Vector& fl):
  ME2_Base(pi,fl)
{
  m_sintt=2|4;
  m_r=fl[2].IsAnti();
  m_cpl=sqr(4.*M_PI*MODEL::s_model->ScalarFunction(std::string("alpha_QED"),
                                                   rpa->gen.CplScale()))
        *sqr(sqr(m_flavs[2].Charge()))
        *(m_flavs[2].Strong()?3.0:1.0);
  m_m2=sqr(m_flavs[2].Mass());
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_oew=2; m_oqcd=0;
  m_cfls[PropID(0,2)].push_back(fl[2]);
  m_cfls[PropID(1,2)].push_back(fl[2]);
  m_cfls[PropID(0,3)].push_back(fl[3]);
  m_cfls[PropID(1,3)].push_back(fl[3]);
}

double XS_PP_ffbar::operator()(const Vec4D_Vector& mom)
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2+m_r]).Abs2();
  double u=(mom[0]-mom[3-m_r]).Abs2();
  //if (s<m_threshold) return 0.;
  double tp(t-m_m2), up(u-m_m2);
  double mtt(2.0*(tp*(up-2.0*m_m2)-4.0*m_m2*m_m2)/(tp*tp));
  double muu(2.0*(up*(tp-2.0*m_m2)-4.0*m_m2*m_m2)/(up*up));
  double mtu(2.0*m_m2*(s-4.0*m_m2)/(tp*up));
  return m_cpl*CouplingFactor(0,2)*(mtt+muu+2.0*mtu);
}


bool XS_PP_ffbar::SetColours(const Vec4D_Vector& mom)
{
  p_colours[2+m_r][0]=p_colours[3-m_r][1]=Flow::Counter();
  return true;
}

DECLARE_TREEME2_GETTER(XS_PP_ffbar,"XS_PP_ffbar")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_PP_ffbar>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo || pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].IsPhoton() && fl[1].IsPhoton() &&
      fl[2].IsFermion() && fl[2].Charge() &&
      fl[3]==fl[2].Bar()) {
    if (pi.m_oqcd==0 && pi.m_oew==2) {
      return new XS_PP_ffbar(pi,fl);
    }
  }
  return NULL;
}

