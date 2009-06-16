#include "EXTRA_XS/NLO/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Exception.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "MODEL/Main/Running_AlphaQED.H"

#define CF 1.
#define CA 1.
#define TR 0.5

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;


namespace EXTRAXS {
  class nunuee_EW_Subtraction : public ME2_Base {
    double        m_born;
    ME2_Base*     p_tree;
    Vec4D_Vector  m_bornmoms;
    Vec4D         Q;
    double        x1,x2,Q2;
    double        m_dalpha;
    double        m_eq2, m_aqed;

  public:
    nunuee_EW_Subtraction(const Process_Info& pi, const Flavour_Vector& flavs,
                   ME2_Base* tree) :
      ME2_Base(pi, flavs), p_tree(tree),
      m_eq2(flavs[3].Charge()),
      m_aqed(MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms()))))

    {
      DEBUG_INFO("Drell Yan EW Subtraction ME opened ...");
      m_bornmoms.reserve(4);
      m_dalpha = 1.;
      double helpd;
      Data_Reader reader(" ",";","!","=");
      if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA")) {
        m_dalpha = helpd;
        msg_Tracking()<<"Set dipole cut alpha="<<m_dalpha<<"."<<std::endl;
      }
// std::cout<<"m_dalpha     =   "<<m_dalpha<<endl;
    }

    ~nunuee_EW_Subtraction() {
      if (p_tree) delete p_tree;
    }

    virtual double operator()(const ATOOLS::Vec4D_Vector& momenta);
    virtual bool SetColours(const ATOOLS::Vec4D_Vector& momenta);
  };
}

double nunuee_EW_Subtraction::operator()(const ATOOLS::Vec4D_Vector& momenta) {

  // yij_k,  hep-ph/9605323 eq(5.4)
  double  m_y423 = momenta[2]*momenta[4]/(momenta[2]*momenta[4]+momenta[4]*momenta[3]+momenta[3]*momenta[2]);
  double  m_y324 = momenta[2]*momenta[3]/(momenta[2]*momenta[4]+momenta[4]*momenta[3]+momenta[3]*momenta[2]);
  m_born=(0.);
  m_bornmoms[0]=momenta[0];
  m_bornmoms[1]=momenta[1];
  Q  = momenta[0]+momenta[1];
  Q2 = Q.Abs2();
  x1 = 2.*momenta[3]*Q/Q2;
  x2 = 2.*momenta[4]*Q/Q2;
  double  zt3 = (momenta[3]*momenta[4])/(momenta[3]*momenta[4] + momenta[2]*momenta[4]);
  double  zt4 = (momenta[3]*momenta[4])/(momenta[3]*momenta[4] + momenta[2]*momenta[3]);
  // 2(p) + 3(q) -> 23(q), spec 4(qb)
  if (m_y324<m_dalpha){
    m_bornmoms[2]=Q-1./x2*momenta[4];
    m_bornmoms[3]=1./x2*momenta[4];
//     m_born+=-(1./(1.-x2)*(2./(2.-x1-x2)-(1.+x1))+(1.-x1)/x2)*(*p_tree)(m_bornmoms);
     m_born+=-(2./(1.-zt3*(1.-m_y324)) - (1.+zt3))*(*p_tree)(m_bornmoms)/(2.*momenta[2]*momenta[3]);
  }
  // 2(p) + 4(qb) -> 23(qb), spec 3(q)
  if (m_y423<m_dalpha){
    m_bornmoms[2]=1./x1*momenta[3];
    m_bornmoms[3]=Q-1./x1*momenta[3];
//     m_born+=-(1./(1.-x1)*(2./(2.-x1-x2)-(1.+x2))+(1.-x2)/x1)*(*p_tree)(m_bornmoms);
     m_born+=-(2./(1.-zt4*(1.-m_y423)) - (1.+zt4))*(*p_tree)(m_bornmoms)/(2.*momenta[2]*momenta[4]);
  }

  m_born*=8.*M_PI*m_aqed*m_eq2*m_eq2;
  return m_born;
}

bool nunuee_EW_Subtraction::SetColours(const ATOOLS::Vec4D_Vector& momenta) {
  return true;
}




DECLARE_ME2_GETTER(nunuee_EW_Subtraction_Getter,"nunuee_EW_Subtraction")
ME2_Base *nunuee_EW_Subtraction_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloewtype==nlo_type::rsub) {
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl.size()!=5) return NULL;
    if ((fl[3].IsLepton() && fl[3]==fl[4].Bar() &&
         fl[0].IsFermion()  && fl[1]==fl[0].Bar() && fl[2].IsPhoton()) ||
        (fl[0].IsLepton() && fl[1]==fl[0].Bar() &&
         fl[3].IsFermion()  && fl[3]==fl[4].Bar() && fl[2].IsPhoton()) &&
        (fl[0]!=fl[3] && fl[0]!=fl[4])) {
      if ((pi.m_oqcd==0 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
        Process_Info lo_pi(pi);
        lo_pi.m_fi.m_ps.erase(lo_pi.m_fi.m_ps.begin());
        lo_pi.m_fi.m_nloewtype=nlo_type::lo;
        ME2_Base* tree_me2 = ME2_Base::GetME2(lo_pi);
        if (!tree_me2) return NULL;
        return new nunuee_EW_Subtraction(pi, fl, tree_me2);
      }
    }
  }
  return NULL;
}
