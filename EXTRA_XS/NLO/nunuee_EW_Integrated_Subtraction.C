#include "EXTRA_XS/NLO/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Exception.H"
#include "MODEL/Main/Running_AlphaQED.H"

#define CF 1.33333333333333333
#define CA 3.
#define TR 0.5

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {
  class nunuee_EW_Integrated_Subtraction : public Virtual_ME2_Base {
    bool   m_isdy;
    ME2_Base* p_tree;
    double m_eq2, m_aqed;
    double m_alpha, m_logalpha;
    double m_dalpha;

  public:
    nunuee_EW_Integrated_Subtraction(const Process_Info& pi, const Flavour_Vector& flavs,
                   ME2_Base* tree) :
      Virtual_ME2_Base(pi, flavs), p_tree(tree),
      m_eq2(flavs[3].Charge()),
      m_aqed(MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms()))))

    {
      m_dalpha = 1.;
      double helpd;
      Data_Reader reader(" ",";","!","=");
      if (reader.ReadFromFile(helpd,"DIPOLE_ALPHA")) {
        m_dalpha = helpd;
        msg_Tracking()<<"Set dipole cut alpha="<<m_dalpha<<"."<<std::endl;
      }
      m_isdy = true;
/*      if (flavs[2].IsLepton() && flavs[3]==flavs[2].Bar() &&
          flavs[0].IsQuark()  && flavs[1]==flavs[0].Bar()) m_isdy = true;
      else m_isdy = false;*/
      m_alpha = 1.;
//       m_logalpha = log(m_alpha);
    }

    ~nunuee_EW_Integrated_Subtraction() {
      if (p_tree) delete p_tree;
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    virtual bool SetColours(const ATOOLS::Vec4D_Vector& momenta);
  };
}

bool nunuee_EW_Integrated_Subtraction::SetColours(const ATOOLS::Vec4D_Vector& momenta) {
  return true;
}

void nunuee_EW_Integrated_Subtraction::Calc(const Vec4D_Vector& momenta) {

  m_born=(*p_tree)(momenta);
  m_born*=m_aqed*m_eq2*m_eq2/(2.*M_PI);

  if (m_isdy == true) {
    // eps_IR^1
    m_res.IR()=3.*m_born;
    // eps_IR^2
    m_res.IR2()=2.*m_born;
    // finite
    m_res.Finite()=(10.-sqr(M_PI))*m_born;
  }
  else {
    // eps_IR^1
    m_res.IR()=3.*m_born;
    // eps_IR^2
    m_res.IR2()=2.*m_born;
    // finite
    m_res.Finite()=(7.-sqr(M_PI)+3.*m_alpha-3.*m_logalpha-2.*m_logalpha*m_logalpha)*m_born;
  }
}


DECLARE_VIRTUALME2_GETTER(nunuee_EW_Integrated_Subtraction_Getter,"nunuee_EW_Integrated_Subtraction")
Virtual_ME2_Base *nunuee_EW_Integrated_Subtraction_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloewtype&nlo_type::vsub) {
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl.size()!=4) return NULL;
    if ((fl[2].IsFermion() && fl[3]==fl[2].Bar() &&
         fl[0].IsFermion() && fl[1]==fl[0].Bar()) &&
        (fl[0].Electromagnetic() || fl[2].Electromagnetic()) &&
        (fl[0]!=fl[2]  &&  fl[0]!=fl[3])) {
      if ((pi.m_oqcd==0 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
        Process_Info tree_pi(pi);
        tree_pi.m_fi.m_nloewtype=nlo_type::lo;
        ME2_Base* tree_me2 = ME2_Base::GetME2(tree_pi);
        if (!tree_me2) return NULL;
        return new nunuee_EW_Integrated_Subtraction(pi, fl, tree_me2);
      }
    }
  }
  return NULL;
}
