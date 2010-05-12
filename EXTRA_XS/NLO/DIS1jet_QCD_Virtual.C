#include "EXTRA_XS/NLO/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {
  class DIS1jet_QCD_Virtual : public Virtual_ME2_Base {
    double m_cpl;
//     ME2_Base* p_tree;
  public:
    DIS1jet_QCD_Virtual(const Process_Info& pi, const Flavour_Vector& flavs) :
      Virtual_ME2_Base(pi, flavs)/*, p_tree(NULL)*/
    {
//       Process_Info tree_pi(pi);
//       tree_pi.m_nloqcdtype=nlo_type::lo;
//       p_tree = ME2_Base::GetME2("DIS1jet", tree_pi);
      m_needsborn=true;

      m_cpl=MODEL::s_model->ScalarFunction(std::string("alpha_S"),
                                           sqr(rpa.gen.Ecms()));
      m_cpl/=2.*M_PI;
      m_cpl*=4.0/3.0;
    }

    ~DIS1jet_QCD_Virtual() {
//       if (p_tree) delete p_tree;
    }

    void Calc(const ATOOLS::Vec4D_Vector& mom);

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
      return 2.*M_PI*m_mur2/(mom[0]*mom[2]);
    }
  };
}


void DIS1jet_QCD_Virtual::Calc(const Vec4D_Vector& mom) {
//   m_born=(*p_tree)(mom);
  m_born*=m_cpl*CouplingFactor(1, 0);
  m_res.IR()=-3.*m_born;
  m_res.IR()=-2.*m_born;
  m_res.Finite()=(-8.)*m_born;
}


DECLARE_VIRTUALME2_GETTER(DIS1jet_QCD_Virtual_Getter,"DIS1jet_QCD_Virtual")
Virtual_ME2_Base *DIS1jet_QCD_Virtual_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype==nlo_type::loop) {
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl.size()!=4) return NULL;
    if (fl[0].IsLepton() && fl[1].IsQuark() && fl[2]==fl[0]  && fl[3]==fl[1]) {
      if ((pi.m_oqcd==1 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
        return new DIS1jet_QCD_Virtual(pi, fl);
      }
    }
  }
  return NULL;
}
