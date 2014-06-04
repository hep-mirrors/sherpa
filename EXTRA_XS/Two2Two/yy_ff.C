#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

#include "EXTRA_XS/Main/ME2_Base.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


namespace EXTRAXS {

  class yy_ff : public ME2_Base {
  private:
    double m_mass;
    int m_charge;
  public:

    yy_ff(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };

  yy_ff::yy_ff(const Process_Info& pi, const Flavour_Vector& fl)
    : ME2_Base(pi, fl)
  {
    m_sintt=1;
    m_oew=0;
    m_oqcd=0;

    m_mass = fl[2].HadMass();
    m_charge = fl[2].Charge();
    if (m_charge == 0) THROW(fatal_error, "You should change a known particle's mass to zero!");
  }

  double yy_ff::operator()(const ATOOLS::Vec4D_Vector& momenta)
  {
    // This is a dummy matrix element to simulate the production of charged
    // boson pairs in EPA
    //std::cout << "HEEEEEEEEEEEEEELLLLLLLLLLLLLLLLLLLLLLOOOOOOOOOO" << std::endl;
    double p2((momenta[0]+momenta[1]).Abs2());
    //std::cout << 4*m_mass*m_mass << "\t" << p2 << std::endl;
    double m2 = m_mass*m_mass;
    if ((4*m2) > p2) return 0.;

    // check whether two particle production threshold is reached
    double xs = 0;
    xs=2 + 8*m2/p2 - 16*m2*m2/(p2*p2);
    xs*=2*log(sqrt(p2)/(2*m_mass) + sqrt(p2/(4*m2) - 1));
    xs-=sqrt(1-4*m2/p2)*(2+8*m2/p2);
    xs *= 4*M_PI*0.0072992701*0.0072992701*pow(m_charge, 4.)/p2;
    std::cout << "xs=" << xs << "\t" << "p2=" << p2 <<std::endl;
    return xs;
  }
}

DECLARE_TREEME2_GETTER(yy_ff,"yy_ff")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,yy_ff>::
operator()(const Process_Info &pi) const
{
  return NULL;
  if (pi.m_fi.NLOType()!=nlo_type::lo && pi.m_fi.NLOType()!=nlo_type::born)
    return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0]==Flavour(kf_photon) && fl[1]==Flavour(kf_photon) &&
      (fl[2].Kfcode()==kf_mu || fl[2].Kfcode()==kf_tau) &&
      fl[3]==fl[2].Bar())
  {
    return new yy_ff(pi, fl);
  }
  return NULL;
}
