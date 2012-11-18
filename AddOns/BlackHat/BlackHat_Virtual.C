#include "AddOns/BlackHat/BlackHat_Virtual.H"

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "blackhat/BH_error.h"

using namespace BLACKHAT;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

BH::BH_interface *BLACKHAT::BlackHat_Virtual::s_interface=NULL;
MODEL::Model_Base *BLACKHAT::BlackHat_Virtual::s_model=NULL;

BlackHat_Virtual::BlackHat_Virtual(const Process_Info& pi,
				   const Flavour_Vector& flavs,
				   BH::BH_Ampl* ampl) :
  Virtual_ME2_Base(pi, flavs), p_ampl(ampl)
{
}

BlackHat_Virtual::~BlackHat_Virtual()
{
  // if (p_ampl) delete p_ampl;
}

void BlackHat_Virtual::Calc(const Vec4D_Vector& momenta) {
  std::vector<std::vector<double> > moms(momenta.size(), std::vector<double>(4, 0.0));
  for (size_t i=0; i<momenta.size(); ++i) {
    for (size_t j=0; j<4; ++j) {
      moms[i][j]=momenta[i][j];
    }
  }
  BH::BHinput input(moms, sqrt(m_mur2));
  s_interface->operator()(input);

  m_res.Finite() = p_ampl->get_finite();
  m_res.IR()     = p_ampl->get_single_pole();
  m_res.IR2()    = p_ampl->get_double_pole();
}

double BlackHat_Virtual::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
  //MSbar
   return 4.*M_PI;
}

double BlackHat_Virtual::ScaleDependenceCoefficient(const int i)
{
  return p_ampl->getScaleVariationCoefficient(i);
}

DECLARE_VIRTUALME2_GETTER(BlackHat_Virtual_Getter,"BlackHat_Virtual")
Virtual_ME2_Base *BlackHat_Virtual_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="BlackHat" &&
      pi.m_loopgenerator!="WhiteHat") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    if (pi.m_fi.m_sv=="FullColor")
      BlackHat_Virtual::Interface()->set("COLOR_MODE",std::string("full_color"));
    else if (pi.m_fi.m_sv=="LeadingColor")
      BlackHat_Virtual::Interface()->set("COLOR_MODE",std::string("leading_color"));
    else if (pi.m_fi.m_sv=="FullMinusLeadingColor")
      BlackHat_Virtual::Interface()->set("COLOR_MODE",std::string("full_minus_leading_color"));
    else if (pi.m_fi.m_sv!="")
      THROW(critical_error,"Invalid option '"+pi.m_fi.m_sv+"'");
    Flavour_Vector fl=pi.ExtractFlavours();
    std::vector<int> kfvector;
    for (size_t i=0; i<fl.size(); ++i) kfvector.push_back(fl[i].HepEvt());
    BH::BH_Ampl* ampl=NULL;
    try {
      msg_Info()<<"Trying BlackHat for "<<kfvector<<" ... "<<std::flush;
      ampl = BlackHat_Virtual::Interface()->new_ampl(kfvector);
    } catch (BH::BHerror err) {
      msg_Info()<<"not found."<<std::endl;
      return NULL;
    }
    if (ampl) {
      if (ampl->get_order_qed()!=pi.m_oew) return NULL;
      if (ampl->get_order_qcd()+1!=pi.m_oqcd) return NULL;
      msg_Info()<<"found."<<std::endl;
      return new BlackHat_Virtual(pi, fl, ampl);
    }
  }
  return NULL;
}
