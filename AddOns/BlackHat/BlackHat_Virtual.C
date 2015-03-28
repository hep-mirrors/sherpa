#include "AddOns/BlackHat/BlackHat_Virtual.H"

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "blackhat/BH_interface.h"
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
  s_interface->set("alpha_S",AlphaQCD());
  s_interface->set("alpha_QED",AlphaQED());
  BH::BHinput input(moms, sqrt(m_mur2));
  s_interface->operator()(input);

#ifdef VIRTUAL_PREFACTOR
  m_res.Finite() = p_ampl->get_finite()*0.5/VIRTUAL_PREFACTOR;
#else
  m_res.Finite() = p_ampl->get_finite();
#endif
  m_res.IR()     = p_ampl->get_single_pole();
  m_res.IR2()    = p_ampl->get_double_pole();
#ifdef INCLUDE_COUPLINGS_IN_VIRTUAL
  double norm=AlphaQCD()/(2.0*M_PI);
  m_res.Finite()/=norm;
  m_res.IR()/=norm;
  m_res.IR2()/=norm;
  m_born=p_ampl->get_born();
#endif
}

double BlackHat_Virtual::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
  //MSbar
   return 4.*M_PI;
}

double BlackHat_Virtual::ScaleDependenceCoefficient(const int i)
{
  return p_ampl->getScaleVariationCoefficient(i);
}

void BlackHat_Virtual::AddCouplings
(const Process_Info &pi,
 std::vector<std::vector<std::pair<std::string,int> > > &couplings,
 std::vector<std::pair<std::string,int> > cpls,size_t i)
{
  if (i==pi.m_mincpl.size()) {
    couplings.push_back(cpls);
    return;
  }
  for (size_t j(pi.m_mincpl[i]);j<=pi.m_maxcpl[i];++j) {
    cpls[i].second=j;
    if (i==0 && (pi.m_fi.m_nloqcdtype&nlo_type::loop)) ++cpls[i].second;
    if (i==1 && (pi.m_fi.m_nloewtype&nlo_type::loop)) ++cpls[i].second;
    AddCouplings(pi,couplings,cpls,i+1);
  }
}

DECLARE_VIRTUALME2_GETTER(BlackHat_Virtual,"BlackHat_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<Virtual_ME2_Base,Process_Info,BlackHat_Virtual>::
operator()(const Process_Info &pi) const
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
    for (size_t i=0; i<fl.size(); ++i) kfvector.push_back((long int) fl[i]);
    BH::BH_Ampl* ampl=NULL;
    try {
      msg_Info()<<"Trying BlackHat for "<<kfvector<<" ... "<<std::flush;
      std::vector<std::pair<std::string,int> > cpls;
      cpls.push_back(std::pair<std::string,int>("alpha_QCD",0));
      cpls.push_back(std::pair<std::string,int>("alpha_QED",0));
      if (MODEL::s_model->Name()=="HEFT")
	cpls.push_back(std::pair<std::string,int>("YUK2",0));
      std::vector<std::vector<std::pair<std::string,int> > > couplings;
      BlackHat_Virtual::AddCouplings(pi,couplings,cpls);
#ifdef INCLUDE_COUPLINGS_IN_VIRTUAL
      ampl = BlackHat_Virtual::Interface()->new_ampl(kfvector,couplings);
#else
      ampl = BlackHat_Virtual::Interface()->new_ampl(kfvector);
#endif
    } catch (BH::BHerror err) {
      msg_Info()<<"not found."<<std::endl;
      return NULL;
    }
    if (ampl) {
      msg_Info()<<"found."<<std::endl;
      return new BlackHat_Virtual(pi, fl, ampl);
    }
  }
  return NULL;
}
