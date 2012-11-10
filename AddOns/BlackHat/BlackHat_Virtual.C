#include "AddOns/BlackHat/BlackHat_Virtual.H"

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "blackhat/BH_error.h"

using namespace WHITEHAT;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

BH::BH_interface *WHITEHAT::BlackHat_Virtual::s_interface=NULL;
MODEL::Model_Base *WHITEHAT::BlackHat_Virtual::s_model=NULL;

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

void BlackHat_Virtual::InitInterface(Model_Base *model)
{
  BlackHat_Virtual::s_model=model;
  if (s_interface==NULL) {
    msg_Info()<<"Initialising BlackHat interface {"<<std::endl;
    Data_Reader reader(" ",";","!","=");
    s_interface=new BH::BH_interface
      (reader.GetValue<std::string>("BH_SETTINGS_FILE",std::string("")));
    s_interface->set("Z_mass",Flavour(kf_Z).Mass());
    s_interface->set("Z_width",Flavour(kf_Z).Width());
    s_interface->set("W_mass",Flavour(kf_Wplus).Mass());
    s_interface->set("W_width",Flavour(kf_Wplus).Width());
    double sin_th_2=model->ScalarConstant(std::string("sin2_thetaW"));
    s_interface->set("sin_th_2",sin_th_2);
    s_interface->set("sin_2th",sin(2.*asin(sqrt(sin_th_2))));
    s_interface->set("alpha_S",model->ScalarFunction(std::string("alpha_S")));
    s_interface->set("alpha_QED",model->ScalarFunction(std::string("alpha_QED")));
    msg_Info()<<"}"<<std::endl;
  }
}

void BlackHat_Virtual::DeleteInterface()
{
  if (s_interface) delete s_interface;
  s_interface=NULL;
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
      msg_Info()<<"found."<<std::endl;
      return new BlackHat_Virtual(pi, fl, ampl);
    }
  }
  return NULL;
}
