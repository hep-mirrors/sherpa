#include "AddOns/BlackHat/BlackHat_Tree.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "blackhat-lib/BH_error.h"

using namespace BLACKHAT;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

BH::BH_interface *BLACKHAT::BlackHat_Tree::s_interface=NULL;
MODEL::Model_Base *BLACKHAT::BlackHat_Tree::s_model=NULL;
namespace BLACKHAT {
#ifdef USING__Threading
  static pthread_mutex_t s_mtx;
#endif
}

BlackHat_Tree::BlackHat_Tree(const Process_Info& pi,
			     const Flavour_Vector& flavs,
			     BH::BH_Ampl* ampl,const int mode) :
  Tree_ME2_Base(pi, flavs), p_ampl(ampl), m_mode(mode)
{
  m_oqcd=ampl->get_order_qcd()+(m_mode?1:0);
  m_oew=ampl->get_order_qed();
#ifdef USING__Threading
  static bool first(true);
  if (first) pthread_mutex_init(&s_mtx,NULL);
  first=false;
#endif
}

BlackHat_Tree::~BlackHat_Tree()
{
  // if (p_ampl) delete p_ampl;
}

void BlackHat_Tree::SetCouplings(const MODEL::Coupling_Map& cpls)
{
  Tree_ME2_Base::SetCouplings(cpls);
  if (p_aqcd) m_asfac=p_aqcd->Default()/s_model->ScalarFunction("alpha_S");
  if (p_aqed) m_afac=p_aqed->Default()/s_model->ScalarFunction("alpha_QED");
}

double BlackHat_Tree::CouplingFactor(const int oqcd,const int oew) const
{
  double fac(1.0);
  if (p_aqcd && oqcd) fac*=pow(m_asfac*p_aqcd->Factor(),oqcd);
  if (p_aqed && oew) fac*=pow(m_afac*p_aqed->Factor(),oew);
  return fac;
}

double BlackHat_Tree::Calc(const Vec4D_Vector& momenta)
{
  std::vector<std::vector<double> > moms
    (momenta.size(), std::vector<double>(4, 0.0));
  for (size_t i=0; i<momenta.size(); ++i) {
    for (size_t j=0; j<4; ++j) {
      moms[i][j]=momenta[i][j];
    }
  }
#ifdef USING__Threading
  pthread_mutex_lock(&s_mtx);
#endif
  BH::BHinput input(moms,-1.0);
  s_interface->operator()(input);
  double res=p_ampl->get_born()*CouplingFactor(m_oqcd,m_oew);
  if (m_mode)
    res*=p_ampl->get_finite()*
      s_model->ScalarFunction("alpha_S")/(2.0*M_PI);
#ifdef USING__Threading
  pthread_mutex_unlock(&s_mtx);
#endif

  return res;
}

DECLARE_TREEME2_GETTER(BlackHat_Tree,"BlackHat_Tree")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,BlackHat_Tree>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="BlackHat" &&
      pi.m_loopgenerator!="WhiteHat") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype==nlo_type::lo ||
      pi.m_fi.m_nloqcdtype==nlo_type::born ||
      pi.m_fi.m_nloqcdtype==nlo_type::real) {
    Flavour_Vector fl=pi.ExtractFlavours();
    std::vector<int> kfvector;
    for (size_t i=0; i<fl.size(); ++i) kfvector.push_back(fl[i].HepEvt());
    int mode=0;
    BH::BH_Ampl* ampl=NULL;
    try {
      msg_Info()<<"Trying BlackHat for "<<kfvector<<" ... "<<std::flush;
      ampl = BlackHat_Tree::Interface()->new_tree_ampl(kfvector);
      msg_Out()<<"Cannot check whether tree process or not with public "
               <<"BlackHat library.\nPlease reenable when suitable BlackHat "
               <<"version is public.\n";
//      if (!ampl->is_born_LO()) {
//	delete ampl;
//	ampl = BlackHat_Tree::Interface()->new_ampl(kfvector);
//	mode=1;
//      }
    } catch (BH::BHerror err) {
      msg_Info()<<"not found."<<std::endl;
      return NULL;
    }
    if (ampl) {
      msg_Info()<<"found."<<std::endl;
      return new BlackHat_Tree(pi, fl, ampl, mode);
    }
  }
  return NULL;
}
