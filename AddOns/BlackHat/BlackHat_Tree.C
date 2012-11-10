#include "AddOns/BlackHat/BlackHat_Tree.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "blackhat/BH_error.h"

using namespace WHITEHAT;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

BH::BH_interface *WHITEHAT::BlackHat_Tree::s_interface=NULL;
MODEL::Model_Base *WHITEHAT::BlackHat_Tree::s_model=NULL;
namespace WHITEHAT {
#ifdef USING__Threading
  static pthread_mutex_t s_mtx;
#endif
}

BlackHat_Tree::BlackHat_Tree(const Process_Info& pi,
			     const Flavour_Vector& flavs,
			     BH::BH_Ampl* ampl) :
  Tree_ME2_Base(pi, flavs), p_ampl(ampl),
  p_aqcd(NULL), p_aqed(NULL)
{
  size_t nqcd(0);
  for (size_t i(0);i<flavs.size();++i)
    if (flavs[i].Strong()) ++nqcd;
  m_oqcd=nqcd-2;
  m_oew=flavs.size()-m_oqcd-2;
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

void BlackHat_Tree::SetCouplings(MODEL::Coupling_Map *const cpls)
{
  if (cpls->find("Alpha_QCD")!=cpls->end()) {
    p_aqcd=cpls->Get("Alpha_QCD");
    m_asfac=p_aqcd->Default()/s_model->ScalarFunction("alpha_S");
  }
  if (cpls->find("Alpha_QED")!=cpls->end()) {
    p_aqed=cpls->Get("Alpha_QED");
    m_afac=p_aqed->Default()/s_model->ScalarFunction("alpha_QED");
  }
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
#ifdef USING__Threading
  pthread_mutex_unlock(&s_mtx);
#endif

  return res;
}

DECLARE_TREEME2_GETTER(BlackHat_Tree_Getter,"BlackHat_Tree")
Tree_ME2_Base *BlackHat_Tree_Getter::operator()(const Process_Info &pi) const
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
    BH::BH_Ampl* ampl=NULL;
    try {
      msg_Info()<<"Trying BlackHat for "<<kfvector<<" ... "<<std::flush;
      ampl = BlackHat_Tree::Interface()->new_tree_ampl(kfvector);
    } catch (BH::BHerror err) {
      msg_Info()<<"not found."<<std::endl;
      return NULL;
    }
    if (ampl) {
      msg_Info()<<"found."<<std::endl;
      return new BlackHat_Tree(pi, fl, ampl);
    }
  }
  return NULL;
}
