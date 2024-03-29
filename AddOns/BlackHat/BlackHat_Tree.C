#include "PHASIC++/Process/External_ME_Args.H"
#include "AddOns/BlackHat/BlackHat_Tree.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "blackhat/BH_interface.h"
#include "blackhat/BH_error.h"

using namespace BLACKHAT;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

BH::BH_interface *BLACKHAT::BlackHat_Tree::s_interface=NULL;
MODEL::Model_Base *BLACKHAT::BlackHat_Tree::s_model=NULL;
namespace BLACKHAT {
}

BlackHat_Tree::BlackHat_Tree(const External_ME_Args& args,
			     BH::BH_Ampl* ampl,const int mode) :
  Tree_ME2_Base(args), p_ampl(ampl), m_mode(mode)
{
  m_oqcd=ampl->get_order_qcd()+(m_mode?2:0);
  m_oew=ampl->get_order_qed();
  WarnForMassiveFlavours(args.Flavours());
}

BlackHat_Tree::~BlackHat_Tree()
{
  // if (p_ampl) delete p_ampl;
}

void BlackHat_Tree::WarnForMassiveFlavours(const Flavour_Vector& flavs) const
{
  for (Flavour_Vector::const_iterator it(flavs.begin());
      it != flavs.end();
      ++it) {
    if (it->IsQuark() && it->IsMassive()) {
      msg_Error() << "WARNING: BlackHat does not support massive quarks."
        << " Will continue and hope for the best." << std::endl;
      break;
    }
  }
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
  s_interface->set("alpha_S",AlphaQCD());
  s_interface->set("alpha_QED",AlphaQED());
  BH::BHinput input(moms,-1.0);
  s_interface->operator()(input);
  double res=p_ampl->get_born();
  if (m_mode) {
    res*=p_ampl->get_finite();
#ifndef INCLUDE_COUPLINGS_IN_VIRTUAL
    res*=2.0*sqr(AlphaQCD()/(4.0*M_PI));
#endif
  }
  return res;
}

int BlackHat_Tree::OrderQCD(const int &id) const
{
  return m_oqcd;
}

int BlackHat_Tree::OrderEW(const int &id) const
{
  return m_oew;
}

void BlackHat_Tree::AddCouplings
(const External_ME_Args &args,
 std::vector<std::vector<std::pair<std::string,int> > > &couplings,
 std::vector<std::pair<std::string,int> > cpls,size_t i)
{
  if (i==args.m_orders.size()) {
    couplings.push_back(cpls);
    return;
  }
  cpls[i].second=args.m_orders[i];
  AddCouplings(args,couplings,cpls,i+1);
}

DECLARE_TREEME2_GETTER(BLACKHAT::BlackHat_Tree,"BlackHat_Tree")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,BLACKHAT::BlackHat_Tree>::
operator()(const PHASIC::External_ME_Args &args) const
{
  if (args.m_source!="" && args.m_source!="BlackHat") return NULL;
  const Flavour_Vector fl = args.Flavours();
  std::vector<int> kfvector;
  for (size_t i=0; i<fl.size(); ++i) kfvector.push_back((long int) fl[i]);
  int mode=0;
  BH::BH_Ampl* ampl=NULL;
  try {
    msg_Info()<<"Trying BlackHat for "<<kfvector<<" ... "<<std::flush;
    std::vector<std::pair<std::string,int> > cpls;
    cpls.push_back(std::pair<std::string,int>("alpha_QCD",0));
    cpls.push_back(std::pair<std::string,int>("alpha_QED",0));
    if (BLACKHAT::BlackHat_Tree::s_model->Name()=="HEFT")
      cpls.push_back(std::pair<std::string,int>("YUK2",0));
    std::vector<std::vector<std::pair<std::string,int> > > couplings;
    BlackHat_Tree::AddCouplings(args,couplings,cpls);
#ifdef INCLUDE_COUPLINGS_IN_VIRTUAL
    ampl = BlackHat_Tree::Interface()->new_tree_ampl(kfvector,couplings);
#else
    ampl = BlackHat_Tree::Interface()->new_tree_ampl(kfvector);
#endif
#ifdef VIRTUAL_PREFACTOR
    if (ampl && !ampl->is_born_LO()) {
      delete ampl;
#ifdef INCLUDE_COUPLINGS_IN_VIRTUAL
      ampl = BlackHat_Tree::Interface()->new_ampl(kfvector,couplings);
#else
      ampl = BlackHat_Tree::Interface()->new_ampl(kfvector);
#endif
      mode=1;
    }
#else
    msg_Out()<<"Cannot check LO process type with current BlackHat library.\n"
	     <<"Please retry with newer version.\n";
#endif
  } catch (BH::BHerror err) {
    msg_Info()<<"not found."<<std::endl;
    return NULL;
  }
  if (ampl) {
    msg_Info()<<"found."<<std::endl;
    return new BlackHat_Tree(args, ampl, mode);
  }
  return NULL;
}
