#include "XS_Selector.H"
#include "XS_4F.H"
#include "XS_QCD.H"
#include "XS_Drell_Yan.H"
#include "Off_Shell_EW.H"
#include "Off_Shell_QCD.H"
#include "Run_Parameter.H"

using namespace EXTRAXS;
using namespace ATOOLS;

XS_Selector::XS_Selector(XS_Base *const owner):
  p_owner(owner),
  m_offshell(false) {}

XS_Base *XS_Selector::GetXS(const size_t nin, const size_t nout,
			    const ATOOLS::Flavour *flavours,
			    const bool seperate_couplings,
			    size_t nqed, size_t nqcd)
{ 
//   std::cout<<"XS_Selector::GetXS nin="<<nin<<" nout="<<nout<<"\n";
//   std::cout<<flavours[0]<<" "<<flavours[1]<<" -> "<<flavours[2]<<" "<<flavours[3]<<"  
//   ("<<seperate_couplings<<") : "<<nqed<<","<<nqcd<<std::endl;
  XS_Base * xs=NULL;
  if (seperate_couplings) {
    for (size_t i=0;i<=2;++i) {
      XS_Base * xst = GetSingleXS(nin,nout,flavours,i,2-i);
      if (xs==0) {
	xs = xst;
	//	if (xs) std::cout<<" new single !!! \n";
      }
      else if (xst!=0) {
	//	std::cout<<" new group !!! \n";
	XS_Group * group = new XS_Group(nin,nout,flavours);
	group->Add(xs);
	group->Add(xst);
	xs=group;
      }
    }
  }
  else {
    nqed=Min(nqed,(size_t)2);
    nqcd=Min(nqcd,(size_t)2);
    for (int j=nqcd;j>=0;--j) {
      for (int i=nqed;i>=0;--i) {
	XS_Base * xst = GetSingleXS(nin,nout,flavours,i,j);
	if (xst!=NULL) return xst;
      }
    }
  }
  //  if (xs!=0) std::cout<<"found"<<std::endl;
  return xs;
}

struct Flavour_Container {
  ATOOLS::Flavour fl[4];
  Flavour_Container(const ATOOLS::Flavour *flav)
  { 
    for (short unsigned int i=0;i<4;++i) fl[i]=flav[i];
  }
};

bool operator<(const Flavour_Container &c1,const Flavour_Container &c2)
{
  for (short unsigned int i=0;i<4;++i) {
    int kf1=c1.fl[i].Kfcode()*(2*c1.fl[i].IsAnti()-1);
    int kf2=c2.fl[i].Kfcode()*(2*c2.fl[i].IsAnti()-1);
    if (kf1<kf2) return true;
    if (kf1>kf2) return false;
  }
  return false;
}

typedef std::map<Flavour_Container,Single_XS::Getter_Function> Getter_Function_Map;

Single_XS *Dummy_Getter(const size_t nin,const size_t nout,
			const ATOOLS::Flavour *flavours, 
			const size_t nqed, const size_t nqcd)
{
  return NULL;
}

Single_XS *XS_Selector::GetSingleXS(const size_t nin,const size_t nout,
				    const ATOOLS::Flavour *flavours,
				    const size_t nqed,const size_t nqcd)
{ 
  Single_XS *xs=NULL;
  static Getter_Function_Map s_gettermap;
  Getter_Function_Map::const_iterator git=s_gettermap.find(Flavour_Container(flavours));
  if (git!=s_gettermap.end()) {
    return git->second(nin,nout,flavours,nqed,nqcd);
  }
  if (m_offshell) { 
    if ((xs=Single_XS::GetProcess<Off_Shell_qqb_llb>(nin,nout,flavours,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<Off_Shell_q1q2b_lnulb>(nin,nout,flavours,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<Off_Shell_q1q2b_q3q4b>(nin,nout,flavours,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<Off_Shell_gg_qqb>(nin,nout,flavours,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<Off_Shell_gg_gg>(nin,nout,flavours,nqed,nqcd))!=NULL);
    else;
    if (xs!=NULL) {
      xs->SetScaleScheme(p_owner->ScaleScheme());
      xs->SetKFactorScheme(p_owner->KFactorScheme());
    }
    return xs;
  }
  if ((xs=Single_XS::GetProcess<XS_ee_ffbar>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_ee_ffbar>;
  else if ((xs=Single_XS::GetProcess<XS_gg_gg>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_gg_gg>;
  else if ((xs=Single_XS::GetProcess<XS_q1g_q1g>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_q1g_q1g>;
  else if ((xs=Single_XS::GetProcess<XS_gg_q1qbar1>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_gg_q1qbar1>;
  else if ((xs=Single_XS::GetProcess<XS_q1qbar1_gg>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_q1qbar1_gg>;
  else if ((xs=Single_XS::GetProcess<XS_q1qbar1_q1qbar1>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_q1qbar1_q1qbar1>;
  else if ((xs=Single_XS::GetProcess<XS_q1qbar1_q2qbar2>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_q1qbar1_q2qbar2>;
  else if ((xs=Single_XS::GetProcess<XS_q1q1_q1q1>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_q1q1_q1q1>;
  else if ((xs=Single_XS::GetProcess<XS_q1q2_q1q2>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_q1q2_q1q2>;
  else if ((xs=Single_XS::GetProcess<XS_f1f1_f1f1>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_f1f1_f1f1>;
  else if ((xs=Single_XS::GetProcess<XS_f1f1b_f1f1b>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_f1f1b_f1f1b>;
  else if ((xs=Single_XS::GetProcess<XS_f1f1b_f2f2b>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_f1f1b_f2f2b>;
  else if ((xs=Single_XS::GetProcess<XS_f1f2_f1f2>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_f1f2_f1f2>;
  else if ((xs=Single_XS::GetProcess<XS_f1f2b_f1f2b>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_f1f2b_f1f2b>;
  else if ((xs=Single_XS::GetProcess<XS_f1f2_f3f4>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_f1f2_f3f4>;
  else if ((xs=Single_XS::GetProcess<XS_f1f2b_f3f4b>(nin,nout,flavours,nqed,nqcd))!=NULL) 
    s_gettermap[Flavour_Container(flavours)]=Single_XS::GetProcess<XS_f1f2b_f3f4b>;
  if (xs!=NULL) {
    xs->SetScaleScheme(p_owner->ScaleScheme());
    xs->SetKFactorScheme(p_owner->KFactorScheme());
    xs->m_order_ew=nqed;
    xs->m_order_strong=nqcd;
  }
  else {
    s_gettermap[Flavour_Container(flavours)]=Dummy_Getter;
  }
  return xs;
}

size_t XS_Selector::FindInGroup(XS_Group *const group,XS_Base *&xs,
				const size_t nin,const size_t nout,
				const ATOOLS::Flavour *fl) 
{
  if (nin==group->NIn() && nout==group->NOut()) {
    XS_Base *xsi=NULL;
    for (size_t i=0;i<group->Size();i++) {
      xsi=dynamic_cast<XS_Base*>((*group)[i]);
      if (xsi->NIn()==nin && xsi->NOut()==nout) {
	size_t pos=i;
	for (size_t j=0;j<nin+nout;j++) {
	  if (xsi->Flavours()[j] != fl[j]) pos=std::string::npos;
	}
	if (pos!=std::string::npos) {
	  xs=xsi;
	  return pos;
	}
      }
    }
  }
  xs=NULL;
  return std::string::npos;
}





