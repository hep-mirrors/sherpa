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

Single_XS *XS_Selector::GetSingleXS(const size_t nin,const size_t nout,
				    const ATOOLS::Flavour *flavours,
				    const size_t nqed,const size_t nqcd)
{ 
  Single_XS *xs=NULL;
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
      xs->SetScaleFactor(p_owner->ScaleFactor());
    }
    return xs;
  }
  if ((xs=Single_XS::GetProcess<XS_ee_ffbar>(nin,nout,flavours,nqed,nqcd))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_gg_gg>(nin,nout,flavours,nqed,nqcd))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1g_q1g>(nin,nout,flavours,nqed,nqcd))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_gg_q1qbar1>(nin,nout,flavours,nqed,nqcd))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1qbar1_gg>(nin,nout,flavours,nqed,nqcd))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1qbar1_q1qbar1>(nin,nout,flavours,nqed,nqcd))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1qbar1_q2qbar2>(nin,nout,flavours,nqed,nqcd))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1q1_q1q1>(nin,nout,flavours,nqed,nqcd))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1q2_q1q2>(nin,nout,flavours,nqed,nqcd))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_f1f1_f1f1>(nin,nout,flavours,nqed,nqcd))!=NULL);
  if (xs!=NULL) {
    xs->SetScaleScheme(p_owner->ScaleScheme());
    xs->SetKFactorScheme(p_owner->KFactorScheme());
    xs->SetScaleFactor(p_owner->ScaleFactor());
    xs->m_order_ew=nqed;
    xs->m_order_strong=nqcd;
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





