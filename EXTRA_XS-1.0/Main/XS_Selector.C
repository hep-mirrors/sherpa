#include "XS_Selector.H"
#include "XS_EW.H"
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

Single_XS *XS_Selector::GetXS(const size_t nin,const size_t nout,
			      const ATOOLS::Flavour *flavours)
{ 
  Single_XS *xs=NULL;
  if (m_offshell) { 
    if ((xs=Single_XS::GetProcess<Off_Shell_qqb_llb>(nin,nout,flavours))!=NULL);
    else if ((xs=Single_XS::GetProcess<Off_Shell_q1q2b_lnulb>(nin,nout,flavours))!=NULL);
    else if ((xs=Single_XS::GetProcess<Off_Shell_q1q2b_q3q4b>(nin,nout,flavours))!=NULL);
    else if ((xs=Single_XS::GetProcess<Off_Shell_gg_qqb>(nin,nout,flavours))!=NULL);
    else if ((xs=Single_XS::GetProcess<Off_Shell_gg_gg>(nin,nout,flavours))!=NULL);
    else;
    if (xs!=NULL) {
      xs->SetScaleScheme(p_owner->ScaleScheme());
      xs->SetKFactorScheme(p_owner->KFactorScheme());
      xs->SetScaleFactor(p_owner->ScaleFactor());
    }
    return xs;
  }
  if ((xs=Single_XS::GetProcess<XS_q1q2b_q3q4b>(nin,nout,flavours))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_ee_ffbar>(nin,nout,flavours))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_gg_gg>(nin,nout,flavours))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1g_q1g>(nin,nout,flavours))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_gg_q1qbar1>(nin,nout,flavours))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1qbar1_gg>(nin,nout,flavours))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1qbar1_q1qbar1>(nin,nout,flavours))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1qbar1_q2qbar2>(nin,nout,flavours))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1q1_q1q1>(nin,nout,flavours))!=NULL);
  else if ((xs=Single_XS::GetProcess<XS_q1q2_q1q2>(nin,nout,flavours))!=NULL);
  if (xs!=NULL) {
    xs->SetScaleScheme(p_owner->ScaleScheme());
    xs->SetKFactorScheme(p_owner->KFactorScheme());
    xs->SetScaleFactor(p_owner->ScaleFactor());
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





