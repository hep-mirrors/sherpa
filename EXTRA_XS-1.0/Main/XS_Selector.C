#include "XS_Selector.H"
#include "XS_QCD.H"
#include "XS_Drell_Yan.H"
#include "Off_Shell_EW.H"
#include "Off_Shell_QCD.H"
#include "Run_Parameter.H"

using std::cout;
using std::cerr;
using std::endl;

using namespace EXTRAXS;
using namespace ATOOLS;

XS_Selector::XS_Selector(XS_Base *const owner):
  p_owner(owner) {}

Single_XS *XS_Selector::GetXS(const size_t nin,const size_t nout,
			      const ATOOLS::Flavour *flavours,const bool offshell)
{ 
  if (offshell) {
    if ((flavours[2].IsLepton() && flavours[3]==flavours[2].Bar() && 
	 flavours[0].IsQuark() && flavours[1]==flavours[0].Bar()) ||
	(flavours[0].IsLepton() && flavours[1]==flavours[0].Bar() && 
	 flavours[2].IsQuark() && flavours[3]==flavours[2].Bar())){ 
      return new Off_Shell_qqb_llb(nin,nout,flavours,p_owner->ScaleScheme(),
				   p_owner->KFactorScheme(),p_owner->ScaleFactor()); 
    }
    if ((flavours[2].IsUptype() && flavours[2].IntCharge()==0 && flavours[3].IsDowntype() && 
	 flavours[0].IsUptype() && flavours[1].IsDowntype()) ||
	(flavours[3].IsUptype() && flavours[3].IntCharge()==0 && flavours[2].IsDowntype() && 
	 flavours[1].IsUptype() && flavours[0].IsDowntype())){ 
      return new Off_Shell_q1q2b_lnulb(nin,nout,flavours,p_owner->ScaleScheme(),
				       p_owner->KFactorScheme(),p_owner->ScaleFactor()); 
    }
    if (flavours[2].IsQuark() && flavours[3]==flavours[2].Bar() && 
	flavours[0].IsGluon() && flavours[1].IsGluon()){ 
      return new Off_Shell_gg_qqb(nin,nout,flavours,p_owner->ScaleScheme(),
				  p_owner->KFactorScheme(),p_owner->ScaleFactor()); 
    }
    ATOOLS::msg.Error()<<"XS_Selector::GetXS("<<nin<<","<<nout<<",["
		       <<flavours[0]<<","<<flavours[1]<<","
		       <<flavours[2]<<","<<flavours[3]<<"]):"<<std::endl
		       <<"   The corresponding process is not impelemented yet ! Abort."<<std::endl;
    exit(174);
  }
  if (flavours[2].IsFermion() && flavours[3]==flavours[2].Bar() &&
      flavours[0].IsPhoton()  && flavours[1]==flavours[0]) { 
    return new XS_pp_ffbar(nin,nout,flavours); 
  }
  if ((flavours[2].IsLepton() && flavours[3]==flavours[2].Bar() && flavours[0].IsQuark() && 
       flavours[1]==flavours[0].Bar()) ||   
      (flavours[0].IsLepton() && flavours[1]==flavours[0].Bar() && flavours[2].IsQuark() && 
       flavours[3]==flavours[2].Bar())) { 
    return new XS_ee_ffbar(nin,nout,flavours); 
  }
  if (((flavours[0].IsQuark() && flavours[1].IsGluon()) ||
       (flavours[1].IsQuark() && flavours[0].IsGluon()) )   &&
      (((flavours[2] == flavours[0]) && (flavours[3]==flavours[1])) ||
       ((flavours[3] == flavours[0]) && (flavours[2]==flavours[1])) ) )  { 
    return new XS_q1g_q1g(nin,nout,flavours); 
  }
  if (flavours[0].IsGluon() && flavours[1].IsGluon()) {
    if (flavours[2].IsQuark() && (flavours[3]==flavours[2].Bar())) { 
      return new XS_gg_q1qbar1(nin,nout,flavours); 
    }
    if (flavours[2].IsGluon() && flavours[3].IsGluon()) { 
      return new XS_gg_gg(nin,nout,flavours); 
    }
  }
  if (flavours[0].IsQuark() && (flavours[1]==flavours[0].Bar())) {
    if (flavours[2].IsGluon() && flavours[3].IsGluon()) { 
      return new XS_q1qbar1_gg(nin,nout,flavours); 
    }
    if ( ((flavours[2]==flavours[0]) && (flavours[3]==flavours[1])) ||
      ((flavours[2]==flavours[1]) && (flavours[2]==flavours[1])) ) { 
      return new XS_q1qbar1_q1qbar1(nin,nout,flavours); 
    }
    if (flavours[2].IsQuark() && (flavours[3]==flavours[2].Bar())) { 
      return new XS_q1qbar1_q2qbar2(nin,nout,flavours); 
    }
  }
  if ( flavours[0].IsQuark() && (flavours[1]==flavours[0]) &&
	(flavours[2]==flavours[0]) && (flavours[3]==flavours[0]) ) { 
    return new XS_q1q1_q1q1(nin,nout,flavours); 
  }
  if ( flavours[0].IsQuark() && flavours[1].IsQuark() &&
       ( ((flavours[2]==flavours[0]) && (flavours[3]==flavours[1])) ||
	 ((flavours[3]==flavours[0]) && (flavours[2]==flavours[1]))) ) { 
    return new XS_q1q2_q1q2(nin,nout,flavours); 
  }

  return 0;
    ATOOLS::msg.Error()<<"XS_Selector::GetXS("<<nin<<","<<nout<<",["
		       <<flavours[0]<<","<<flavours[1]<<","
		       <<flavours[2]<<","<<flavours[3]<<"]):"<<std::endl
		       <<"   The corresponding process is not impelemented yet ! Abort."<<std::endl;
    exit(174);
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
	if (pos==std::string::npos) {
	  xs=xsi;
	  return pos;
	}
      }
    }
  }
  xs=NULL;
  return std::string::npos;
}





