#include "XS_Selector.H"
#include "XS_QCD.H"
#include "XS_Drell_Yan.H"
#include "Run_Parameter.H"

using std::cout;
using std::cerr;
using std::endl;

using namespace EXTRAXS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

Single_XS * XS_Selector::GetXS(int nin,int nout,Flavour * fl)
{ 
  if (rpa.gen.Tracking()) {
    cout<<"looking for: ";
    for (int i=0; i<nin; ++i) cout<<fl[i]<<' ';
    cout<<"-> ";
    for (int i=0; i<nout; ++i) cout<<fl[i+nin]<<' ';
    cout<<endl;
  }

  Single_XS * xs = 0;
  if (nin !=2 && nout !=2) {
    if (AORGTOOLS::rpa.gen.Error()) {
      cout<<"Such a XS is not available as FastFunc!"<<endl;
      cout<<"nin, nout = "<<nin<<", "<<nout<<endl;
    }
    return xs;
  }



  if (fl[2].isfermion() && fl[3]==fl[2].bar() &&
      fl[0].isphoton()  && fl[1]==fl[0])    { return new XS_pp_ffbar(nin,nout,fl); }

  if (fl[2].islepton() && fl[3]==fl[2].bar() &&
      fl[0].isquark()  && fl[1]==fl[0].bar())    { return new XS_ee_ffbar(nin,nout,fl); }
  if (fl[0].islepton() && fl[1]==fl[0].bar() &&
      fl[2].isquark()  && fl[3]==fl[2].bar())    { return new XS_ee_ffbar(nin,nout,fl); }

  if (((fl[0].isquark() && fl[1].isgluon()) ||
       (fl[1].isquark() && fl[0].isgluon()) )   &&
      (((fl[2] == fl[0]) && (fl[3]==fl[1])) ||
       ((fl[3] == fl[0]) && (fl[2]==fl[1])) ) )  { return new XS_q1g_q1g(nin,nout,fl); }
  if (fl[0].isgluon() && fl[1].isgluon()) {
    if (fl[2].isquark() && (fl[3]==fl[2].bar())) { return new XS_gg_q1qbar1(nin,nout,fl); }
    if (fl[2].isgluon() && fl[3].isgluon())      { return new XS_gg_gg(nin,nout,fl); }
  }
  if (fl[0].isquark() && (fl[1]==fl[0].bar())) {
    if (fl[2].isgluon() && fl[3].isgluon())      { return new XS_q1qbar1_gg(nin,nout,fl); }
    if ( ((fl[2]==fl[0]) && (fl[3]==fl[1])) ||
      ((fl[2]==fl[1]) && (fl[2]==fl[1])) )       { return new XS_q1qbar1_q1qbar1(nin,nout,fl); }
    if (fl[2].isquark() && (fl[3]==fl[2].bar())) { return new XS_q1qbar1_q2qbar2(nin,nout,fl); }
  }
  if ( fl[0].isquark() && (fl[1]==fl[0]) &&
	(fl[2]==fl[0]) && (fl[3]==fl[0]) )       { return new XS_q1q1_q1q1(nin,nout,fl); }
  if ( fl[0].isquark() && fl[1].isquark() &&
       ( ((fl[2]==fl[0]) && (fl[3]==fl[1])) ||
	 ((fl[3]==fl[0]) && (fl[2]==fl[1]))) )   { return new XS_q1q2_q1q2(nin,nout,fl); }


  //  if (rpa.gen.Error()) {
    cout<<"looking for: ";
    for (int i=0; i<nin; ++i) cout<<fl[i]<<' ';
    cout<<"-> ";
    for (int i=0; i<nout; ++i) cout<<fl[i+nin]<<' ';
    cout<<endl;
    //  }
  AORGTOOLS::msg.Error()<<"Such a XS is not yet available as FastFunc!"<<endl;
  return 0;
}



bool XS_Selector::FindInGroup(XS_Group * group,XS_Base *& xs,
			      int nin,int nout,Flavour * fl) {
  if ((nin != group->Nin()) || (nout != group->Nout())) {
    xs = 0; 
    return 0; 
  }
  XS_Base * xsi = 0;
  for (int i=0;i<group->Size();i++) {
    xsi = (*group)[i];
    if ( (xsi->Nin() == nin) && (xsi->Nout() == nout)) {
      bool found = 1;
      for (int j=0;j<nin+nout;j++) {
	if (xsi->Flavs()[j] != fl[j]) {
	  found = 0; 
	  break;
	}
      }
      if (found==1) {
	xs = xsi;
	return 1;
      }
    }
  }
  xs = 0;
  return 0;
}





