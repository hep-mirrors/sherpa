#include "XS_Selector.H"
#include "XS_4F.H"
#include "XS_QCD.H"
#include "XS_SQCD.H"
#include "XS_Drell_Yan.H"
#include "Run_Parameter.H"

using namespace EXTRAXS;
using namespace ATOOLS;

XS_Selector::XS_Selector(XS_Base *const owner):
  p_owner(owner) {}

XS_Base *XS_Selector::GetXS(const size_t nin, const size_t nout,
			    const ATOOLS::Flavour *flavours,
			    const bool seperate_couplings,
			    size_t nqed, size_t nqcd,
			    const bool &sort)
{ 
  msg_Debugging()<<METHOD<<"(): '"<<flavours[0].IDName()
		 <<" "<<flavours[1].IDName()<<" ->";
  for (size_t i(2);i<nin+nout;++i)
    msg_Debugging()<<" "<<flavours[i].IDName();
  msg_Debugging()<<"', sort = "<<sort<<", nqed = "
		 <<nqed<<", nqcd = "<<nqcd<<" => ";
  XS_Base::SetSortFlavours(sort);
  XS_Base * xs=NULL;
  if (seperate_couplings) {
    for (size_t i=0;i<=nout;++i) {
      XS_Base * xst = GetSingleXS(nin,nout,flavours,i,nout-i);
      if (xs==0) {
	xs = xst;
	//	if (xs) std::cout<<" new single !!! \n";
      }
      else if (xst!=0) {
	//	std::cout<<" new group !!! \n";
	XS_Group * group = new XS_Group(nin,nout,flavours,NULL);
	group->Add(xs);
	group->Add(xst);
	xs=group;
      }
    }
  }
  else {
    nqed=Min(nqed,(size_t)nin+nout-2);
    nqcd=Min(nqcd,(size_t)nin+nout-2);
    for (int j=nqcd;j>=0;--j) {
      for (int i=nqed;i>=0;--i) {
	XS_Base * xst = GetSingleXS(nin,nout,flavours,i,j);
	if (xst!=NULL) {
	  msg_Debugging()<<"'"<<xst->Name()<<"'\n";
	  return xst;
	}
      }
    }
  }
  msg_Debugging()<<"'"<<(xs==NULL?"NULL":xs->Name())<<"'\n";
  return xs;
}

struct Flavour_Container {
  size_t nqed, nqcd;
  ATOOLS::Flavour fl[4];
  Flavour_Container(const ATOOLS::Flavour *flav,const size_t ne,const size_t nc) :
    nqed(ne), nqcd(nc)
  { 
    for (short unsigned int i=0;i<4;++i) fl[i]=flav[i];
  }
};

bool operator<(const Flavour_Container &c1,const Flavour_Container &c2)
{
  if (c1.nqed<c2.nqed) return true;
  if (c1.nqed>c2.nqed) return false;
  if (c1.nqcd<c2.nqcd) return true;
  if (c1.nqcd>c2.nqcd) return false;

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
  if (nin+nout==4) {
    Getter_Function_Map::const_iterator git=
      s_gettermap.find(Flavour_Container(flavours,nqed,nqcd));
    if (git!=s_gettermap.end()) {
      xs=git->second(nin,nout,flavours,nqed,nqcd);
      if (xs==NULL) return xs;
      xs->m_orderEW=nqed;
      xs->m_orderQCD=nqcd;
      return xs;
    }
  }
  Flavour_Container flc(flavours,nqed,nqcd);
  XS_Model_Base *model(p_owner->GetModel());
    if ((xs=Single_XS::GetProcess<XS_ee_ffbar>(nin,nout,flavours,nqed,nqcd))!=NULL) 
      s_gettermap[flc]=&Single_XS::GetProcess<XS_ee_ffbar>;
    else if ((xs=Single_XS::GetProcess<XS_pp_q1qbar1>
	      (nin,nout,flavours,model,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<XS_gg_gg>
	      (nin,nout,flavours,model,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<XS_q1g_q1g>
	      (nin,nout,flavours,model,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<XS_gg_q1qbar1>
	      (nin,nout,flavours,model,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<XS_q1qbar1_gg>
	      (nin,nout,flavours,model,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<XS_q1qbar1_q1qbar1>
	      (nin,nout,flavours,model,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<XS_q1qbar1_q2qbar2>
	      (nin,nout,flavours,model,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<XS_q1q1_q1q1>
	      (nin,nout,flavours,model,nqed,nqcd))!=NULL); 
    else if ((xs=Single_XS::GetProcess<XS_q1q2_q1q2>
	      (nin,nout,flavours,model,nqed,nqcd))!=NULL);
    else if ((xs=Single_XS::GetProcess<XS_f1f1_f1f1>(nin,nout,flavours,nqed,nqcd))!=NULL) 
      s_gettermap[flc]=&Single_XS::GetProcess<XS_f1f1_f1f1>;
    else if ((xs=Single_XS::GetProcess<XS_f1f1b_f1f1b>(nin,nout,flavours,nqed,nqcd))!=NULL) 
      s_gettermap[flc]=&Single_XS::GetProcess<XS_f1f1b_f1f1b>;
    else if ((xs=Single_XS::GetProcess<XS_f1f1b_f2f2b>(nin,nout,flavours,nqed,nqcd))!=NULL) 
      s_gettermap[flc]=&Single_XS::GetProcess<XS_f1f1b_f2f2b>;
    else if ((xs=Single_XS::GetProcess<XS_f1f2_f1f2>(nin,nout,flavours,nqed,nqcd))!=NULL) 
      s_gettermap[flc]=&Single_XS::GetProcess<XS_f1f2_f1f2>;
    else if ((xs=Single_XS::GetProcess<XS_f1f2b_f1f2b>(nin,nout,flavours,nqed,nqcd))!=NULL) 
      s_gettermap[flc]=&Single_XS::GetProcess<XS_f1f2b_f1f2b>;
    else if ((xs=Single_XS::GetProcess<XS_f1f2_f3f4>(nin,nout,flavours,nqed,nqcd))!=NULL) 
      s_gettermap[flc]=&Single_XS::GetProcess<XS_f1f2_f3f4>;
    else if ((xs=Single_XS::GetProcess<XS_f1f2b_f3f4b>(nin,nout,flavours,nqed,nqcd))!=NULL) 
      s_gettermap[flc]=&Single_XS::GetProcess<XS_f1f2b_f3f4b>;
    //SUSY-QCD processes
    else if (rpa.gen.ModelName()==std::string("MSSM")) {
      if ((xs=Single_XS::GetProcess<XS_gg_GluinoGluino>(nin,nout,flavours,model,nqed,nqcd))!=NULL);
      else if ((xs=Single_XS::GetProcess<XS_q1g_sQ1Gluino>
		(nin,nout,flavours,model,nqed,nqcd))!=NULL);
      else if ((xs=Single_XS::GetProcess<XS_gg_sQ1sQbar1>
		(nin,nout,flavours,model,nqed,nqcd))!=NULL);
      else if ((xs=Single_XS::GetProcess<XS_q1qbar1_GluinoGluino>
		(nin,nout,flavours,model,nqed,nqcd))!=NULL);
      else if ((xs=Single_XS::GetProcess<XS_q1qbar1_sQ1sQbar1>
		(nin,nout,flavours,model,nqed,nqcd))!=NULL);
      else if ((xs=Single_XS::GetProcess<XS_q1qbar1_sQ2sQbar2>
		(nin,nout,flavours,model,nqed,nqcd))!=NULL);
      else if ((xs=Single_XS::GetProcess<XS_q1q1_sQ1sQ1>
		(nin,nout,flavours,model,nqed,nqcd))!=NULL); 
      else if ((xs=Single_XS::GetProcess<XS_q1q1_sQ1LsQ1R>
		(nin,nout,flavours,model,nqed,nqcd))!=NULL); 
      else if ((xs=Single_XS::GetProcess<XS_q1q2_sQ1sQ2>
		(nin,nout,flavours,model,nqed,nqcd))!=NULL);
      else if ((xs=Single_XS::GetProcess<XS_q1qbar2_sQ1sQbar2>
		(nin,nout,flavours,model,nqed,nqcd))!=NULL);
      else if ((xs=Single_XS::GetProcess<XS_q1q2_sQ1LsQ2R>
		(nin,nout,flavours,model,nqed,nqcd))!=NULL);
    }
  if (xs!=NULL) {
    xs->m_orderEW=nqed;
    xs->m_orderQCD=nqcd;
  }
  else {
    s_gettermap[flc]=Dummy_Getter;
  }
//   if (xs)  std::cout<<" found "<<xs->Name()<<"\n";
//   else std::cout<<" found 0";
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





