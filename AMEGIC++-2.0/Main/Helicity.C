#include "Helicity.H"
#include "MyComplex.H"
#include "Vector.H"
#include "MathTools.H"
#include "Message.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace AMEGIC;

Helicity::Helicity(int Nin,int Nout,Flavour* fl,Pol_Info* pl)
{
  int N=Nin+Nout;
  int i,j,k;
  int * pnum = new int[N];


  p_type = new char[N+1];
  angle  = new double[N];
  for(i=0;i<N;i++){
    p_type[i] = pl[i].p_type;
    angle[i]  = pl[i].angle;
  }
  p_type[N] = 0;

  nsign = 1;
  for(i=0;i<N;i++) nsign *= pl[i].num;

  //list for calculation of physical fermion helicities
  fnsign       = 1;
  fermion_hels = 0;
  /*
  for(i=0;i<N;i++){
    pnum[i]=fl[i].isfermion() ? 2 : pl[i].num;
    fnsign*=pnum[i];
    if(p_type[i]=='h')fermion_hels=1;
  }
  */
  Slist = new signlist[nsign];
  msg.Debugging()<<"*****Helicity "<<p_type<<":"<<std::endl;
  for (i=0;i<nsign;i++){
      Slist[i].s = new int[N];
      Slist[i].on = 1;
      Slist[i].Mult = 1;
      Slist[i].polfactor = 1.;
      Slist[i].partner = -1;
      msg.Debugging()<<i<<":";
      for (j=0;j<N;j++){
	int div = 1;
	for (k=0;k<j;k++) div *= pl[k].num;
	int l = (i/div)%pl[j].num;
	//msg.Debugging()<<div<<","<<l<<";";
	Slist[i].s[j] = pl[j].type[l];
	if(j<Nin)Slist[i].polfactor *= pl[j].factor[l];
	msg.Debugging()<<Slist[i].s[j]<<";";
      }
      msg.Debugging()<<" "<<Slist[i].polfactor<<std::endl;
  }

  if(fermion_hels){
    FSlist = new signlist[fnsign];
    msg.Debugging()<<"*****FHelicity "<<p_type<<":"<<std::endl;
    for (i=0;i<fnsign;i++){
      FSlist[i].s         = new int[N];
      FSlist[i].on        = 1;
      FSlist[i].Mult      = 1;
      FSlist[i].polfactor = 1.;
      msg.Debugging()<<i<<":";
      for (j=0;j<N;j++){
	int div = 1;
	for (k=0;k<j;k++)div *= pnum[k];
	int l = (i/div)%pnum[j];
	//msg.Debugging()<<div<<","<<l<<";";
	if (fl[j].isfermion()&&l>0) FSlist[i].s[j] = mt::p_m;
	                       else FSlist[i].s[j] = pl[j].type[l];
	
	msg.Debugging()<<FSlist[i].s[j]<<";";
      }
      msg.Debugging()<<std::endl;
    }
  }
  else FSlist=Slist;

  delete [] pnum;
}










