#include "Helicity.H"
#include "MyComplex.H"
#include "Vector.H"
#include "MathTools.H"
#include "Message.H"
#include "Pol_Info.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace AMEGIC;
using namespace std;

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
  Slist = new signlist[nsign];
  for (i=0;i<nsign;i++){
      Slist[i].s         = new int[N];
      Slist[i].on        = 1;
      Slist[i].Mult      = 1;
      Slist[i].polfactor = 1.;
      Slist[i].partner   = -1;
      for (j=0;j<N;j++){
	int div = 1;
	for (k=0;k<j;k++) div *= pl[k].num;
	int l = (i/div)%pl[j].num;
	Slist[i].s[j] = pl[j].type[l];
	if(j<Nin)Slist[i].polfactor *= pl[j].factor[l];
	Tensor_Struc ts;
	Slist[i].polfactor*=ts.GetTfactor(pl[j].type[l]);          //extra factor for spin2 polarisation tensor 
      } 
  }

  if(fermion_hels){
    FSlist = new signlist[fnsign];
    for (i=0;i<fnsign;i++){
      FSlist[i].s         = new int[N];
      FSlist[i].on        = 1;
      FSlist[i].Mult      = 1;
      FSlist[i].polfactor = 1.;
      for (j=0;j<N;j++){
	int div = 1;
	for (k=0;k<j;k++)div *= pnum[k];
	int l = (i/div)%pnum[j];
	if (fl[j].IsFermion()&&l>0) FSlist[i].s[j] = mt::p_m;
	                       else FSlist[i].s[j] = pl[j].type[l];
	
      }
    }
  }
  else FSlist=Slist;

  delete [] pnum;
}










