#include "Helicity.H"
#include "MyComplex.H"
#include "Vector.H"
#include "MathTools.H"
#include "Message.H"
#include "Pol_Info.H"

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;

Helicity::Helicity(int Nin,int Nout,Flavour* fl,Pol_Info* pl)
{
  int N=Nin+Nout;
  int i,j,k;
  int * pnum = new int[N];


  p_pol_types = new char[N+1];
  p_angles    = new double[N];
  for(i=0;i<N;i++){
    p_pol_types[i] = pl[i].pol_type;
    p_angles[i]    = pl[i].angle;
    // mark polarised particles
    if ((pl[i].pol_type=='h' || pl[i].pol_type=='c') &&
	(pl[i].num==1)) {
      if (pl[i].type[0]==mt::p_p) p_pol_types[i] = '+';
      if (pl[i].type[0]==mt::p_m) p_pol_types[i] = '-';
    }
  }
  p_pol_types[N] = 0;  // end of char-string

  m_nsign = 1;
  for(i=0;i<N;i++) m_nsign *= pl[i].num;

  //list for calculation of physical fermion helicities
  m_fnsign       = 1;
  m_fermion_hels = 0;
  p_slist = new Sign_List[m_nsign];
  ATOOLS::msg.Debugging()<<"*****Helicity "<<p_pol_types<<":"<<endl;
  for (i=0;i<m_nsign;i++) {
      p_slist[i].s         = new int[N];
      p_slist[i].on        = 1;
      p_slist[i].multi     = 1;
      p_slist[i].polfactor = 1.;
      p_slist[i].partner   = -1;
      ATOOLS::msg.Debugging()<<i<<":";
       for (j=0;j<N;j++){
	int div = 1;
	for (k=0;k<j;k++) div *= pl[k].num;
	int l = (i/div)%pl[j].num;
	p_slist[i].s[j] = pl[j].type[l];
	if(j<Nin)p_slist[i].polfactor *= pl[j].factor[l];
	Tensor_Struc ts;
	p_slist[i].polfactor*=ts.GetTfactor(pl[j].type[l]);          //extra factor for spin2 polarisation tensor 
	ATOOLS::msg.Debugging()<<p_slist[i].s[j]<<";";
      } 
      ATOOLS::msg.Debugging()<<" "<<p_slist[i].polfactor<<endl;
  }

  if (m_fermion_hels) {
    p_fslist = new Sign_List[m_fnsign];
    for (i=0;i<m_fnsign;i++){
      p_fslist[i].s         = new int[N];
      p_fslist[i].on        = 1;
      p_fslist[i].multi     = 1;
      p_fslist[i].polfactor = 1.;
      for (j=0;j<N;j++){
	int div = 1;
	for (k=0;k<j;k++)div *= pnum[k];
	int l = (i/div)%pnum[j];
	if (fl[j].IsFermion()&&l>0) p_fslist[i].s[j] = mt::p_m;
	                       else p_fslist[i].s[j] = pl[j].type[l];
	
      }
    }
  }
  else p_fslist=p_slist;

  delete [] pnum;
}




Helicity::~Helicity() 
{
  if (p_slist!=p_fslist) {
    if (p_fslist) delete [] p_fslist;
  }
  if (p_slist) delete [] p_fslist;
  delete [] p_pol_types;
  delete [] p_angles;
}

bool Helicity::IsContrib(int i,int* pm,int length)
{
  if (!pm) return 1;

  for (int j=0;j<length;j++) if(pm[j]<99 && p_slist[i].s[j]!=pm[j]) return 0;
  return 1;
}

int Helicity::Compare(Helicity* h_cmp, int N)
{
  if (MaxHel() != h_cmp->MaxHel()) return 0;
  for (int i=0;i<MaxHel();i++) {
    for (int j=0;j<N;j++) 
      if (p_slist[i].s[j] != (*h_cmp)[i][j]) return 0;
  }
  return 1;
}

