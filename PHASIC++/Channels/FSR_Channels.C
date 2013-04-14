#include "PHASIC++/Channels/FSR_Channels.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Channels/FSR_Channel.H"
#include "PHASIC++/Channels/Rambo.H"
#include "PHASIC++/Channels/RamboKK.H"
#include "PHASIC++/Channels/Sarge.H"
#include "PHASIC++/Channels/VHAAG.H"
#include "PHASIC++/Channels/VHAAG_ND.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Channels/VHAAG_res.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Permutation.H"

using namespace PHASIC;
using namespace ATOOLS;

bool FSR_Channels::Initialize()
{
  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(rpa->GetPath());
  dr.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
  m_inttype=dr.GetValue<int>("INTEGRATOR",6);
  nin=p_psh->Process()->NIn();
  nout=p_psh->Process()->NOut();
  int m_nin(nin), m_nout(nout);
  if (nin==1) {
    if (nout==2) m_inttype = 0;
    if (m_inttype<4)  m_inttype = 0;
    else m_inttype = 4;
  }
  if (nin==2) { 
    if (nout==2&&m_inttype==3) m_inttype=7;
    if (p_psh->Process()->Process()->Info().m_fi.m_nloqcdtype&nlo_type::real) {
      if (m_inttype==2) m_inttype=3;
      if (m_inttype>=4&&m_inttype<7) m_inttype=7;
    }
    if (m_inttype<4||m_inttype>20) DropAllChannels();
  }
  if (p_psh->Process()->NOut()==1) {
    if (!p_psh->Process()->Process()->InitIntegrator(p_psh))
      THROW(critical_error,"InitIntegrator failed");
    return false;
  }
  bool sintegrator(false);
  switch (m_inttype) {
  case 0:
    {
      bool kk_fs=false;
      for (int i=0;i<m_nout;i++){
	if (p_psh->Flavs()[i+m_nin].IsKK()) kk_fs=true;
      }
      if (kk_fs) {
	Add(new RamboKK(m_nin,m_nout,&p_psh->Flavs().front()));
	break;
      }
    }
      
    if (m_nin==1 && m_nout==2) Add(new Decay2Channel(m_nin,m_nout,&p_psh->Flavs().front()));
    else Add(new Rambo(m_nin,m_nout,&p_psh->Flavs().front(),
                       p_psh->Process()->Process()->Generator()));
    break;
  case 1: 
    {
	VHAAG_res *firsthaag=NULL,*hlp=NULL;
	Permutation pp(m_nin+m_nout-3);
	for (int j=0;j<pp.MaxNumber();j++) {
	  Add(hlp=new VHAAG_res(m_nin,m_nout,2*j,firsthaag));
	  if (!firsthaag) firsthaag=hlp;
	  Add(hlp=new VHAAG_res(m_nin,m_nout,2*j+1,firsthaag));
	  if (!firsthaag) firsthaag=hlp;
 	}
    }
    break;
  case 2:
    {
      if (m_nin==2 && m_nout==2) {
	Add(new S1Channel(m_nin,m_nout,(Flavour*)&p_psh->Flavs().front()));
	Add(new T1Channel(m_nin,m_nout,(Flavour*)&p_psh->Flavs().front()));
	Add(new U1Channel(m_nin,m_nout,(Flavour*)&p_psh->Flavs().front()));
      }
      else {
      VHAAG *firsthaag=NULL,*hlp=NULL;
      Permutation pp(m_nin+m_nout-1);
      for (int j=0;j<pp.MaxNumber();j++) {
	int* pm = pp.Get(j);
	if (pm[1]==0||pm[m_nin+m_nout-3]==0) 
	  Add(hlp=new VHAAG(m_nin,m_nout,j,firsthaag));
	if (!firsthaag) firsthaag=hlp;
      }
      }
    }
    break;
  case 3: 
    {
      VHAAG_ND *firsthaag=NULL,*hlp=NULL;
      Permutation pp(m_nin+m_nout-1);
      for (int j=0;j<pp.MaxNumber();j++) {
	int* pm = pp.Get(j);
	if (pm[1]!=0&&pm[m_nin+m_nout-3]!=0) 
	  Add(hlp=new VHAAG_ND(m_nin,m_nout,j,firsthaag));
	if (!firsthaag) firsthaag=hlp;
      }
    }
    break;
  case 4:case 5:case 6:case 7: 
    DropRedundantChannels();
    sintegrator=p_psh->Process()->Process()->FillIntegrator(p_psh);
    break;
  default:
    msg_Error()<<"Wrong phasespace integration switch ! Using RAMBO as default."<<std::endl;
    Add(new Rambo(m_nin,m_nout,&p_psh->Flavs().front(),
                  p_psh->Process()->Process()->Generator()));
  }  
  if (!p_psh->Process()->Process()->InitIntegrator(p_psh))
    THROW(critical_error,"InitIntegrator failed");
  return sintegrator;
}

void FSR_Channels::DropRedundantChannels()
{
  Reset();
  int number = Number();
  if (number<2) return;
  int *marker = new int[number];  
  for (short int i=0;i<number;i++) marker[i] = 0; 
  /*Vec4D** perm_vec = new Vec4D*[number]; 
  for (short int i=0;i<number;i++) perm_vec[i] = new Vec4D[m_nin+m_nout+1];
  // Create Momenta
  int rannum   = 1 + 2 + 3*(m_nout-2);
  double * rans = new double[rannum];
  for (short int i=0;i<rannum;i++) rans[i] = ran->Get();  
  // Init markers for deletion and results to compare.
  double * res    = new double[number];
  for (short int i=0;i<number;i++) { marker[i] = 0;res[i] = 0.; }
  for (short int i=0;i<number;i++) {
    perm_vec[i][0] = Vec4D(rpa->gen.Ecms()/2.,0.,0.,rpa->gen.Ecms()/2.);
    perm_vec[i][1] = Vec4D(rpa->gen.Ecms()/2.,0.,0.,-rpa->gen.Ecms()/2.); 
    p_fsrchannels->GeneratePoint(i,perm_vec[i],p_process->Cuts(),rans);
    p_fsrchannels->GenerateWeight(i,perm_vec[i],p_process->Cuts());
    res[i] = p_fsrchannels->Weight();
    if (res[i]==0.) marker[i] = 1;
  }
  delete[] rans;*/
  // kick identicals & permutations
  for (short int i=0;i<number;i++) {
    if (marker[i]==0) {
      for (short int j=i+1;j<number;j++) {
	if (marker[j]==0) {
	  //if ( (Compare(perm_vec[i],perm_vec[j])) && 
	  //     (ATOOLS::IsEqual(res[i],res[j])) ) {
	  if (CompareCh(ChID(i),ChID(j))) {
	    marker[j] = 1; 
	  }
	}
      }
    }
  }
  // kick non-resonants
  /*
    int max_reson    = 0;
    Flavour * fl_res = 0;
    int * reson      = new int[number];
    for (short int i=0;i<number;i++) {
    if (marker[i]==0) {
    reson[i]     = fsrchannels->CountResonances(i,fl_res);
    if (reson[i]!=0) {
    //shorten
    int hit    = 0;
    for (short int j=0;j<reson[i];j++) {
    if (sqr(fl_res[j].Mass())>ycut*sqr(rpa->gen.Ecms()) &&
    sqr(fl_res[j].Mass())<sqr(rpa->gen.Ecms())) 
    hit++;
    }
    reson[i] = hit;
    if (reson[i]>max_reson) max_reson = reson[i];
    }
    else reson[i] = -1;
    }
    else reson[i] = -1;
    }
    //Drop them
    for (short int i=0;i<number;i++) {
    if (reson[i]<max_reson && reson[i]!=-1) marker[i] = 1;
    }
    delete [] reson;
  */
  int count = 0;
  for (short int i=0;i<number;i++) {
    if (marker[i]) {
      DropChannel(i-count);
      count++;
    }
  }
  //delete [] res;
  delete [] marker; 
  //for (short int i=0;i<number;i++) delete [] perm_vec[i];
  //delete [] perm_vec; 
}
  
bool FSR_Channels::CompareCh(std::string C1,std::string C2)
{
  int l=Min(C1.length(),C1.length());
  for (int i=0;i<l;i++) {
    if (C1[i]!=C2[i]) return 0;
    if (C1[i]=='Z') return 1;
  }
  return 1;
}


bool FSR_Channels::Compare(const Vec4D *p1,const Vec4D *p2)
{
  int m_nin(p_psh->Process()->NIn());
  int m_nout(p_psh->Process()->NOut());
  if (m_nout==2) {
    for (short int i=0;i<m_nout;i++) { 
      if (p1[m_nin+i] != p2[m_nin+i]) return 0;
    }
    return 1;
  }
  else {
    //Identicals
    for (short int i=0;i<m_nout;i++) {
      if (p1[i+m_nin] != p2[i+m_nin]) return 0;
    }
    return 1;
    //Permutations - not reached right now.
    int * perm = new int[m_nout];
    for (short int i=0;i<m_nout;i++) perm[i] = 0; 
    
    int over = 0;
    int hit,sw1;
    for(;;) {
      sw1 = 1;
      for(short int i=0;i<m_nout;i++) {
	for (short int j=i+1;j<m_nout;j++) 
	  if (perm[i]==perm[j]) {sw1 = 0; break;}
      }    
      if (sw1) {
	hit = 1;
	for (short int i=0;i<m_nout;i++) {
	  if (p1[i+m_nin] != p2[perm[i]+m_nin]) {
	    hit = 0;
	    break;
	  }
	}
	if (hit) return 1;
      }
      for (short int j=m_nout-1;j>=0;j--) {
	if ((perm[j]+1)<m_nout) {
	  perm[j]++;            
	  break;
	}
	else {
	  perm[j] = 0;
	  if (j==0) over = 1;
	}
      }
      if (over) break;
    }
    delete[] perm;
    return 0;
  }
}

