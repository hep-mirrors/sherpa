#include "PHASIC++/Channels/FSR_Channels.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Channels/Channel_Generator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PHASIC;
using namespace ATOOLS;

bool FSR_Channels::Initialize()
{
  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  dr.AddWordSeparator(",");
  dr.SetInputPath(rpa->GetPath());
  dr.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
  std::vector<std::string> inttypes;
  dr.VectorFromFile(inttypes,"INTEGRATOR");
  if (inttypes.empty()) inttypes.push_back("Default");
  if (p_psh->Process()->Process()->Info().m_integrator!="") {
    dr.SetString(p_psh->Process()->Process()->Info().m_integrator);
    dr.VectorFromString(inttypes,"");
  }
  nin=p_psh->Process()->NIn();
  nout=p_psh->Process()->NOut();
  int m_nin(nin), m_nout(nout);
  int sintegrator=0;
  for (size_t i(0);i<inttypes.size();++i) {
    Channel_Generator *cg=
      Channel_Generator::Getter_Function::GetObject
      (inttypes[i],Channel_Generator_Key
       (p_psh->Process()->Process(),this));
    if (cg==NULL) {
      s_loader->AddPath(rpa->gen.Variable("SHERPA_LIB_PATH"));
      if (s_loader->LoadLibrary("Proc_"+inttypes[i])) {
	cg=Channel_Generator::Getter_Function::GetObject
	  (inttypes[i],Channel_Generator_Key
	   (p_psh->Process()->Process(),this));
      }
      if (cg==NULL) {
	THROW(fatal_error,"Channel generator '"
	      +inttypes[i]+"' not found.");
      }
    }
    sintegrator|=cg->GenerateChannels();
    delete cg;
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

