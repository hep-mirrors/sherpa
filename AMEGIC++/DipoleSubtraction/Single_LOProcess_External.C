#include "AMEGIC++/DipoleSubtraction/Single_LOProcess_External.H"
#include "AMEGIC++/Amplitude/FullAmplitude_External.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "AMEGIC++/Phasespace/Phase_Space_Generator.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

#include <unistd.h>

using namespace AMEGIC;
using namespace PHASIC;
using namespace MODEL;
using namespace PDF;
using namespace BEAM;
using namespace ATOOLS;
using namespace std;

/*-------------------------------------------------------------------------------

  Constructors

  ------------------------------------------------------------------------------- */


Single_LOProcess_External::Single_LOProcess_External(const Process_Info &pi,
                                           BEAM::Beam_Spectra_Handler *const beam,
                                           PDF::ISR_Handler *const isr) :  
  Single_LOProcess(pi, beam, isr)
{
  m_emitgluon = false;
}

Single_LOProcess_External::~Single_LOProcess_External()
{
  if (p_MHVamp) delete p_MHVamp;
}


/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/



int Single_LOProcess_External::InitAmplitude(Model_Base * model,Topology* top,
					vector<Process_Base *> & links,
					vector<Process_Base *> & errs)
{
  DEBUG_FUNC("");
  m_type = 21;
  if (!model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;
  model->GetCouplings(m_cpls);
  
  m_partonlist.clear();
  for (size_t i=0;i<m_nin;i++) if (m_flavs[i].Strong()) m_partonlist.push_back(i);
  for (size_t i=m_nin;i<m_nin+m_nout;i++) if (m_flavs[i].Strong()) m_partonlist.push_back(i);

  if (m_gen_str>1) {
    ATOOLS::MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename); 
  }
  string newpath=rpa->gen.Variable("SHERPA_CPP_PATH");
  ATOOLS::MakeDir(newpath);
  if (!FileExists(newpath+"/makelibs")) {
    Copy(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/makelibs",
	     newpath+"/makelibs");
  }

  p_hel    = new Helicity(m_nin,m_nout,&m_flavs.front(),p_pl);
  
  //////////////////////////////////////////////// 

  Process_Info pi(m_pinfo);
  pi.m_fi.m_nloqcdtype=nlo_type::born;
  p_MHVamp = new FullAmplitude_External(pi,model,&m_cpls,p_hel,0,0); 
  if (p_MHVamp->Status()==0) {
    msg_Tracking()<<"Single_LOProcess_External::InitAmplitude : No process for "<<m_name<<"."<<endl;
    return 0;
  }
  m_oew=p_MHVamp->OrderEW();
  m_oqcd=p_MHVamp->OrderQCD();
  p_MHVamp->Calc()->FillCombinations(m_ccombs,m_cflavs);

  //////////////////////////////////////////////

  switch (Tests()) {
  case 1 :
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (m_allowmap && ATOOLS::IsEqual(links[j]->Result(),Result())) {
	if (CheckMapping(links[j])) {
	  msg_Tracking()<<"Single_LOProcess_External::InitAmplitude : "<<std::endl
			<<"   Found a partner for process "<<m_name<<" : "<<links[j]->Name()<<std::endl;
	  p_partner   = (Single_LOProcess_External*)links[j];
	  m_pslibname = links[j]->PSLibName();
	  break;
	}
      } 
    }
    if (p_partner==this) links.push_back(this);
    msg_Info()<<".";
    
    return 1;
  case -3: return 0;
  default :
    msg_Error()<<"ERROR in Single_Fin_Process_External::InitAmplitude : "<<std::endl
	       <<"   Failed for "<<m_name<<"."<<endl;
    errs.push_back(this);
    return 0;
  }
  return 1;
}


int Single_LOProcess_External::InitAmplitude(Model_Base * model,Topology* top,
					vector<Process_Base *> & links,
					vector<Process_Base *> & errs,
					std::vector<ATOOLS::Vec4D>* epol,std::vector<double> * pfactors)
{
  m_type = 11;
  if (!model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;
  model->GetCouplings(m_cpls);
  
  if (m_gen_str>1) {
    ATOOLS::MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename); 
  }
  string newpath=rpa->gen.Variable("SHERPA_CPP_PATH");
  ATOOLS::MakeDir(newpath);
  if (!FileExists(newpath+"/makelibs")) {
    Copy(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/makelibs",
	     newpath+"/makelibs");
  }
  int cnt=0;
  for (size_t i(0);i<m_pinfo.m_ii.m_ps.size();++i) {
    if (m_pinfo.m_ii.m_ps[i].m_tag==-1) {
      m_emit=i;
      cnt++;
    }
    if (m_pinfo.m_ii.m_ps[i].m_tag==-2) {
      m_spect=i;
      cnt+=10;
    }
  }
  for (size_t i(0);i<m_pinfo.m_fi.m_ps.size();++i) {
    if (m_pinfo.m_fi.m_ps[i].m_tag==-1) {
      m_emit=i+NIn();
      cnt++;
    }
    if (m_pinfo.m_fi.m_ps[i].m_tag==-2) {
      m_spect=i+NIn();
      cnt+=10;
    }
  }
  DEBUG_VAR(m_pinfo);
  if (cnt!=11) THROW(critical_error,"mistagged process "+m_name);
  m_emitgluon = m_flavs[m_emit].IsGluon();
  m_name+= "_S"+ToString((int)m_emit)+"_"+ToString((int)m_spect);

  if (m_flavs[m_emit].IsGluon()) {
    p_pl[m_emit]=Pol_Info();
    p_pl[m_emit].Init(2);
    p_pl[m_emit].pol_type = 'e'; 
    p_pl[m_emit].type[0] = 90;
    p_pl[m_emit].type[1] = 91;
    p_pl[m_emit].factor[0] = 1.;
    p_pl[m_emit].factor[1] = 1.;
  }

  Flavour* fl = &m_flavs.front();
  p_hel    = new Helicity(m_nin,m_nout,fl,p_pl);
  p_epol   = epol;
  //////////////////////////////////////////////// 

  Process_Info pi(m_pinfo);
  pi.m_fi.m_nloqcdtype=nlo_type::born;
  p_MHVamp = new FullAmplitude_External(pi,model,&m_cpls,p_hel,m_emit,m_spect); 
  if (p_MHVamp->Status()==0) {
    msg_Tracking()<<"Single_LOProcess_External::InitAmplitude : No process for "<<m_name<<"."<<endl;
    return 0;
  }
  m_oew=p_MHVamp->OrderEW();
  m_oqcd=p_MHVamp->OrderQCD();

  //////////////////////////////////////////////

  int tr=Tests(pfactors);
//   PRINT_INFO("Tests Result: "<<tr);
  switch (tr) {
  case 1 :
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (ATOOLS::IsEqual(links[j]->Result(),Result())) {
	if (CompareTestMoms(links[j]->GetTestMoms())) {
	  msg_Tracking()<<"Single_LOProcess_External::InitAmplitude : "<<std::endl
			<<"   Found a partner for process "<<m_name<<" : "<<links[j]->Name()<<std::endl;
	  p_partner   = (Single_LOProcess_External*)links[j];
	  m_pslibname = links[j]->PSLibName();
	  break;
	}
      } 
    }
    if (p_partner==this) links.push_back(this);
    Minimize();
   
    return 1;
  case -3: return 0;
  default :
    msg_Error()<<"ERROR in Single_LOProcess_External::InitAmplitude : "<<std::endl
	       <<"   Failed for "<<m_name<<"."<<endl;
//     errs.push_back(this);
    return 0;
  }
  return 1;
}



int Single_LOProcess_External::Tests(std::vector<double> * pfactors) {

  int number      = 1;
  int gauge_test  = 1;
  int string_test = 1;

  /* ---------------------------------------------------
     
     The reference result for momenta moms

     --------------------------------------------------- */

  string testname = string("");
  if (FoundMappingFile(testname,m_pslibname)) {
    if (testname != string("")) {
      gauge_test = string_test = 0;
    }
  }
  
  double M2 = 0.;
  double helvalue;

  if (gauge_test) {

    msg_Tracking()<<"Single_LOProcess_External::Tests for "<<m_name<<std::endl
		  <<"   Prepare gauge test and init helicity amplitudes. This may take some time."
		  <<std::endl;
    if (m_emitgluon) p_MHVamp->SetSqMatrix((*pfactors)[1],p_testmoms[GetEmit()],(*p_epol)[0]);
    p_MHVamp->Calc(p_testmoms);
    for (size_t i=0;i<p_hel->MaxHel();i++) { 
      if (p_hel->On(i) && p_hel->GetEPol(i)==90) {
	helvalue = p_MHVamp->MSquare(i)*p_hel->PolarizationFactor(i); 
	M2      +=  helvalue;
      } 
    }
     
    m_iresult  = M2;
  }
  /* ---------------------------------------------------
     
  First test : gauge test
  
  --------------------------------------------------- */
  number++;


  double M2g = 0.;
  double * M_doub = new double[p_hel->MaxHel()];
 for (size_t i=0; i<p_hel->MaxHel(); ++i) M_doub[i]=0.;
 if (m_emitgluon) p_MHVamp->SetSqMatrix((*pfactors)[1],p_testmoms[GetEmit()],(*p_epol)[0]);
 p_MHVamp->Calc(p_testmoms);
 for (size_t i=0; i<p_hel->MaxHel(); ++i) { 
     if (p_hel->On(i) && p_hel->GetEPol(i)==90) {
       M_doub[i]  = p_MHVamp->MSquare(i)*p_hel->PolarizationFactor(i); 
	 M2g       += M_doub[i];
     } 
 }


  //shorten helicities
  int switchhit = 0;
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (M_doub[i]==0. || dabs(M_doub[i]/M2g)<(ATOOLS::Accu()*1.e-2)) {
      p_hel->SwitchOff(i);
      switchhit++;
    }
  }
  msg_Tracking()<<"Single_LOProcess_External::Tests for "<<m_name<<std::endl
		<<"   Switched off or mapped "<<switchhit<<" helicities."<<std::endl;

  m_iresult  = M2g;

  if (gauge_test) {
    if (!ATOOLS::IsEqual(M2,M2g)) {
      msg_Out()<<"WARNING:  Gauge test not satisfied: "
	       <<M2<<" vs. "<<M2g<<" : "<<dabs(M2/M2g-1.)*100.<<"%"<<endl
	       <<"Gauge(1): "<<abs(M2)<<endl
	       <<"Gauge(2): "<<abs(M2g)<<endl;
    }
  }
  
  m_libname    = testname;

  /* ---------------------------------------------------
     
     Second test : string test

     --------------------------------------------------- */

  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (p_hel->On(i)) {
      for (size_t j=i+1;j<p_hel->MaxHel();j++) {
	if (p_hel->On(j)) {
	  if (ATOOLS::IsEqual(M_doub[i],M_doub[j])) {
	    p_hel->SwitchOff(j);
	    p_hel->SetPartner(i,j);
	    p_hel->IncMultiplicity(i);
	  }
	}
      }
    }
  }
  delete[] M_doub;

  return 1;
}



double Single_LOProcess_External::operator()(const ATOOLS::Vec4D_Vector &labmom,const ATOOLS::Vec4D *mom,
					std::vector<double> * pfactors,std::vector<ATOOLS::Vec4D>* epol,const int mode)
{
  DEBUG_FUNC(m_name);
  if (p_partner!=this) {
    if (m_lookup) {
      m_lastxs = p_partner->LastXS()*m_sfactor;
      if (m_lastxs!=0.) return m_lastxs;
    }
    return m_lastxs = p_partner->operator()(labmom,mom,pfactors,epol,mode)*m_sfactor;
  }
  p_int->SetMomenta(labmom);
  p_scale->CalculateScale(labmom,mode);

  double M2(0.);

  if (m_emitgluon) p_MHVamp->SetSqMatrix((*pfactors)[1],mom[GetEmit()],(*p_epol)[0]);
  p_MHVamp->Calc(mom);
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (p_hel->On(i) && p_hel->GetEPol(i)==90) {
      double mh=p_MHVamp->MSquare(i);
      mh *= p_hel->Multiplicity(i) * p_hel->PolarizationFactor(i);
      M2 += mh;
    }
  }

  m_lastxs = M2;
  return m_lastxs;
}



void Single_LOProcess_External::Calc_AllXS
(const ATOOLS::Vec4D_Vector &labmom,const ATOOLS::Vec4D *mom,
 std::vector<std::vector<double> > &dsij,const int mode) 
{
  p_int->SetMomenta(labmom);
  p_scale->CalculateScale(labmom,mode);

  dsij[0][0] =0.;
  for (size_t i=0;i<m_partonlist.size();i++) {
    for (size_t k=i+1;k<m_partonlist.size();k++) {
      dsij[k][i] =0.;
    }
  }

  p_MHVamp->Calc(mom);
  for (size_t h=0;h<p_hel->MaxHel();h++) {
    if (p_hel->On(h)) {
      double fac = p_hel->Multiplicity(h) * p_hel->PolarizationFactor(h);
      dsij[0][0] += p_MHVamp->MSquare(h,0,0)*fac;
      for (size_t i=0;i<m_partonlist.size();i++) {
	for (size_t k=i+1;k<m_partonlist.size();k++) {
	  dsij[i][k] = dsij[k][i] += p_MHVamp->MSquare(h,m_partonlist[i],m_partonlist[k])*fac;
	}
      }
    }
  }
}

int AMEGIC::Single_LOProcess_External::NumberOfDiagrams()
{ 
  return 0;
}
