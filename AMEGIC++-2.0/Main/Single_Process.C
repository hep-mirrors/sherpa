#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "Single_Process.H"
#include "Run_Parameter.H"
#include "ISR_Base.H"

#include "Running_AlphaS.H"
#include "Combined_Selector.H"

#include "prof.hh"
#include "Test_Selector.H"

using namespace AMEGIC;
using namespace PHASIC;
using namespace ISR;
using namespace BEAM;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

int fak(int N)
{
  if (N == 0) return 1;
  if (N < 0) return 0;  
  int res =1;
  for (int i=1;i<=N;i++) res *= i;
  return res;
}

/*-------------------------------------------------------------------------------

  Constructors

  ------------------------------------------------------------------------------- */

Single_Process::Single_Process(int _nin,int _nout,Flavour * _fl,
			       ISR_Handler * _isr,Beam_Spectra_Handler * _beam,Selector_Data * _seldata,
			       int _gen_str,int _orderQCD, int _orderEW,
			       int _kfactorscheme, int _scalescheme,double _scalefactor,
			       Pol_Info * _pl,int _nex,Flavour * _ex_fl) :
  Process_Base(_nin,_nout,_fl,_isr,_beam,_gen_str,_orderQCD,_orderEW,
	       _scalescheme,_kfactorscheme,_scalefactor,_pl,_nex,_ex_fl),
  p_partner(this)
{
  GenerateNames(m_nin,p_flin,p_plin,m_nout,p_flout,p_plout,m_name,m_ptypename,m_libname);

  PolarizationNorm();
  InitCuts();
  if (_seldata) p_sel = new Combined_Selector(m_nin,m_nout,p_fl,_seldata);
  else {
    msg.Error()<<"Potential Error in Single_Process "<<m_name<<endl
	       <<"   No selection cuts specified. Init No_Selector !"<<endl;
    p_sel = new No_Selector();
  }

  p_ps   = new Phase_Space_Handler(this,p_isr,p_beam);
  
  // making directory
  int  mode_dir = 448;
  mkdir((string("Process/")+m_ptypename).c_str(),mode_dir); 
  
  msg.Tracking()<<"Initialized Single_Process : "<<m_name<<", "<<m_nvec<<", 1/norm = "<<1./m_Norm<<endl
		<<"String handling : "<<m_gen_str<<endl;;
}


Single_Process::~Single_Process()
{
  if (p_hel)      {delete p_hel; p_hel=0;}
  if (p_BS)       {delete p_BS;   p_BS=0;}
  if (p_shand)    {delete p_shand;p_shand=0;}
  if (p_ampl)     {delete p_ampl; p_ampl=0;}
}

/*------------------------------------------------------------------------------
  
  Generic stuff for initialization of Single_Processes
      
  ------------------------------------------------------------------------------*/

void Single_Process::PolarizationNorm() {
  int                   is_massless_pol  = m_pol.Massless_Vectors(m_nin,p_flin);
  if (!is_massless_pol) is_massless_pol  = m_pol.Massless_Vectors(m_nout,p_flout);
  int                   nmassive_pols    = m_pol.Massive_Vectors(m_nin,p_flin);
  nmassive_pols                         += m_pol.Massive_Vectors(m_nout,p_flout);

  // arrange Flavours
  m_nvec = m_nin+m_nout+is_massless_pol+nmassive_pols;
  p_fl   = new Flavour[m_nvec];
  p_pl   = new Pol_Info[m_nvec];
  p_b    = new int[m_nvec];
  for (short int i=0;i<m_nin;i++)             { p_fl[i] = p_flin[i]       ; p_pl[i] = p_plin[i]       ; p_b[i] = -1; }
  for (short int i=m_nin;i<m_nin+m_nout;i++)  { p_fl[i] = p_flout[i-m_nin]; p_pl[i] = p_plout[i-m_nin]; p_b[i] = 1; } 
  for (short int i=m_nin+m_nout;i<m_nvec;i++) { p_fl[i] = Flavour(kf::pol); p_b[i]  = 1; }

  m_Norm = SymmetryFactors() * m_pol.Spin_Average(m_nin,p_flin);
#ifndef Explicit_Pols
  m_pol.Attach(m_nin+m_nout,p_fl);
  m_Norm *= m_pol.Massive_Norm();

  msg.Debugging()<<"In Single_Process::PolarizationNorm for "<<m_name<<" : "<<is_massless_pol<<" "
		 <<m_nin<<" -> "<<m_nout<<" ("<<m_nvec<<")"<<endl;
#endif
}

double Single_Process::SymmetryFactors()
{
  double sym = 1.;
  Fl_Iter fli;
  for (Flavour hflav=fli.first();hflav!=Flavour(kf::none);hflav = fli.next()) {
    // no need to compare Hadrons: if (hflav.IsHadron()) break;
    if (hflav==Flavour(kf::pi)) break; 
    int cp  = 0;
    int cap = 0;
    for (int j=0;j<m_nout;j++) {
      if (p_flout[j]==hflav)                                      cp++;
      else {
	if ((p_flout[j]==hflav.Bar()) && (hflav != hflav.Bar()))  cap++;
      }
    }
    if (cp>1)  sym *= double(fak(cp));
    if (cap>1) sym *= double(fak(cap));
  } 
  return 1./sym;
}

void Single_Process::FixISRThreshold()
{
  double m_mass_in  = 0.;
  double m_mass_out = 0.;
  
  for (int i = 0;i<m_nin;i++)  m_mass_in  += p_flin[i].Mass(); 
  for (int i = 0;i<m_nout;i++) m_mass_out += p_flout[i].Mass(); 
  
  double isrth = AMATOOLS::Max(m_mass_in,m_mass_out);
  
  SetISRThreshold(isrth);

}

/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/



int Single_Process::InitAmplitude(Interaction_Model_Base * model,Topology* top,Vec4D *& _testmoms,
				  vector<double> & results,vector<Single_Process *> & links)
{
  if (_testmoms==0) {
    msg.Debugging()<<"Init moms : "<<_testmoms<<" : "<<m_nin+m_nout<<endl;
    _testmoms = new Vec4D[m_nvec];
    p_ps->TestPoint(_testmoms);
    for (int i=0;i<m_nin+m_nout;i++)
      msg.Debugging()<<i<<" th mom : "<<_testmoms[i]<<" ("<<_testmoms[i].Abs2()<<")"<<endl;
  }
  if (p_moms) { delete [] p_moms; }
  p_moms = new Vec4D[m_nvec]; 
  for (int i=0;i<m_nin+m_nout;i++) p_moms[i] = _testmoms[i];

  double result = 0.;
  p_hel    = new Helicity(m_nin,m_nout,p_fl,p_pl);
  p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nvec,p_fl,p_b);  
  p_shand  = new String_Handler(m_gen_str,p_BS);

  p_ampl   = new Amplitude_Handler(m_nin+m_nout,p_fl,p_b,&m_pol,model,top,m_orderQCD,m_orderEW,
				   p_BS,p_shand,m_ptypename+string("/")+m_name);
  if (p_ampl->GetGraphNumber()==0) {
    msg.Error()<<"Single_Process::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return -1;
  }
  m_pol.Add_Extern_Polarisations(p_BS,p_fl,p_hel);
  p_BS->Initialize();
  p_shand->Initialize(p_ampl->GetRealGraphNumber(),p_hel->Max_Hel());

  switch (Tests(result)) {
  case 2 : 
    for (int j=0;j<results.size();j++) {
      if (AMATOOLS::IsZero((results[j]-result)/(results[j]+result))) {
	msg.Tracking()<<"Test : 2.  Can map "<<m_name<<" on "<<links[j]->Name()<<endl;
	p_partner = links[j];
      }
    }
    if (p_partner==this) {
      results.push_back(result);
      links.push_back(this);
    }
    return 1;
  case 1 :
    for (int j=0;j<results.size();j++) {
      if (AMATOOLS::IsZero((results[j]-result)/(results[j]+result))) {
	msg.Tracking()<<"Test : 1.  Can map "<<m_name<<" on "<<links[j]->Name()<<endl;
	p_partner = links[j];
      }
    }
    if (p_partner==this) {
      results.push_back(result);
      links.push_back(this);
    }
    return InitLibrary(result);
  default :
    msg.Error()<<"Error in Single_Process::InitAmplitude : Failed for "<<m_name<<"."<<endl;
    return -2;
  }
}

void Single_Process::InitDecay(Topology* top)
{
  /*  
  top          = new Topology(m_nin+m_nout);
 
  p_hel    = new Helicity(m_nin,m_nout,fl,pl);
  p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nvec,fl,b);  
  p_shand  = new String_Handler(m_gen_str,p_BS);

  p_ampl   = new Amplitude_Handler(m_nin+m_nout,fl,b,&m_pol,top,p_BS,p_shand,m_ptypename+string("/")+m_name);
  if (p_ampl->GetGraphNumber()==0) {
    msg.Error()<<"Single_Process::InitDecays : No diagrams for "<<name<<"."<<endl;
    return;
  }
  m_pol.Add_Extern_Polarisations(p_BS,fl,p_hel);
  p_BS->Initialize();
  p_shand->Initialize(p_ampl->GetGraphNumber(),p_hel->Max_Hel());
  */
}

int Single_Process::Tests(double & result) {
  int number      = 1;
  int gauge_test  = 1;
  int string_test = 1;

  /* ---------------------------------------------------
     
     The reference result for momenta moms

     --------------------------------------------------- */

  string testname = string("");
  if (FoundMappingFile(testname)) {
    if (testname != string("")) {
      gauge_test = string_test = 0;
    }
  }
  
  p_ampl->SetStringOff();

  double M2 = 0.;
  if (gauge_test) {
    m_pol.Set_Gauge_Vectors(m_nin+m_nout,p_moms,Vec4D(sqrt(3.),1.,1.,-1.));
    p_BS->Setk0(0);
    p_BS->CalcEtaMu(p_moms);  
    p_BS->InitGaugeTest(.9);

    msg.Debugging()<<number<<" :";AORGTOOLS::msg.Debugging().flush();
    for (short int i=0;i<p_hel->Max_Hel();i++) { 
      if (p_hel->On(i)) {
	M2 +=  p_ampl->Differential(i,(*p_hel)[i])*p_hel->PolarizationFactor(i); 
	msg.Debugging()<<"*";msg.Debugging().flush();
      } 
      else {
	msg.Debugging()<<"0";msg.Debugging().flush();
      }
    }
    msg.Debugging()<<endl;
    M2     *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
    result  = M2;
  }

    p_ampl->ClearCalcList();
    // To prepare for the string test.
    p_ampl->SetStringOn();
    (p_shand->Get_Generator())->Reset();
    p_ampl->FillCoupling(p_shand);


  /* ---------------------------------------------------
     
     First test : gauge test

     --------------------------------------------------- */
#ifndef Explicit_Pols 
  Vec4D gauge(sqrt(3.),1.,1.,1.);
  if (m_nout==4) gauge = p_moms[4];
  if (m_nout==5) gauge = p_moms[4];
  //?????????
  if (m_nout==6) gauge = p_moms[4];  
  m_pol.Reset_Gauge_Vectors(m_nin+m_nout,p_moms,gauge);
#else
  p_BS->Setk0(1);
#endif

  p_BS->CalcEtaMu(p_moms);
  number++;
  msg.Debugging()<<number<<" :";

  if (!gauge_test) p_ampl->SetStringOff();  //second test without string production 

  double M2g = 0.;
  double * M_doub = new double[p_hel->Max_Hel()];

  for (short int i=0;i<p_hel->Max_Hel();i++) { 
    if (p_hel->On(i)) {
      M_doub[i]  = p_ampl->Differential(i,(*p_hel)[i])*p_hel->PolarizationFactor(i);  
      M2g       += M_doub[i];
      msg.Debugging()<<"*";msg.Debugging().flush();
    }
  }
  msg.Debugging()<<endl;

  //shorten helicities
  for (short int i=0;i<p_hel->Max_Hel();i++) {
    if (M_doub[i]/M2g<1.e-30) {
      p_hel->switch_off(i);
      msg.Debugging()<<"Switch off zero helicity "<<i<<" : "
		     <</*p_ampl->Differential(i,(*p_hel)[i])<<"/"<<*/M_doub[i]/M2g<<endl;
    }
    /*
      //may cause problems due to numerics
    if (p_hel->On(i)) {
      for (short int j=i+1;j<p_hel->Max_Hel();j++) {
	if (p_hel->On(j)) {
	  if (AMATOOLS::IsEqual(M_doub[i]/M2g,M_doub[j]/M2g)) {
	    msg.Debugging()<<"Switch off helicity "<<j<<" -> "<<i<<" : "
		   	   <<p_ampl->Differential(i,(*p_hel)[i])<<"/"
			  <<M_doub[i]/M2g<<" -> "<<M_doub[j]/M2g<<endl;
			  
	    p_hel->switch_off(j);
	    p_hel->SetPartner(i,j);
	    p_hel->Inc_Mult(i);
	  }
	}
      }
    }*/
  }
  M2g    *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
  result  = M2g;

  delete[] M_doub;
  p_ampl->ClearCalcList();  
  p_BS->StartPrecalc();

  if (gauge_test) {
    msg.Tracking()<<"Gauge(1): "<<abs(M2)<<endl
		  <<"Gauge(2): "<<abs(M2g)<<endl
		  <<"Gauge test: "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
    if (!AMATOOLS::IsZero(abs(M2/M2g-1.))) {
      msg.Tracking()<<"Gauge test not satisfied: "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
    }
  }
  else {
    number++;
    if (p_shand->SearchValues(m_gen_str,testname,p_BS)) {
      m_pol.Set_Gauge_Vectors(m_nin+m_nout,p_moms,Vec4D(sqrt(3.),1.,1.,-1.));
      p_BS->CalcEtaMu(p_moms);  
      p_shand->Initialize(p_ampl->GetRealGraphNumber(),p_hel->Max_Hel());
      (p_shand->Get_Generator())->Reset();
      p_ampl->FillCoupling(p_shand);
      p_shand->Complete(p_hel);
      M2 = operator()(p_moms);
      gauge_test = string_test = 0;
    }
    msg.Debugging()<<"Mapping file(1) : "<<abs(M2)<<endl
		   <<"Original    (2) : "<<abs(M2g)<<endl
		   <<"Cross check     : "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
    if (!AMATOOLS::IsZero(abs(M2/M2g-1.))) {
      msg.Tracking()<<"Cross check not satisfied: "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
      return 0;
    }
    m_libname    = testname;
    return 2;
  }

  /* ---------------------------------------------------
     
     Second test : string test

     --------------------------------------------------- */

  {
    PROFILE_LOCAL("Shand.Complete()");
    p_shand->Complete(p_hel);
  }
  if (string_test) {
    //String-Test
    if (p_shand->Is_String()) {
      double  M2S = 0.;
      p_shand->Calculate();
      msg.Debugging()<<"3:";msg.Debugging().flush();
      for (short int i=0;i<p_hel->Max_Hel();i++) {
	if (p_hel->On(i)) {
	  msg.Debugging()<<"*";msg.Debugging().flush();
	  M2S      += p_ampl->Differential(i)*p_hel->PolarizationFactor(i)*p_hel->Multiplicity(i);
	}
      }
      msg.Debugging()<<endl;
      M2S *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
      msg.Tracking()<<"String test: "<<abs(M2g/M2S-1.)*100.<<"%"<<endl;
      if (!AMATOOLS::IsZero(abs(M2g/M2S-1.))) {
	msg.Tracking()<<"String test not satisfied!!"<<endl;
	return 0;
      }
      return 1;
    }
    msg.Tracking()<<"no strings created."<<endl;
    return 0;
  }
  msg.Tracking()<<"string_test switched off !"<<endl;
  return 0;
}

int Single_Process::InitLibrary(double result) {
  if (p_shand->IsLibrary()) {
    msg.Tracking()<<"No library needs to be initialized. Already done."<<endl;
    return 1;
  }

  char help[20];
  sprintf(help,"%i",p_ampl->GetGraphNumber());
  m_libname += string("_");
  m_libname += string(help);
  sprintf(help,"%i",p_shand->NumberOfCouplings());
  m_libname += string("_");
  m_libname += string(help);
  sprintf(help,"%i",p_shand->NumberOfZfuncs());
  m_libname += string("_");
  m_libname += string(help);
  msg.Debugging()<<"In Init library for "<<m_name<<endl;


  String_Handler * shand1;
  shand1      = new String_Handler(p_shand->Get_Generator());
  
  int number  = 0;
  string proc = string("Process/")+m_ptypename+string("/V");
  string testname;
  double M2s;

  for (;;) {
    ++number;
    sprintf(help,"%i",number);
    testname  = m_libname+string("_")+string(help);
    if (shand1->SearchValues(m_gen_str,testname,p_BS)) {
      shand1->Initialize(p_ampl->GetRealGraphNumber(),p_hel->Max_Hel());
      (shand1->Get_Generator())->Reset();
      p_ampl->FillCoupling(shand1);
      shand1->Calculate();
      
      M2s = 0.;
      AORGTOOLS::msg.Debugging()<<"Check "<<number<<" :";AORGTOOLS::msg.Debugging().flush();
      for (short int i=0;i<p_hel->Max_Hel();i++) {
	M2s     += p_ampl->Differential(shand1,i) * p_hel->PolarizationFactor(i) *
	  p_hel->Multiplicity(i);
	msg.Debugging()<<"*";AORGTOOLS::msg.Debugging().flush();
      }
      msg.Debugging()<<endl;
      M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
      msg.Debugging()<<"Cross check: "<<abs(M2s/result-1.)*100.<<"%"<<"  : "
		     <<M2s<<"/"<<result<<endl;
      if (AMATOOLS::IsZero(abs((M2s-result)/(M2s+result)))) {
	msg.Tracking()<<"Found a suitable string."<<endl;
	m_libname = testname;
	if (p_shand->SearchValues(m_gen_str,testname,p_BS)) {
	  p_shand->Initialize(p_ampl->GetRealGraphNumber(),p_hel->Max_Hel());
	  (p_shand->Get_Generator())->Reset();
	  p_ampl->FillCoupling(p_shand);
	  
	  p_shand->Calculate();
	  M2s = 0.;
	  AORGTOOLS::msg.Debugging()<<number<<" :";AORGTOOLS::msg.Debugging().flush();
	  for (short int i=0;i<p_hel->Max_Hel();i++) {
	    M2s     += p_ampl->Differential(p_shand,i) * p_hel->PolarizationFactor(i) * 
	      p_hel->Multiplicity(i);
	    msg.Debugging()<<"*";AORGTOOLS::msg.Debugging().flush();
	  }
	  msg.Debugging()<<endl;
	  M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
	  msg.Tracking()<<"Cross check: "<<abs(M2s/result-1.)*100.<<"%"<<endl;
	  
	  CreateMappingFile();
	  if (shand1) { delete shand1; shand1 = 0; }
	  return 1;
	}
	else {
	  msg.Error()<<"Error in Single_Process::InitLibrary for "<<testname<<endl;
	  abort();
	}
      }
    }
    else break;
  }
  
  for (;;) {
    sprintf(help,"%i",number);
    testname    = m_libname+string("_")+string(help);
    if (!(IsFile(string("Process/")+m_ptypename+string("/")+testname+string("/V.H")))) break;
    ++number;
  }
  m_libname = testname;
  if (p_partner==this) {
    msg.Debugging()<<"Write Library for "<<m_name<<" = "<<m_libname<<endl;
    int  mode_dir = 448;
    msg.Debugging()<<" m_ptypename = "<<m_ptypename<<endl<<" m_libname = "<<m_libname<<endl;
    mkdir((string("Process/")+m_ptypename+string("/")+m_libname).c_str(),mode_dir); 
    p_shand->Output(p_hel,m_ptypename+string("/")+m_libname);
    CreateMappingFile();
    if (shand1) { delete shand1; shand1 = 0; }
    return 0;
  }
  else {
    if (p_partner->p_shand->IsLibrary()) {
      if (p_partner->p_shand->SearchValues(m_gen_str,p_partner->m_libname,p_BS)) {
	p_partner->p_shand->Initialize(p_ampl->GetRealGraphNumber(),p_hel->Max_Hel());
	(p_partner->p_shand->Get_Generator())->Reset();
	p_ampl->FillCoupling(p_partner->p_shand);
	p_partner->p_shand->Calculate();
	
	M2s = 0.;
	AORGTOOLS::msg.Debugging()<<"Check "<<number<<" :";AORGTOOLS::msg.Debugging().flush();
	for (short int i=0;i<p_hel->Max_Hel();i++) {
	  M2s     += p_ampl->Differential(p_partner->p_shand,i) * p_hel->PolarizationFactor(i) *
	    p_hel->Multiplicity(i);
	  msg.Debugging()<<"*";AORGTOOLS::msg.Debugging().flush();
	}
	msg.Debugging()<<endl;
	M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
	msg.Debugging()<<"Cross check: "<<abs(M2s/result-1.)*100.<<"%"<<"  : "
		       <<M2s<<"/"<<result<<endl;
	if (AMATOOLS::IsZero(abs((M2s-result)/(M2s+result)))) {
	  m_libname = p_partner->m_libname;
	  CreateMappingFile();
	}
	else {
	  m_libname = testname;
	  msg.Debugging()<<"Write Library for "<<m_name<<" = "<<m_libname<<endl;
	  int  mode_dir = 448;
	  msg.Debugging()<<" m_ptypename = "<<m_ptypename<<endl<<" m_libname = "<<m_libname<<endl;
	  mkdir((string("Process/")+m_ptypename+string("/")+m_libname).c_str(),mode_dir); 
	  p_shand->Output(p_hel,m_ptypename+string("/")+m_libname);
	  CreateMappingFile();
	  if (shand1) { delete shand1; shand1 = 0; }
	  return 0;
	}
      }
      else {
	msg.Error()<<"Partner with no library."<<endl;
	abort();
      }
    }
    return 0;
  }
}

void Single_Process::CreateMappingFile() {
  char outname[100];
  sprintf(outname,"%s.map",(string("Process/")+m_ptypename+string("/")+m_name).c_str());
  if (IsFile(outname)) {
    msg.Debugging()<<"Mapping file already exists. Check identity."<<endl;
    ifstream from;
    from.open(outname);
    string tempname;
    from>>tempname;
    if (tempname != m_libname) {
      msg.Error()<<"Files do not coincide. Maybe changed input data ?"<<endl;
      abort();
    }
    else return;
  }
  msg.Debugging()<<" Write Mapping Information to : "<<outname<<endl;

  ofstream to;
  to.open(outname,ios::out);
  to<<m_libname<<endl;
  to.close();
  msg.Debugging()<<"File "<<outname<<" saved."<<endl;  

}

bool Single_Process::FoundMappingFile(std::string & tempname) {
  char outname[100];
  sprintf(outname,"%s.map",(string("Process/")+m_ptypename+string("/")+m_name).c_str());
  if (IsFile(outname)) {
    msg.Debugging()<<"Mapping file already exists. Check identity."<<endl;
    ifstream from;
    from.open(outname);
    from>>tempname;
    msg.Tracking()<<"Found "<<tempname<<endl;
    return 1;
  }
  return 0;
}


void Single_Process::InitAnalysis(std::vector<APHYTOOLS::Primitive_Observable_Base *> _obs) {
  p_analysis = new APHYTOOLS::Primitive_Analysis(this->Name());
  for (int i=0;i<_obs.size();i++) {
    p_analysis->AddObservable(_obs[i]->GetCopy());
  }
  m_analyse  = 1;
}

/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/


bool Single_Process::SetUpIntegrator() {  
  msg.Debugging()<<"In  Single_Process::SetUpIntegrator(),  partner = "<<p_partner<<endl;
  p_sel->BuildCuts(p_cuts);
  if (m_nin==2) {
    if ( (p_fl[0].Mass() != p_isr->Flav(0).Mass()) ||
	 (p_fl[1].Mass() != p_isr->Flav(1).Mass()) ) p_isr->SetPartonMasses(p_fl);
  }
  if (CreateChannelLibrary()) {
    if (p_ps->CreateIntegrators()) return 1;
  }
  return 0;
}

bool Single_Process::CreateChannelLibrary()
{
  msg.Tracking()<<"Creating Multichannel for phasespace integration for "<<m_name<<endl;
  p_psgen     = new Phase_Space_Generator(m_nin,m_nout);
  bool newch  = 0;
  msg.Tracking()<<"Try to Construct("<<p_ps<<"->"<<p_ps->FSRIntegrator()<<","<<m_ptypename<<","<<m_libname<<")"<<endl;
  if (m_nin>1)  newch = p_psgen->Construct(p_ps->FSRIntegrator(),m_ptypename,m_libname,p_fl,this); 

  if (newch) {
    msg.Error()<<p_ps->NumberOfFSRIntegrators()<<" new Channels produced for "<<m_libname<<" ! "<<endl
	       <<"After program termination please enter \"make install\" and rerun !"<<endl;
    return 0;
  }
  else {
    msg.Tracking()<<"No new Channels produced for "<<m_libname<<" ! "<<endl
		  <<" added the following channels to the fs multi-channel : "<<endl;
    for (short int i=0;i<p_ps->NumberOfFSRIntegrators();i++)
      msg.Tracking()<<"     "<<(p_ps->FSRIntegrator()->Channel(i))->Name()<<endl;
    msg.Debugging()<<"Program continues."<<endl;
    return 1;
  }
}

/*------------------------------------------------------------------------------
  
  Process management
  
  ------------------------------------------------------------------------------*/

void Single_Process::Empty() {
  if (p_cuts)        { delete p_cuts; p_cuts = 0; }
  if (p_ps)          { delete p_ps; p_ps = 0; } 
  if (p_partner != this) {
    msg.Debugging()<<"Emptied "<<m_name<<" completely"<<endl;
    return;
  }
  msg.Debugging()<<"Emptied "<<m_name<<" partially"<<endl; 
}

void Single_Process::SetTotalXS(int _tables)  { 
  if (m_analyse) p_analysis->FinishAnalysis(m_resdir+string("/Tab")+m_name,_tables);
  m_totalxs  = m_totalsum/m_n; 
  m_totalerr = sqrt( (m_totalsumsqr/m_n - 
		    (AMATOOLS::sqr(m_totalsum)-m_totalsumsqr)/(m_n*(m_n-1.)) )  / m_n); 
  if (m_nin==2) {
    AORGTOOLS::msg.Events()<<"      xs for "<<m_name<<" : "
			   <<m_totalxs*AORGTOOLS::rpa.Picobarn()<<" pb"
			   <<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
			   <<"       max : "<<m_max<<endl;
  }
  if (m_nin==1) {
    AORGTOOLS::msg.Events()<<"      xs for "<<m_name<<" : "
			   <<m_totalxs<<" GeV"
			   <<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
			   <<"       max : "<<m_max<<endl;
 
  }
}

/*------------------------------------------------------------------------------

  Calculating total cross sections
  
  ------------------------------------------------------------------------------*/

bool Single_Process::CalculateTotalXSec(std::string _resdir) { 
  msg.Events()<<"In Single_Process::CalculateTotalXSec("<<_resdir<<") for "<<m_name<<endl; 
  
  string _name;
  double _totalxs,_totalerr,_max;
  char filename[100];
  sprintf(filename,"%s.xstotal",(_resdir+string("/")+m_name).c_str());
  if (_resdir!=string("")) {
    if (IsFile(filename)) {
      ifstream from;
      from.open(filename,ios::in);
      from>>_name>>_totalxs>>_totalerr>>_max;
      if (_name==m_name) {
	m_totalxs  = _totalxs;
	m_totalerr = _totalerr; 
	m_max      = _max;
      }
      msg.Events()<<"Found result : xs for "<<m_name<<" : "
		  <<m_totalxs*AORGTOOLS::rpa.Picobarn()<<" pb"
		  <<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
		  <<"       max : "<<m_max<<endl;
      from.close();
      if (m_totalxs>0.) return 1;
      return 0;
    }
  }
  m_totalxs = p_ps->Integrate();
  if (m_nin==2) m_totalxs /= AORGTOOLS::rpa.Picobarn();
  if (!(AMATOOLS::IsZero((m_n*m_totalxs-m_totalsum)/(m_n*m_totalxs+m_totalsum)))) {
    msg.Error()<<"Result of PS-Integrator and internal summation to not coincide!"<<endl
	       <<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<m_totalsum/m_n<<endl;
  }
  SetTotalXS(0);
  if (m_totalxs>0.) {
    if (_resdir!=string("")) {
      ofstream to;
      to.open(filename,ios::out);
      WriteOutXSecs(to);
      msg.Events()<<"Store result : xs for "<<m_name<<" : "
		  <<m_totalxs*AORGTOOLS::rpa.Picobarn()<<" pb"
		  <<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
		  <<"       max : "<<m_max<<endl;
      to.close();
    }
    return 1;
  }
  return 0;      
}

void Single_Process::WriteOutXSecs(ofstream & _to)    
{ 
  msg.Debugging()<<"Write out xsec for "<<m_name<<endl;
  _to<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<endl; 
}

bool Single_Process::Find(std::string _name,Process_Base *& _proc)  
{ 
  if (_name==m_name) {
    _proc = this;
    return 1;
  }
  return 0;
}

bool Single_Process::LookUpXSec(double ycut,bool calc,string obs) { 
  string filename = (m_resdir+string("/Tab")+m_name+string("/")+obs).c_str();
  if (IsFile(filename)) {
    Histogram * histo = new Histogram(filename);
    double    * res   = new double[histo->Depth()];
    histo->Extrapolate(ycut,res,1);
    m_totalxs = res[0];
    m_max     = res[1];
    msg.Events()<<m_name<<" : Set total xsec and max at ycut = "<<ycut
		<<" : "<<endl<<"   "<<m_totalxs<<" / "<<m_max<<endl;
    delete histo;
    delete res;

    if (calc) {
      p_ps->ReadIn(m_resdir+string("/MC_")+m_name);
      if (p_ps->BeamIntegrator() != 0) p_ps->BeamIntegrator()->Print();
      if (p_ps->ISRIntegrator()  != 0) p_ps->ISRIntegrator()->Print();
      if (p_ps->FSRIntegrator()  != 0) p_ps->FSRIntegrator()->Print();
    }
    return 1;
  }
  else {
    if (!calc)                return 0;
    if (!PrepareXSecTables()) return 0;
    return 1;
  }
}

bool Single_Process::PrepareXSecTables() { 
  msg.Events()<<"In Single_Process::PrepareXSecTables() for "<<m_name<<endl; 

  string filename = (m_resdir+string("/Tab")+m_name+string("dY_cut")).c_str();
  if (IsFile(filename)) {
    msg.Events()<<"Found "<<filename<<endl;
  }

  m_totalxs = p_ps->Integrate();
  if (m_nin==2) m_totalxs /= AORGTOOLS::rpa.Picobarn();

  if (!(AMATOOLS::IsZero((m_n*m_totalxs-m_totalsum)/(m_n*m_totalxs+m_totalsum)))) {
    msg.Error()<<"Result of PS-Integrator and internal summation to not coincide!"<<endl;
    msg.Error()<<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<m_totalsum/m_n<<endl;
  }
  SetTotalXS(1);
  p_ps->WriteOut(m_resdir+string("/MC_")+m_name);
  if (m_totalxs>0.) return 1;
  return 0;
}

void Single_Process::AddPoint(const double value) {
  m_n++;
  m_totalsum    += value;
  m_totalsumsqr += value*value;
  if (value>m_max) m_max = value;
  if (m_analyse) p_analysis->DoAnalysis(value*rpa.Picobarn());
}

double Single_Process::Differential(AMATOOLS::Vec4D* _moms) { return DSigma(_moms,0); }

double Single_Process::Differential2() { 
  if (p_isr->On()==0) return 0.;
  return DSigma2(); 
}

double Single_Process::DSigma(AMATOOLS::Vec4D* _moms,bool lookup)
{
  m_last = m_lastdxs = 0.;
  if (p_partner == this) {
    m_lastdxs = operator()(_moms);
  }
  else {
    if (lookup) m_lastdxs = p_partner->LastXS();
           else m_lastdxs = p_partner->operator()(_moms);
  }
  for (int i=0;i<m_nin+m_nout;i++) {
    if (_moms[i][0] < p_fl[i].PSMass()) return m_last = 0.;
  }
  if (m_lastdxs <= 0.)                  return m_lastdxs = m_last = 0.;
  if (m_nin==2) m_lastlumi = p_isr->Weight(p_flin);
          else  m_lastlumi = 1.;

  return m_last = m_Norm * m_lastdxs * m_lastlumi;
}

double Single_Process::DSigma2() { 
  if ((p_flin[0]==p_flin[1]) || (p_isr->On()==0) ) return 0.;
  if (p_partner == this) {
  }
  double tmp = m_Norm * m_lastdxs * p_isr->Weight2(p_flin); 
  m_last    += tmp;
  return tmp;
}

double Single_Process::operator()(AMATOOLS::Vec4D * mom)
{
  double M2 = 0.;

#ifndef Explicit_Pols   
  Vec4D gauge(sqrt(3.),1.,1.,1.);
  
  if (m_nout==4) gauge = mom[4];
  if (m_nout==5) gauge = mom[4];
  if (m_nout==6) gauge = mom[4];
  
  m_pol.Set_Gauge_Vectors(m_nin+m_nout,mom,gauge);
#endif
  p_BS->CalcEtaMu(mom);
  
  if (p_shand->Is_String()) {
    p_shand->Calculate();
    for (short int i=0;i<p_hel->Max_Hel();i++) {
      if (p_hel->On(i)) 
	M2 += p_ampl->Differential(i) * p_hel->Multiplicity(i) * p_hel->PolarizationFactor(i);
    }
  }
  else {
    for (short int i=0;i<p_hel->Max_Hel();i++)
      if (p_hel->On(i)) 
	M2 += p_ampl->Differential(i,(*p_hel)[i]) * p_hel->PolarizationFactor(i);
    p_shand->Complete(p_hel);
    p_ampl->ClearCalcList();
  }

  return M2 * sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
}


bool Single_Process::OneEvent()        { return (p_ps->OneEvent()); }
bool Single_Process::SameEvent()       { return (p_ps->SameEvent()); }
double Single_Process::WeightedEvent() { return (p_ps->WeightedEvent()); }




int Single_Process::NumberOfDiagrams() { 
  if (p_partner==this) return p_ampl->GetGraphNumber(); 
  return p_partner->NumberOfDiagrams();
}

Point * Single_Process::Diagram(int i) { return p_ampl->GetPointlist(i); } 



/*------------------------------------------------------------------------------
  
  Helpers
  
  ------------------------------------------------------------------------------*/

void Single_Process::PrintDifferential()
{
  if (!(AORGTOOLS::rpa.gen.Debugging())) return;
  AORGTOOLS::msg.Out()<<m_name<<" : "<<m_last<<" -> "
		      <<m_lastdxs<<" @ "<<m_lastlumi<<", "<<endl;
}
