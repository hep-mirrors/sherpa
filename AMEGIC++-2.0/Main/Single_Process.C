#include "Single_Process.H"
#include "Run_Parameter.H"
#include "ISR_Base.H"

#include "Running_AlphaS.H"
#include "Combined_Selector.H"

#include "Histogram.H"

#include "Random.H"
#include "prof.hh"

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>

using namespace AMEGIC;
using namespace PHASIC;
using namespace PDF;
using namespace BEAM;
using namespace ATOOLS;
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
			       int _kfactorscheme, int _scalescheme,double _scalefactor,double _scale,
			       Pol_Info * _pl,int _nex,Flavour * _ex_fl) :
  Process_Base(_nin,_nout,_fl,_isr,_beam,_gen_str,_orderQCD,_orderEW,
	       _scalescheme,_kfactorscheme,_scalefactor,_scale,_pl,_nex,_ex_fl),
  p_hel(0), p_BS(0), p_ampl(0), p_shand(0), p_partner(this), 
  m_helsample(false), m_inithelsample(false), m_throws(0), m_helresult(0.), m_helresult2(0.)
{
  m_newlib   = false;
  m_libnumb  = 0;
  m_save_max = 0.;
  GenerateNames(m_nin,p_flin,p_plin,m_nout,p_flout,p_plout,m_name,m_ptypename,m_libname);
  m_pslibname = m_libname;

  PolarizationNorm();
  if (_seldata) p_selector = new Combined_Selector(m_nin,m_nout,p_flavours,_seldata);
  else {
    if (m_nout>2)
      msg.Error()<<"Potential Error in Single_Process "<<m_name<<endl
		 <<"   No selection cuts specified. Init No_Selector !"<<endl;
    p_selector = new No_Selector();
  }

  double sum_massin = 0.,sum_massout = 0.;
  for (int i=0;i<m_nin;i++)  sum_massin  += p_flin[i].Mass();
  for (int i=0;i<m_nout;i++) sum_massout += p_flout[i].Mass();
  m_threshold = ATOOLS::Max(sum_massin,sum_massout);

  p_pshandler = new Phase_Space_Handler(this,p_isrhandler,p_beamhandler);
  
  // making directory
  if (m_gen_str>1) {
    unsigned int  mode_dir = 0755;
    mkdir((string("Process/")+m_ptypename).c_str(),mode_dir); 
  }
  msg.Tracking()<<"Initialized Single_Process : "<<m_name<<", "<<m_nvector<<", 1/norm = "<<1./m_Norm<<endl;
}


Single_Process::~Single_Process()
{
  if (p_hel)      {delete p_hel; p_hel=0;}
  if (p_BS)       {delete p_BS;   p_BS=0;}
  if (p_shand)    {delete p_shand;p_shand=0;}
  if (p_ampl)     {delete p_ampl; p_ampl=0;}
  if (p_psgen)    {delete p_psgen; p_psgen=0;}
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
  m_nvector = m_nin+m_nout+is_massless_pol+nmassive_pols;
  p_flavours   = new Flavour[m_nvector];
  p_pl   = new Pol_Info[m_nvector];
  p_b    = new int[m_nvector];
  for (short int i=0;i<m_nin;i++)             { p_flavours[i] = p_flin[i]       ; p_pl[i] = p_plin[i]       ; p_b[i] = -1; }
  for (short int i=m_nin;i<m_nin+m_nout;i++)  { p_flavours[i] = p_flout[i-m_nin]; p_pl[i] = p_plout[i-m_nin]; p_b[i] = 1; } 
  for (short int i=m_nin+m_nout;i<m_nvector;i++) { p_flavours[i] = Flavour(kf::pol); p_b[i]  = 1; }

  m_Norm = SymmetryFactors() * m_pol.Spin_Average(m_nin,p_flin);
#ifndef Explicit_Pols
  m_pol.Attach(m_nin+m_nout,p_flavours);
  m_Norm *= m_pol.Massive_Norm();

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
  
  double isrth = ATOOLS::Max(m_mass_in,m_mass_out);
  
  SetISRThreshold(isrth);

}

/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/



int Single_Process::InitAmplitude(Interaction_Model_Base * model,Topology* top,Vec4D *& _testmoms,
				  vector<Single_Process *> & links,vector<Single_Process *> & errs,
				  int & totalsize, int & procs)
{
  if (_testmoms==0) {
    string model_name = model->Name();
    if (model_name==string("ADD")) {
      double ms=model->ScalarConstant("M_s");
      double ecms=rpa.gen.Ecms();
      if (ms<0.95*ecms) {
	rpa.gen.SetEcms(0.5*ms);
	_testmoms = new Vec4D[m_nvector];
	p_pshandler->TestPoint(_testmoms);
	rpa.gen.SetEcms(ecms);    
	Vec4D * dummys = new Vec4D[m_nvector];
	p_pshandler->TestPoint(dummys);
	delete [] dummys;
      }
      else {
	_testmoms = new Vec4D[m_nvector];
	p_pshandler->TestPoint(_testmoms);
      }
    }
    else {
      _testmoms = new Vec4D[m_nvector];
      p_pshandler->TestPoint(_testmoms);
    }
  }
  if (p_momenta) { delete [] p_momenta; }
  p_momenta = new Vec4D[m_nvector]; 
  for (int i=0;i<m_nin+m_nout;i++) p_momenta[i] = _testmoms[i];

  p_hel    = new Helicity(m_nin,m_nout,p_flavours,p_pl);
  p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nvector,p_flavours,p_b);  
  p_shand  = new String_Handler(m_gen_str,p_BS);

  p_ampl   = new Amplitude_Handler(m_nin+m_nout,p_flavours,p_b,model,top,m_orderQCD,m_orderEW,
				   p_BS,p_shand);
  if (p_ampl->GetGraphNumber()==0) {
    msg.Tracking()<<"Single_Process::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return -1;
  }
  procs++;
  for (int j=0;j<links.size();j++) {
    if (p_ampl->CompareAmplitudes(links[j]->GetAmplitudeHandler())) {
      if (p_hel->Compare(links[j]->GetHelicity(),m_nin+m_nout)) {
	msg.Tracking()<<"Found compatible Process: "<<links[j]->Name()<<endl;

	if (!FoundMappingFile(m_libname,m_pslibname)) {
	  system((string("cp Process/")+m_ptypename+string("/")+links[j]->Name()+string(".map ")
		  +string("Process/")+m_ptypename+string("/")+Name()+string(".map")).c_str());
	  system((string("cp Process/")+m_ptypename+string("/")+links[j]->Name()+string(".col ")
		  +string("Process/")+m_ptypename+string("/")+Name()+string(".col")).c_str());
	}
	
	p_partner = links[j];
	Minimize();
	return 1;
      }
    }
  }
  p_ampl->CompleteAmplitudes(m_nin+m_nout,p_flavours,p_b,&m_pol,
			     top,p_BS,m_ptypename+string("/")+m_name);

  m_pol.Add_Extern_Polarisations(p_BS,p_flavours,p_hel);
  p_BS->Initialize();

  switch (Tests()) {
  case 2 : 
    for (int j=0;j<links.size();j++) {
      if (ATOOLS::IsEqual(links[j]->Result(),Result())) {
	msg.Tracking()<<"Test : 2.  Can map "<<m_name<<" on "<<links[j]->Name()<<endl;
	p_partner = links[j];
	break;
      } 
    }
    if (p_partner==this) {
      links.push_back(this);
      totalsize++;
    }
    Minimize();
    return 1;
  case 1 :
    for (int j=0;j<links.size();j++) {
      if (ATOOLS::IsEqual(links[j]->Result(),Result())) {
	msg.Tracking()<<"Test : 1.  Can map "<<m_name<<" on "<<links[j]->Name()<<endl;
	p_partner = links[j];
	m_pslibname = links[j]->PSLibName();
	break;
      } 
    }
    if (p_partner==this) links.push_back(this);
    
    if (CheckLibraries()) return 1;
    for (int j=0;j<links.size();j++) {
      //if (ATOOLS::IsEqual(links[j]->Result(),Result())) {
      if (links[j]->NewLibs()) {
	if (CheckStrings(links[j])) return 1;	
      }      
    }
    if (p_partner!=this) links.push_back(this);
    
    if (m_gen_str<2) return 1;
    totalsize++;
    WriteLibrary();
    if (p_partner==this && Result()>0.) SetUpIntegrator();
    return 0;
  case -3: return -3;
  default :
    msg.Error()<<"Error in Single_Process::InitAmplitude : Failed for "<<m_name<<"."<<endl;
    errs.push_back(this);
    return 1;
  }
}


int Single_Process::InitAmplitude(Interaction_Model_Base * model,Topology * top)
{
  if (p_momenta) { delete [] p_momenta; }
  p_momenta   = new Vec4D[m_nvector]; 
  p_pshandler->TestPoint(p_momenta);

  p_hel    = new Helicity(m_nin,m_nout,p_flavours,p_pl);
  p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nvector,p_flavours,p_b);  
  p_shand  = new String_Handler(m_gen_str,p_BS);
  p_ampl   = new Amplitude_Handler(m_nin+m_nout,p_flavours,p_b,model,top,m_orderQCD,m_orderEW,
				   p_BS,p_shand);
  if (p_ampl->GetGraphNumber()==0) {
    msg.Tracking()<<"Single_Process::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return -1;
  }
  p_ampl->CompleteAmplitudes(m_nin+m_nout,p_flavours,p_b,&m_pol,
			     top,p_BS,m_ptypename+string("/")+m_name);

  m_pol.Add_Extern_Polarisations(p_BS,p_flavours,p_hel);
  p_BS->Initialize();

  switch (Tests()) {
  case 2 : return 1;
  case 1 : 
    if (CheckLibraries()) return 1;
    WriteLibrary();
    return 0;
  default :
    msg.Error()<<"Error in Single_Process::InitAmplitude : Failed for "<<m_name<<"."<<endl;
    return -2;
  }
}


void Single_Process::InitDecay(Topology* top) { }

int Single_Process::Tests() {
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
  else p_shand->Initialize(p_ampl->GetRealGraphNumber(),p_hel->MaxHel());

  p_ampl->SetStringOff();

  double M2 = 0.;
  double helvalue;

  if (gauge_test) {
    m_pol.Set_Gauge_Vectors(m_nin+m_nout,p_momenta,Vec4D(sqrt(3.),1.,1.,-1.));
    p_BS->Setk0(0);
    p_BS->CalcEtaMu(p_momenta);  
    p_BS->InitGaugeTest(.9);

    msg.Debugging()<<number<<" :";ATOOLS::msg.Debugging().flush();
    for (short int i=0;i<p_hel->MaxHel();i++) { 
      if (p_hel->On(i)) {
	helvalue = p_ampl->Differential(i,(*p_hel)[i])*p_hel->PolarizationFactor(i); 
	//cout<<i<<". :"<<helvalue<<endl;
	M2      +=  helvalue;
	msg.Debugging()<<"*";msg.Debugging().flush();
      } 
      else {
	msg.Debugging()<<"0";msg.Debugging().flush();
      }
    }
    msg.Debugging()<<endl;
    M2     *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS));
    m_iresult  = M2;
  }
  p_ampl->ClearCalcList();
  // To prepare for the string test.
  p_ampl->SetStringOn();
  (p_shand->Get_Generator())->Reset(1);
  /* ---------------------------------------------------
     
  First test : gauge test
  
  --------------------------------------------------- */
#ifndef Explicit_Pols 
  Vec4D gauge(sqrt(3.),1.,1.,1.);
  if (m_nout==4) gauge = p_momenta[4];
  if (m_nout==5) gauge = p_momenta[4];
  //?????????
  if (m_nout==6) gauge = p_momenta[4];  
  m_pol.Reset_Gauge_Vectors(m_nin+m_nout,p_momenta,gauge);
#else
  p_BS->Setk0(1);
#endif
  p_BS->CalcEtaMu(p_momenta);
  number++;
  msg.Debugging()<<number<<" :";

  if (!gauge_test) p_ampl->SetStringOff();  //second test without string production 

  double M2g = 0.;
  double * M_doub = new double[p_hel->MaxHel()];

  for (short int i=0;i<p_hel->MaxHel();i++) { 
    if (p_hel->On(i)) {
      M_doub[i]  = p_ampl->Differential(i,(*p_hel)[i])*p_hel->PolarizationFactor(i);  
      //cout<<i<<". :"<<M_doub[i]<<endl;
      M2g       += M_doub[i];
      msg.Debugging()<<"*";msg.Debugging().flush();
    }
  }
  msg.Debugging()<<endl;

  //shorten helicities
  for (short int i=0;i<p_hel->MaxHel();i++) {
    if (M_doub[i]==0. || M_doub[i]/M2g<(ATOOLS::Accu()*1.e-2)) {
      p_hel->SwitchOff(i);
      msg.Debugging()<<"Switch off zero helicity "<<i<<" : "
		     <<M_doub[i]<<"/"<<M_doub[i]/M2g<<endl;
    }
  }

  M2g    *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS));
  m_iresult  = M2g;

  p_ampl->ClearCalcList();  
  p_ampl->FillCoupling(p_shand);
  p_ampl->KillZList();  
  p_BS->StartPrecalc();

  if (gauge_test) {
    if (!ATOOLS::IsEqual(M2,M2g)) {
      msg.Out()<<"Gauge(1): "<<abs(M2)<<endl
	       <<"Gauge(2): "<<abs(M2g)<<endl;
      msg.Out()<<"WARNING:  Gauge test not satisfied: "
	       <<M2<<" vs. "<<M2g<<" : "<<dabs(M2/M2g-1.)*100.<<"%"<<endl;
    }
    else {
      msg.Debugging()<<"Gauge(1): "<<abs(M2)<<endl
		     <<"Gauge(2): "<<abs(M2g)<<endl;
      if (M2g!=0.)
	msg.Debugging()<<"Gauge test: "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
      else
	msg.Debugging()<<"Gauge test: "<<0.<<"%"<<endl;
    }
  }
  else {
    delete[] M_doub;
    number++;
    if (p_shand->SearchValues(m_gen_str,testname,p_BS)) {
      p_shand->Initialize(p_ampl->GetRealGraphNumber(),p_hel->MaxHel());
      (p_shand->Get_Generator())->Reset();
  
      M2 = operator()(p_momenta);
      gauge_test = string_test = 0;
    }
    if (!ATOOLS::IsEqual(M2,M2g)) {
      msg.Out()<<"Mapping file(1) : "<<abs(M2)<<endl
	       <<"Original    (2) : "<<abs(M2g)<<endl
	       <<"Cross check (T) : "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
      msg.Out()<<"WARNING: Library cross check not satisfied: "
	       <<M2<<" vs. "<<M2g<<"  difference:"<<abs(M2/M2g-1.)*100.<<"%"<<endl;
      if (abs(M2/M2g-1.)>rpa.gen.Accu()) return 0;
      msg.Out()<<"         assuming numerical reasons, continuing "<<endl;
    } 
    else {
      msg.Debugging()<<"Mapping file(1) : "<<abs(M2)<<endl
		     <<"Original    (2) : "<<abs(M2g)<<endl;
      if (M2g!=0.)
	msg.Debugging()<<"Cross check (T) : "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
      else {
	msg.Debugging()<<"Cross check (T) : "<<0.<<"%"<<endl;
	m_libname    = testname;
	return -3;
      }
    }

    m_libname    = testname;
    return 2;
  }

  /* ---------------------------------------------------
     
     Second test : string test

     --------------------------------------------------- */

  if (string_test) {
    //String-Test
    for (short int i=0;i<p_hel->MaxHel();i++) {
      if (p_hel->On(i)) {
	for (short int j=i+1;j<p_hel->MaxHel();j++) {
	  if (p_hel->On(j)) {
	    if (ATOOLS::IsEqual(M_doub[i],M_doub[j])) {
	      msg.Debugging()<<"Mapping equal helicities "<<j<<" -> "<<i<<endl;
	      p_hel->SwitchOff(j);
	      p_hel->SetPartner(i,j);
	      p_hel->IncMultiplicity(i);
	    }
	  }
	}
      }
    }
    delete[] M_doub;
    PROFILE_LOCAL("Shand.Complete()");
    p_shand->Complete(p_hel);

    if (p_shand->Is_String()) {
      double  M2S = 0.;
      p_shand->Calculate();
      
      msg.Debugging()<<"3:";msg.Debugging().flush();
      for (short int i=0;i<p_hel->MaxHel();i++) {
	if (p_hel->On(i)) {
	  msg.Debugging()<<"*";msg.Debugging().flush();
	  M2S += p_ampl->Differential(i)*p_hel->PolarizationFactor(i)*p_hel->Multiplicity(i);
	}
      }
      msg.Debugging()<<endl;
      M2S *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS));
      if (!ATOOLS::IsEqual(M2g,M2S)) {
	msg.Out()<<"WARNING: String test not satisfied: "
		 <<M2g<<" vs. "<<M2S<<"  difference:"<<abs(M2g/M2S-1.)*100.<<"%"<<endl;
	if (abs(M2g/M2S-1.)>rpa.gen.Accu()) return 0;
	msg.Out()<<"         assuming numerical reasons, continuing "<<endl;
      }
      else {
	if (M2S!=0.)
	  msg.Debugging()<<"String test: "<<abs(M2g/M2S-1.)*100.<<"%"<<endl;
	else
	  msg.Debugging()<<"String test: "<<0.<<"%"<<endl;
      
      }
      return 1;
    }
    return 1;
  }
  delete[] M_doub;
  return 0;
}

int Single_Process::CheckLibraries() {
  if (m_gen_str==0) return 1;
  if (p_shand->IsLibrary()) return 1;

  char help[20];
  String_Handler * shand1;
  shand1      = new String_Handler(p_shand->Get_Generator());
  
  m_libnumb  = 0;
  string proc = string("Process/")+m_ptypename+string("/V");
  string testname;
  double M2s, helvalue;

  for (;;) {
    sprintf(help,"%i",m_libnumb);
    testname  = CreateLibName()+string("_")+string(help);
    if (shand1->SearchValues(m_gen_str,testname,p_BS)) {

      shand1->Calculate();
      
      M2s = 0.;
      ATOOLS::msg.Debugging()<<"Check "<<m_libnumb<<" :";ATOOLS::msg.Debugging().flush();
      for (short int i=0;i<p_hel->MaxHel();i++) {
	helvalue = p_ampl->Differential(shand1,i) * p_hel->PolarizationFactor(i) *
	  p_hel->Multiplicity(i);
	M2s     += helvalue;
	msg.Debugging()<<"*";ATOOLS::msg.Debugging().flush();
      }
      msg.Debugging()<<endl;
      M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS));
      if (Result()!=0.)
	msg.Debugging()<<"Cross check (1): "<<abs(M2s/Result()-1.)*100.<<"%"<<"  : "
		       <<M2s<<"/"<<Result()<<endl;
      else msg.Debugging()<<"Cross check (1): "<<M2s<<"/"<<Result()<<endl;
      if (ATOOLS::IsEqual(M2s,Result())) {
	msg.Tracking()<<"Found a suitable Library."<<endl;
	m_libname = testname;
	m_pslibname = testname;
	if (shand1) { delete shand1; shand1 = 0; }
	//Clean p_shand!!!!
	Minimize();
	CreateMappingFile();
	return 1;
      }
    } 
    else break;
    ++m_libnumb;
  }
  if (shand1) { delete shand1; shand1 = 0; }
  return 0;
}

int Single_Process::CheckStrings(Single_Process* tproc)
{
  if (tproc->LibName().find(CreateLibName())!=0) return 0;

  String_Handler * shand1;
  shand1 = new String_Handler(p_shand->Get_Generator(),
			      (tproc->GetStringHandler())->GetSKnots());
  (shand1->Get_Generator())->ReplaceZXlist((tproc->GetStringHandler())->Get_Generator());
  double M2s, helvalue;
  shand1->Calculate();

  M2s = 0.;
  ATOOLS::msg.Debugging()<<"Check "<<tproc->Name()<<" :";ATOOLS::msg.Debugging().flush();
  for (short int i=0;i<p_hel->MaxHel();i++) {
    helvalue = p_ampl->Differential(shand1,i) * p_hel->PolarizationFactor(i) *
      p_hel->Multiplicity(i);
    M2s     += helvalue;
    msg.Debugging()<<"*";ATOOLS::msg.Debugging().flush();
  }
  msg.Debugging()<<endl;
  M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS));
  if (Result()!=0.)
    msg.Debugging()<<"Cross check (2): "<<abs(M2s/Result()-1.)*100.<<"%"<<"  : "
		   <<M2s<<"/"<<Result()<<endl;
  else msg.Debugging()<<"Cross check (2): "<<M2s<<"/"<<Result()<<endl;  
  (shand1->Get_Generator())->ReStore();
  delete shand1;

  if (ATOOLS::IsEqual(M2s,Result())) {
    msg.Tracking()<<"Found a suitable string."<<endl;
    m_libname = tproc->LibName();
    m_pslibname = tproc->PSLibName();
    Minimize();
    CreateMappingFile();
    return 1;
  }
  return 0;
}
  
void Single_Process::WriteLibrary() 
{
  if (m_gen_str<2) return;
  char help[20];
  string testname;
  for (;;) {
    sprintf(help,"%i",m_libnumb);
    testname    = CreateLibName()+string("_")+string(help);
    if (!(IsFile(string("Process/")+m_ptypename+string("/")+testname+string("/V.H")))) break;
    ++m_libnumb;
  }
  m_libname = testname;
  if (p_partner==this) m_pslibname = m_libname;
                  else m_pslibname = p_partner->PSLibName();
  msg.Debugging()<<"Write Library for "<<m_name<<" = "<<m_libname<<", = case."<<endl;
  int  mode_dir = 448;
  msg.Debugging()<<" m_ptypename = "<<m_ptypename<<endl<<" m_libname = "<<m_libname<<endl;
  mkdir((string("Process/")+m_ptypename+string("/")+m_libname).c_str(),mode_dir); 
  p_shand->Output(p_hel,m_ptypename+string("/")+m_libname);
  CreateMappingFile();
  m_newlib=true;
}

std::string  Single_Process::CreateLibName()
{
  string name=m_ptypename;
  char help[20];
  sprintf(help,"%i",p_ampl->GetGraphNumber());
  name += string("_");
  name += string(help);
  sprintf(help,"%i",p_shand->NumberOfCouplings());
  name += string("_");
  name += string(help);
  sprintf(help,"%i",p_shand->NumberOfZfuncs());
  name += string("_");
  name += string(help);
  sprintf(help,"%i",p_hel->MaxHel());
  name += string("_");
  name += string(help);
  sprintf(help,"%i",p_BS->MomlistSize());
  name += string("_");
  name += string(help);
  return name;
}

void Single_Process::CreateMappingFile() {
  if (m_gen_str<2) return;
  char outname[100];
  sprintf(outname,"%s.map",(string("Process/")+m_ptypename+string("/")+m_name).c_str());
  if (IsFile(outname)) {
    string MEname,PSname;
    FoundMappingFile(MEname,PSname);
    if (MEname != m_libname || PSname != m_pslibname) {
      msg.Error()<<"In Single_Process::CreateMappingFile() : Files do not coincide. Maybe changed input data ?"<<endl;
      abort();
    }
    else return;
  }

  std::ofstream to;
  to.open(outname,ios::out);
  to<<"ME: "<<m_libname<<endl
    <<"PS: "<<m_pslibname<<endl;
  to.close();
}

bool Single_Process::FoundMappingFile(std::string & MEname, std::string & PSname) {
  char outname[100];
  char buffer[100];
  string buf;
  int pos;
  sprintf(outname,"%s.map",(string("Process/")+m_ptypename+string("/")+m_name).c_str());
  if (IsFile(outname)) {
    ifstream from;
    from.open(outname);
    from.getline(buffer,100);
    buf = string(buffer);
    pos = buf.find(string("ME:"));
    if (pos==-1) MEname = PSname = buf;
    else {
      MEname = buf.substr(pos+4);
      from.getline(buffer,100);
      buf = string(buffer);
      pos = buf.find(string("PS:"));
      if (pos==-1) PSname = MEname;
      else PSname = buf.substr(pos+4);
      if (PSname==string("")) PSname = MEname;
    }
    return 1;
  }
  return 0;
}



/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/


bool Single_Process::SetUpIntegrator() 
{  
  if (m_nin==2) {
    if ( (p_flavours[0].Mass() != p_isrhandler->Flav(0).Mass()) ||
	 (p_flavours[1].Mass() != p_isrhandler->Flav(1).Mass()) ) p_isrhandler->SetPartonMasses(p_flavours);
    if (CreateChannelLibrary()) {
      if (p_pshandler->CreateIntegrators()) return 1;
    }
  }
  if (m_nin==1) return p_pshandler->CreateIntegrators();
  return 0;
}

bool Single_Process::CreateChannelLibrary()
{
  p_psgen     = new Phase_Space_Generator(m_nin,m_nout);
  bool newch  = 0;
  if (m_nin>1)  newch = p_psgen->Construct(p_pshandler->FSRIntegrator(),m_ptypename,m_pslibname,p_flavours,this); 

  if (newch>0) {
    msg.Error()<<p_pshandler->NumberOfFSRIntegrators()<<" new Channels produced for "<<m_pslibname<<" ! "<<endl
	       <<"After program termination please enter \"make install\" and rerun !"<<endl;
    return 0;
  }
  else {
    return 1;
  }
}

/*------------------------------------------------------------------------------
  
  Process management
  
  ------------------------------------------------------------------------------*/
void Single_Process::Minimize()
{
  if (p_partner==this) return;
  if (p_hel)      {delete p_hel; p_hel=0;}
  if (p_BS)       {delete p_BS;   p_BS=0;}
  if (p_shand)    {delete p_shand;p_shand=0;}
  if (p_ampl)     {delete p_ampl; p_ampl=0;}
  if (p_psgen)    {delete p_psgen; p_psgen=0;}

  if (p_momenta)     { delete [] p_momenta;  p_momenta     = 0; }
  if (p_selector)      { delete p_selector;      p_selector      = 0; }
  if (p_pshandler)       { delete p_pshandler;       p_pshandler       = 0; }
  p_selected=p_partner;
}

void Single_Process::Empty() {
  if (p_pshandler)          { delete p_pshandler; p_pshandler = 0; } 
  if (p_partner != this) {
    return;
  }
}

void Single_Process::SetTotalXS(int _tables)  { 
  if (_tables!=2) {
    m_totalxs  = m_totalsum/m_n; 
    m_totalerr = sqrt( (m_totalsumsqr/m_n - 
			(ATOOLS::sqr(m_totalsum)-m_totalsumsqr)/(m_n*(m_n-1.)) )  / m_n); 
    if ((m_nin==1 && m_nout==2) || m_n==1) m_totalerr = 0.;
  }
  else {
    //   _tables==2  means  check xs with sum of subprocesses
    //   nothing to do for a Single_Process
  }
  if (m_nin==2) {
    msg.Events()<<"      xs for "<<om::bold<<m_name<<om::reset<<" : "
		<<om::blue<<om::bold<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"<<om::reset
		<<" +/- "<<om::reset<<om::blue<<m_totalerr/m_totalxs*100.<<" %,"<<om::reset<<endl
		<<"       max : "<<m_max<<endl
		<<om::bold<<"   exp. eff: "<<om::red<<(100.*m_totalxs/m_max)<<" %"<<om::reset<<endl;
  }
  if (m_nin==1) {
    msg.Events()<<"      width for "<<m_name<<" : "
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
  double _totalxs,_totalerr,_max,sum,sqrsum;
  long int n;
  char filename[100];
  sprintf(filename,"%s.xstotal",(_resdir+string("/")+m_name).c_str());
  if (_resdir!=string("")) {
    if (IsFile(filename)) {
      ifstream from;
      from.open(filename,ios::in);
      from>>_name>>_totalxs>>_max>>_totalerr>>sum>>sqrsum>>n;
      if (_name==m_name) {
	m_totalxs  = _totalxs;
	m_totalerr = _totalerr; 
	m_max      = _max;
	m_n        = n;
	m_totalsum = sum;
	m_totalsumsqr = sqrsum;
      }
      msg.Events()<<"Found result : xs for "<<m_name<<" : "
		  <<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
		  <<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
		  <<"       max : "<<m_max<<endl;
      from.close();
      p_pshandler->ReadIn(_resdir+string("/MC_")+m_name);
      if (p_pshandler->BeamIntegrator() != 0) p_pshandler->BeamIntegrator()->Print();
      if (p_pshandler->ISRIntegrator()  != 0) p_pshandler->ISRIntegrator()->Print();
      if (p_pshandler->FSRIntegrator()  != 0) p_pshandler->FSRIntegrator()->Print();
      p_pshandler->InitIncoming();
    }
  }
  m_resultpath=_resdir;
  m_resultfile=filename;
  ATOOLS::Exception_Handler::AddTerminatorObject(this);
  m_totalxs = p_pshandler->Integrate();
  if (m_nin==2) m_totalxs /= ATOOLS::rpa.Picobarn();
  if (!(ATOOLS::IsZero((m_n*m_totalxs-m_totalsum)/(m_n*m_totalxs+m_totalsum)))) {
    msg.Error()<<"Result of PS-Integrator and internal summation to not coincide!"<<endl
	       <<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<m_totalsum/m_n<<endl;
  }
  SetTotalXS(0);
  if (m_totalxs>=0.) {
    if (_resdir!=string("")) {
      std::ofstream to;
      to.open(filename,ios::out);
      WriteOutXSecs(to);
      msg.Events()<<"Store result : xs for "<<m_name<<" : "
		  <<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
		  <<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
		  <<"       max : "<<m_max<<endl;
      p_pshandler->WriteOut(_resdir+string("/MC_")+m_name);
      to.close();
    }
    ATOOLS::Exception_Handler::RemoveTerminatorObject(this);
    return 1;
  }
  ATOOLS::Exception_Handler::RemoveTerminatorObject(this);
  return 0;      
}

void Single_Process::PrepareTerminate()
{
  if (m_resultpath.length()==0 && m_resultfile.length()==0) return;
  std::ofstream to;
  to.open(m_resultfile.c_str(),ios::out);
  WriteOutXSecs(to);
  msg.Events()<<"Store result : xs for "<<m_name<<" : "
	      <<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
	      <<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
	      <<"       max : "<<m_max<<endl;
  p_pshandler->WriteOut(m_resultpath+string("/MC_")+m_name);
  to.close();
}

void Single_Process::WriteOutXSecs(std::ofstream & _to)    
{ 
  _to<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<" "
     <<m_totalsum<<" "<<m_totalsumsqr<<" "<<m_n<<endl; 
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
      p_pshandler->ReadIn(m_resdir+string("/MC_")+m_name);
      if (p_pshandler->BeamIntegrator() != 0) p_pshandler->BeamIntegrator()->Print();
      if (p_pshandler->ISRIntegrator()  != 0) p_pshandler->ISRIntegrator()->Print();
      if (p_pshandler->FSRIntegrator()  != 0) p_pshandler->FSRIntegrator()->Print();
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

  m_totalxs = p_pshandler->Integrate();
  if (m_nin==2) m_totalxs /= ATOOLS::rpa.Picobarn();

  if (!(ATOOLS::IsZero((m_n*m_totalxs-m_totalsum)/(m_n*m_totalxs+m_totalsum)))) {
    msg.Error()<<"Result of PS-Integrator and internal summation to not coincide!"<<endl;
    msg.Error()<<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<m_totalsum/m_n<<endl;
  }
  SetTotalXS(1);
  p_pshandler->WriteOut(m_resdir+string("/MC_")+m_name);
  if (m_totalxs>0.) return 1;
  return 0;
}

void Single_Process::AddPoint(const double value) {
  m_n++;
  m_totalsum    += value;
  m_totalsumsqr += value*value;
  if (value>m_max) m_max = value;
  int iter=1;
  if (m_n<200000)        iter =  40000;
  else if  (m_n<400000)  iter = 200000;
  // ==== uncomment for better unweighting efficiency ===
  //   if (iter!=1 && m_n%iter==0) {
  //     m_max = m_save_max;
  //     m_save_max = 0.;
  //   }
  if (value>m_save_max) m_save_max = value;
}

double Single_Process::Differential(const ATOOLS::Vec4D* _moms) { return DSigma(_moms,0); }

double Single_Process::Differential2() { 
  if (p_isrhandler->On()==0) return 0.;
  return DSigma2(); 
}


double Single_Process::DSigma(const ATOOLS::Vec4D* _moms,bool lookup)
{
  m_last = m_lastdxs = 0.;
  for (int i=0;i<m_nin+m_nout;i++) {
    if (_moms[i][0] < p_flavours[i].PSMass()) return m_last = 0.;
  }
  if (p_partner == this) {
    if (m_helsample) {
      if (!m_inithelsample) InitializeHelicityWeights();
      m_lastdxs = operator()(_moms,SelectedHelicity());
    }
    else m_lastdxs = operator()(_moms);
  }
  else {
    if (lookup) m_lastdxs = p_partner->LastXS();
           else m_lastdxs = p_partner->operator()(_moms);
  }
  if (m_lastdxs <= 0.) return m_lastdxs = m_last = 0.;
  if (m_nin==2) {
    m_lastlumi = p_isrhandler->Weight(p_flin);
    int    pols[2] = {p_pl[0].type[0],p_pl[1].type[0]};
    double dofs[2] = {p_pl[0].factor[0],p_pl[1].factor[0]};
    if (p_pl[0].num>1) pols[0] = 99;
    if (p_pl[1].num>1) pols[1] = 99;
    m_lastlumi *= p_beamhandler->Weight(pols,dofs);
  }
  else  m_lastlumi = 1.;

  return m_last = m_Norm * m_lastdxs * m_lastlumi;
}

double Single_Process::DSigma2() { 
  if ((p_flin[0]==p_flin[1]) || (p_isrhandler->On()==0) ) return 0.;
  if (p_partner == this) {
  }
  double tmp = m_Norm * m_lastdxs * p_isrhandler->Weight2(p_flin); 
  m_last    += tmp;
  return tmp;
}

double Single_Process::operator()(const ATOOLS::Vec4D * mom)
{
  double M2 = 0.;

#ifndef Explicit_Pols   
  Vec4D gauge(sqrt(3.),1.,1.,1.);
  
  if (m_nout==4) gauge = mom[4];
  if (m_nout==5) gauge = mom[4];
  if (m_nout==6) gauge = mom[4];
  
  m_pol.Set_Gauge_Vectors(m_nin+m_nout,mom,gauge);
#endif
  ATOOLS::Vec4D *cpmom = new ATOOLS::Vec4D[m_nvector];
  for (size_t i=0;i<m_nvector;++i) cpmom[i]=mom[i];
  p_BS->CalcEtaMu(cpmom);
  delete [] cpmom;

  double helvalue;
  if (p_shand->Is_String()) {
    p_shand->Calculate();
    for (short int i=0;i<p_hel->MaxHel();i++) {
      if (p_hel->On(i)) {
	helvalue = p_ampl->Differential(i) * p_hel->Multiplicity(i) * p_hel->PolarizationFactor(i);
	M2      += helvalue;
      }
    }
  }
  else {
    for (short int i=0;i<p_hel->MaxHel();i++) {
      if (p_hel->On(i)) {
	helvalue = p_ampl->Differential(i,(*p_hel)[i]) * p_hel->PolarizationFactor(i);
	M2 += helvalue;
      }
    }
    p_shand->Complete(p_hel);
    p_ampl->ClearCalcList();
  }

  return M2 * sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS));
}


bool   Single_Process::OneEvent(double _mass) { return p_pshandler->OneEvent(_mass); }
bool   Single_Process::SameEvent()            { return p_pshandler->SameEvent(); }
ATOOLS::Blob_Data_Base * Single_Process::WeightedEvent()     { return p_pshandler->WeightedEvent(); }
ATOOLS::Blob_Data_Base * Single_Process::SameWeightedEvent() { return p_pshandler->SameWeightedEvent(); }




int Single_Process::NumberOfDiagrams() { 
  if (p_partner==this) return p_ampl->GetGraphNumber(); 
  return p_partner->NumberOfDiagrams();
}

Point * Single_Process::Diagram(int i) { 
  if (p_partner==this) return p_ampl->GetPointlist(i); 
  return p_partner->Diagram(i);
} 



/*------------------------------------------------------------------------------
  
  Helpers
  
  ------------------------------------------------------------------------------*/

void Single_Process::PrintDifferential()
{
  if (!(ATOOLS::msg.LevelIsDebugging())) return;
  ATOOLS::msg.Out()<<m_name<<" : "<<m_last<<" -> "
		      <<m_lastdxs<<" @ "<<m_lastlumi<<", "<<endl;
}


/*------------------------------------------------------------------------------
  
  Stuff to sample over helicities

  ------------------------------------------------------------------------------*/

void Single_Process::InitializeHelicityWeights()
{
  int activehels = 0, active = 0;
  for (int i=0;i<p_hel->MaxHel();i++) {
    if (p_hel->On(i)) {
      active++,
      activehels += p_hel->Multiplicity(i);
    }
  }

  double alpha_start = 1./activehels;
  double alpha;
  for (int i=0;i<p_hel->MaxHel();i++) {
    if (p_hel->On(i)) {
      alpha = alpha_start * p_hel->Multiplicity(i);
      m_helnumbers.push_back(i);
      m_helalphas.push_back(alpha);
      m_helalphasaves.push_back(alpha);
      m_helresults.push_back(0.);
      m_helresults2.push_back(0.);
      m_helthrows.push_back(0);
    }
  }
  m_throws        = 0;
  m_helnumber     = m_helnumbers.size();
  m_inithelsample = true;
  if (msg.Level()>1) {
    msg.Out()<<"Initialize sampling over helicities."<<std::endl
	     <<"   Found "<<active<<" active helicities with the following weights: "<<std::endl;
    for (int i=0;i<m_helnumber;i++) {
      msg.Out()<<"   "<<i<<"("<<m_helnumber<<") : "<<m_helnumbers[i]<<" * "
	       <<p_hel->Multiplicity(m_helnumbers[i])<<" -> "<<m_helalphas[i]<<std::endl;
    }
  }
}

const int Single_Process::SelectedHelicity()
{
  if (m_throws>0 && (!(m_throws%1000))) OptimizeHelicityWeights();
  double disc = ran.Get();
  int hel;
  for (hel=0;hel<m_helnumber;hel++) {
    disc -= m_helalphas[hel];
    if (disc<=0.) break;
  }
  if (hel>=m_helnumber) {
    msg.Error()<<"Warning in Single_Process::SelectedHelicity() after "<<m_throws<<std::endl;
    hel = m_helnumber-1;
  }
  return hel;
}

void Single_Process::OptimizeHelicityWeights()
{
  short int i;

  double aptot = 0.;
  for (i=0;i<m_helnumber;i++) {
    if (m_helthrows[i]>0) m_helresults2[i]  = sqrt(m_helresults[i]/m_helthrows[i]);
                     else m_helresults2[i]  = 0.;
    aptot              += m_helalphas[i]*m_helresults2[i];
  }

  double s1x = 0.;  
  for (i=0;i<m_helnumber;i++) {
    if (dabs(aptot-m_helresults2[i])>s1x) s1x = dabs(aptot-m_helresults2[i]);
    if (m_helthrows[i]>0)        m_helalphas[i] *= m_helresults2[i]/aptot;
    if (m_helalphas[i] < 1.e-8 ) m_helalphas[i] = 0.;
  }

  double norm = 0;
  for (i=0;i<m_helnumber;i++) norm += m_helalphas[i];
  for (i=0;i<m_helnumber;i++) {
    m_helresults[i] = 0.;
    m_helthrows[i]  = 0;
    m_helalphas[i] /= norm;
  }

  msg.Tracking()<<"After Optimize Helicity Weights for "<<m_name<<" after "<<m_throws<<" throws: "<<std::endl; 
  for (i=0;i<m_helnumber;i++) {
    msg.Tracking()<<i<<" th helicity: "<<m_helnumbers[i]<<"("<<p_hel->Multiplicity(m_helnumbers[i])
		  <<"), alpha = "<<m_helalphas[i]<<std::endl;
  }

  double variance = (m_throws*m_helresult2)/((m_throws-1)*ATOOLS::sqr(m_helresult)) - 1./(m_throws-1);
  if (variance>0.) variance = m_helresult/m_throws * sqrt(variance);
  else {
    msg.Error()<<"Negative variance."<<std::endl;
    variance = m_helresult/m_throws * sqrt(-variance);
  }

  msg.Tracking()<<"S1X: "<<s1x<<", variance : "<<variance<<std::endl
		<<"result,result2,n: "<<m_helresult<<", "
		<<m_helresult2<<", "<<m_throws<<std::endl
		<<"-----------------------------------------------"<<endl;
}


double Single_Process::operator()(const ATOOLS::Vec4D * mom,const int hel)
{
  if (!p_shand->Is_String()) {
    msg.Error()<<"Error in Single_Process::operator()(ATOOLS::Vec4D * p,int hel)"<<std::endl
	       <<"   Sampling over helicities in processes implemented only for libs."<<std::endl
	       <<"   Will abort the run. Check for libraries."<<std::endl;
    abort();
  }
#ifndef Explicit_Pols   
  Vec4D gauge(sqrt(3.),1.,1.,1.);
  
  if (m_nout==4) gauge = mom[4];
  if (m_nout==5) gauge = mom[4];
  if (m_nout==6) gauge = mom[4];
  
  m_pol.Set_Gauge_Vectors(m_nin+m_nout,mom,gauge);
#endif

  ATOOLS::Vec4D *cpmom = new ATOOLS::Vec4D[m_nvector];
  for (size_t i=0;i<m_nvector;++i) cpmom[i]=mom[i];
  p_BS->CalcEtaMu(cpmom);
  delete [] cpmom;
  p_shand->Calculate();

  int acthel = m_helnumbers[hel];

  double M2  = p_ampl->Differential(acthel) * 
               p_hel->PolarizationFactor(acthel) * 
               sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS))/m_helalphas[hel]; 

  AddToHelicity(M2,hel);
  return M2;
}


void Single_Process::AddToHelicity(const double M2,const int hel)
{
  m_throws++;
  m_helresult        += M2;
  m_helresult2       += M2*M2;

  m_helthrows[hel]++;
  m_helresults[hel]  += M2;
}

