//#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "Single_Process.H"
#include "Run_Parameter.H"
#include "ISR_Base.H"

#include "Running_AlphaS.H"
#include "Combined_Selector.H"
#include "prof.hh"

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
  p_hel(0), p_BS(0), p_ampl(0), p_shand(0), p_partner(this)
{
  m_newlib=false;
  m_save_max=0.;
  GenerateNames(m_nin,p_flin,p_plin,m_nout,p_flout,p_plout,m_name,m_ptypename,m_libname);
  m_pslibname = m_libname;

  PolarizationNorm();
  InitCuts();
  if (_seldata) p_sel = new Combined_Selector(m_nin,m_nout,p_fl,_seldata);
  else {
    if (m_nout>2)
      msg.Error()<<"Potential Error in Single_Process "<<m_name<<endl
		 <<"   No selection cuts specified. Init No_Selector !"<<endl;
    p_sel = new No_Selector();
  }

  double sum_massin = 0.,sum_massout = 0.;
  for (int i=0;i<m_nin;i++)  sum_massin  += p_flin[i].Mass();
  for (int i=0;i<m_nout;i++) sum_massout += p_flout[i].Mass();
  m_isrthreshold = ATOOLS::Max(sum_massin,sum_massout);

  p_ps   = new Phase_Space_Handler(this,p_isr,p_beam);
  
  // making directory
  if (m_gen_str>1) {
    unsigned int  mode_dir = 0755;
    mkdir((string("Process/")+m_ptypename).c_str(),mode_dir); 
  }
  msg.Tracking()<<"Initialized Single_Process : "<<m_name<<", "<<m_nvec<<", 1/norm = "<<1./m_Norm<<endl;
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
	_testmoms = new Vec4D[m_nvec];
	p_ps->TestPoint(_testmoms);
	rpa.gen.SetEcms(ecms);    
	Vec4D * dummys = new Vec4D[m_nvec];
	p_ps->TestPoint(dummys);
	delete [] dummys;
      }
      else {
	_testmoms = new Vec4D[m_nvec];
	p_ps->TestPoint(_testmoms);
      }
    }
    else {
      _testmoms = new Vec4D[m_nvec];
      p_ps->TestPoint(_testmoms);
    }
  }
  if (p_moms) { delete [] p_moms; }
  p_moms = new Vec4D[m_nvec]; 
  for (int i=0;i<m_nin+m_nout;i++) p_moms[i] = _testmoms[i];

  p_hel    = new Helicity(m_nin,m_nout,p_fl,p_pl);
  p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nvec,p_fl,p_b);  
  p_shand  = new String_Handler(m_gen_str,p_BS);

  p_ampl   = new Amplitude_Handler(m_nin+m_nout,p_fl,p_b,&m_pol,model,top,m_orderQCD,m_orderEW,
				   p_BS,p_shand,m_ptypename+string("/")+m_name);
  if (p_ampl->GetGraphNumber()==0) {
    msg.Tracking()<<"Single_Process::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return -1;
  }
  procs++;
  m_pol.Add_Extern_Polarisations(p_BS,p_fl,p_hel);
  p_BS->Initialize();

  switch (Tests()) {
  case 2 : 
    for (int j=0;j<links.size();j++) {
      if (ATOOLS::IsZero((links[j]->Result()-Result())/(links[j]->Result()+Result()))) {
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
      if (ATOOLS::IsZero((links[j]->Result()-Result())/(links[j]->Result()+Result()))) {
	msg.Tracking()<<"Test : 1.  Can map "<<m_name<<" on "<<links[j]->Name()<<endl;
	p_partner = links[j];
	m_pslibname = links[j]->PSLibName();
	break;
      } 
    }
    if (p_partner==this) {
      links.push_back(this);
      totalsize++;
    }
    if (CheckLibraries()) return 1;
    for (int j=0;j<links.size();j++) {
      if (ATOOLS::IsZero((links[j]->Result()-Result())/(links[j]->Result()+Result()))) {
	if (links[j]->NewLibs()) {
	  if (CheckStrings(links[j])) return 1;
	}
      }
    }
    if (p_partner!=this) {
      links.push_back(this);
      totalsize++;
    }
    if (m_gen_str<2) return 1;
    WriteLibrary();
    if (p_partner==this) SetUpIntegrator();
    return 0;
  default :
    msg.Error()<<"Error in Single_Process::InitAmplitude : Failed for "<<m_name<<"."<<endl;
    errs.push_back(this);
    return 1;
  }
}


int Single_Process::InitAmplitude(Interaction_Model_Base * model,Topology * top)
{
  if (p_moms) { delete [] p_moms; }
  p_moms   = new Vec4D[m_nvec]; 
  p_ps->TestPoint(p_moms);

  p_hel    = new Helicity(m_nin,m_nout,p_fl,p_pl);
  p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nvec,p_fl,p_b);  
  p_shand  = new String_Handler(m_gen_str,p_BS);
  p_ampl   = new Amplitude_Handler(m_nin+m_nout,p_fl,p_b,&m_pol,model,top,m_orderQCD,m_orderEW,
				   p_BS,p_shand,m_ptypename+string("/")+m_name);
  if (p_ampl->GetGraphNumber()==0) {
    msg.Tracking()<<"Single_Process::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return -1;
  }
  m_pol.Add_Extern_Polarisations(p_BS,p_fl,p_hel);
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
    m_pol.Set_Gauge_Vectors(m_nin+m_nout,p_moms,Vec4D(sqrt(3.),1.,1.,-1.));
    p_BS->Setk0(0);
    p_BS->CalcEtaMu(p_moms);  
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
    M2     *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
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
    if (M_doub[i]/M2g<1.e-20) {
      p_hel->SwitchOff(i);
      msg.Debugging()<<"Switch off zero helicity "<<i<<" : "
		     <<M_doub[i]<<"/"<<M_doub[i]/M2g<<endl;
    }
  }

  M2g    *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
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
  
      M2 = operator()(p_moms);
      gauge_test = string_test = 0;
    }
    if (!ATOOLS::IsEqual(M2,M2g)) {
      msg.Out()<<"Mapping file(1) : "<<abs(M2)<<endl
	       <<"Original    (2) : "<<abs(M2g)<<endl
	       <<"Cross check (T) : "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
      msg.Out()<<"WARNING: Library cross check not satisfied: "
	       <<M2<<" vs. "<<M2g<<"  difference:"<<abs(M2/M2g-1.)*100.<<"%"<<endl;
      return 0;
      msg.Out()<<"         assuming numerical reasons, continuing "<<endl;
    } 
    else {
      msg.Debugging()<<"Mapping file(1) : "<<abs(M2)<<endl
		     <<"Original    (2) : "<<abs(M2g)<<endl;
      if (M2g!=0.)
	msg.Debugging()<<"Cross check (T) : "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
      else
	msg.Debugging()<<"Cross check (T) : "<<0.<<"%"<<endl;
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
	  M2S      += p_ampl->Differential(i)*p_hel->PolarizationFactor(i)*p_hel->Multiplicity(i);
	}
      }
      msg.Debugging()<<endl;
      M2S *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
      if (!ATOOLS::IsEqual(M2g,M2S)) {
	msg.Out()<<"WARNING: String test not satisfied: "
		 <<M2g<<" vs. "<<M2S<<"  difference:"<<abs(M2g/M2S-1.)*100.<<"%"<<endl;
	if (!ATOOLS::IsZero(M2g)) return 0;
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
  
  int number  = 0;
  string proc = string("Process/")+m_ptypename+string("/V");
  string testname;
  double M2s, helvalue;

  for (;;) {
    ++number;
    sprintf(help,"%i",number);
    testname  = CreateLibName()+string("_")+string(help);
    if (shand1->SearchValues(m_gen_str,testname,p_BS)) {

      shand1->Calculate();
      
      M2s = 0.;
      ATOOLS::msg.Debugging()<<"Check "<<number<<" :";ATOOLS::msg.Debugging().flush();
      for (short int i=0;i<p_hel->MaxHel();i++) {
	helvalue = p_ampl->Differential(shand1,i) * p_hel->PolarizationFactor(i) *
	  p_hel->Multiplicity(i);
	M2s     += helvalue;
	msg.Debugging()<<"*";ATOOLS::msg.Debugging().flush();
      }
      msg.Debugging()<<endl;
      M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
      msg.Debugging()<<"Cross check (1): "<<abs(M2s/Result()-1.)*100.<<"%"<<"  : "
		     <<M2s<<"/"<<Result()<<endl;
      if (ATOOLS::IsZero(abs((M2s-Result())/(M2s+Result())))) {
	msg.Tracking()<<"Found a suitable Library."<<endl;
	m_libname = testname;
	if (shand1) { delete shand1; shand1 = 0; }
	//Clean p_shand!!!!
	Minimize();
	CreateMappingFile();
	return 1;
      }
    } 
    else break;
  }
  if (shand1) { delete shand1; shand1 = 0; }
  return 0;
}

int Single_Process::CheckStrings(Single_Process* tproc)
{
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
  M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
  msg.Debugging()<<"Cross check (2): "<<abs(M2s/Result()-1.)*100.<<"%"<<"  : "
		 <<M2s<<"/"<<Result()<<endl;
  (shand1->Get_Generator())->ReStore();
  delete shand1;

  if (ATOOLS::IsZero(abs((M2s-Result())/(M2s+Result())))) {
    msg.Tracking()<<"Found a suitable string."<<endl;
    m_libname = tproc->LibName();
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
  int number  = 0;
  string testname;
  for (;;) {
    sprintf(help,"%i",number);
    testname    = CreateLibName()+string("_")+string(help);
    if (!(IsFile(string("Process/")+m_ptypename+string("/")+testname+string("/V.H")))) break;
    ++number;
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


void Single_Process::InitAnalysis(std::vector<ATOOLS::Primitive_Observable_Base *> _obs) {
  p_analysis = new ATOOLS::Primitive_Analysis(this->Name());
  for (int i=0;i<_obs.size();i++) {
    p_analysis->AddObservable(_obs[i]->GetCopy());
  }
  m_analyse  = 1;
}

/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/


bool Single_Process::SetUpIntegrator() {  
  p_sel->BuildCuts(p_cuts);
  if (m_nin==2) {
    if ( (p_fl[0].Mass() != p_isr->Flav(0).Mass()) ||
	 (p_fl[1].Mass() != p_isr->Flav(1).Mass()) ) p_isr->SetPartonMasses(p_fl);
    if (CreateChannelLibrary()) {
      if (p_ps->CreateIntegrators()) return 1;
    }
  }
  if (m_nin==1) return p_ps->CreateIntegrators();
  return 0;
}

bool Single_Process::CreateChannelLibrary()
{
  p_psgen     = new Phase_Space_Generator(m_nin,m_nout);
  bool newch  = 0;
  if (m_nin>1)  newch = p_psgen->Construct(p_ps->FSRIntegrator(),m_ptypename,m_pslibname,p_fl,this); 

  if (newch) {
    msg.Error()<<p_ps->NumberOfFSRIntegrators()<<" new Channels produced for "<<m_pslibname<<" ! "<<endl
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

  if (p_moms)     { delete [] p_moms;  p_moms     = 0; }
  if (p_sel)      { delete p_sel;      p_sel      = 0; }
  if (p_cuts)     { delete p_cuts;     p_cuts     = 0; }
  if (p_ps)       { delete p_ps;       p_ps       = 0; }
  if (p_analysis) { delete p_analysis; p_analysis = 0; }
}

void Single_Process::Empty() {
  if (p_cuts)        { delete p_cuts; p_cuts = 0; }
  if (p_ps)          { delete p_ps; p_ps = 0; } 
  if (p_partner != this) {
    return;
  }
}

void Single_Process::SetTotalXS(int _tables)  { 
  if (_tables!=2) {
    if (m_analyse) p_analysis->FinishAnalysis(m_resdir+string("/Tab")+m_name,_tables);
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
    msg.Events()<<"      xs for "<<m_name<<" : "
		<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
		<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
		<<"       max : "<<m_max<<endl
		<<"   exp. eff: "<<(100.*m_totalxs/m_max)<<"%"<<endl;
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
		  <<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
		  <<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
		  <<"       max : "<<m_max<<endl;
      from.close();
      p_ps->ReadIn(m_resdir+string("/MC_")+m_name);
      if (p_ps->BeamIntegrator() != 0) p_ps->BeamIntegrator()->Print();
      if (p_ps->ISRIntegrator()  != 0) p_ps->ISRIntegrator()->Print();
      if (p_ps->FSRIntegrator()  != 0) p_ps->FSRIntegrator()->Print();

      if (m_totalxs>=0.) return 1;
      return 0;
    }
  }
  m_totalxs = p_ps->Integrate();
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
      p_ps->WriteOut(_resdir+string("/MC_")+m_name);
      to.close();
    }
    return 1;
  }
  return 0;      
}

void Single_Process::WriteOutXSecs(std::ofstream & _to)    
{ 
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
  if (m_nin==2) m_totalxs /= ATOOLS::rpa.Picobarn();

  if (!(ATOOLS::IsZero((m_n*m_totalxs-m_totalsum)/(m_n*m_totalxs+m_totalsum)))) {
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
  int iter=1;
  if (m_n<200000)        iter =  40000;
  else if  (m_n<400000)  iter = 200000;
  // ==== uncomment for better unweighting efficiency ===
  //   if (iter!=1 && m_n%iter==0) {
  //     m_max = m_save_max;
  //     m_save_max = 0.;
  //   }
  if (value>m_save_max) m_save_max = value;
  if (m_analyse) p_analysis->DoAnalysis(value*rpa.Picobarn());
}

double Single_Process::Differential(ATOOLS::Vec4D* _moms) { return DSigma(_moms,0); }

double Single_Process::Differential2() { 
  if (p_isr->On()==0) return 0.;
  return DSigma2(); 
}


double Single_Process::DSigma(ATOOLS::Vec4D* _moms,bool lookup)
{
  m_last = m_lastdxs = 0.;
  for (int i=0;i<m_nin+m_nout;i++) {
    if (_moms[i][0] < p_fl[i].PSMass()) return m_last = 0.;
  }
  if (p_partner == this) {
    m_lastdxs = operator()(_moms);
  }
  else {
    if (lookup) m_lastdxs = p_partner->LastXS();
           else m_lastdxs = p_partner->operator()(_moms);
  }
  if (m_lastdxs <= 0.)                  return m_lastdxs = m_last = 0.;
  if (m_nin==2) {
    m_lastlumi = p_isr->Weight(p_flin);
    int    pols[2] = {p_pl[0].type[0],p_pl[1].type[0]};
    double dofs[2] = {p_pl[0].factor[0],p_pl[1].factor[0]};
    if (p_pl[0].num>1) pols[0] = 99;
    if (p_pl[1].num>1) pols[1] = 99;
    m_lastlumi *= p_beam->Weight(pols,dofs);
  }
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

double Single_Process::operator()(ATOOLS::Vec4D * mom)
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

  return M2 * sqr(m_pol.Massless_Norm(m_nin+m_nout,p_fl,p_BS));
}


bool   Single_Process::OneEvent(double _mass) { return (p_ps->OneEvent(_mass)); }
bool   Single_Process::SameEvent()            { return (p_ps->SameEvent()); }
double Single_Process::WeightedEvent()        { return (p_ps->WeightedEvent()); }
double Single_Process::SameWeightedEvent()    { return (p_ps->SameWeightedEvent()); }




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
  if (!(ATOOLS::rpa.gen.Tracking())) return;
  ATOOLS::msg.Out()<<m_name<<" : "<<m_last<<" -> "
		      <<m_lastdxs<<" @ "<<m_lastlumi<<", "<<endl;
}
