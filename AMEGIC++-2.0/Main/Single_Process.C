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

#include "XS_Selector.H"
#include "Single_XS.H"

using namespace AMEGIC;
using namespace PHASIC;
using namespace EXTRAXS;
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

Single_Process::Single_Process(int _nin,int _nout,Flavour* _fl,
			       ISR::ISR_Handler * _isr,BEAM::Beam_Handler * _beam,
			       APHYTOOLS::Selector_Data * _seldata,
			       int _gen_str,int _kfactorscheme, int _scalescheme, 
			       Pol_Info * _pl, int _runmode) 
{
  nin = _nin; nout = _nout; isr = _isr; beam = _beam; gen_str = _gen_str;
  kfactorscheme = _kfactorscheme;
  scalescheme   = _scalescheme;
  runmode       = _runmode; 

  flin  = new Flavour[nin];
  flout = new Flavour[nout];  
  plin  = new Pol_Info[nin];
  plout = new Pol_Info[nout]; 
  for (short int i=0;i<nin;i++) {
    flin[i]  = _fl[i];
    if (_pl!=0) plin[i] = _pl[i];
           else plin[i] = Pol_Info(flin[i]); 

  }
  for (short int i=0;i<nout;i++) { 
    flout[i] = _fl[i+nin]; 
    if (_pl!=0) plout[i] = _pl[i+nin];
           else plout[i] = Pol_Info(flout[i]); 
  }

  nstrong = 0;
  neweak  = 0;
  for (int i=0;i<nin;i++) {
    if (flin[i].Strong())      nstrong++;
    if (flin[i].IntCharge()!=0)  neweak++;
  }
  for (int i=0;i<nout;i++) {
    if (flout[i].Strong())     nstrong++;
    if (flout[i].IntCharge()!=0) neweak++;
  }

 xsflag = 0;
  if (_runmode!=AMPLITUDE_MODE) {
    if ((xsflag = FindXS()) != 0) {
      msg.Debugging()<<"In Single_Process : xs is available as fast function "<<std::endl;
    }
    if ((_runmode==XS_MODE) && (xsflag == 0)) {
      msg.Debugging()<<"Error in Single_Process : xs not available as fast function ! "<<std::endl
		     <<"                          Delete this process ! "<<std::endl;
      return;
    }
  }

  GenerateNames(nin,flin,plin,nout,flout,plout,name,ptypename,libname);

  PolarizationNorm();
  Initialize(_seldata);

  // making directory
  int  mode_dir = 448;

  mkdir((string("Process/")+ptypename).c_str(),mode_dir); 

  msg.Tracking()<<"Initialized Single_Process : "<<name
		<<", "<<nvec<<", 1/norm = "<<1./Norm<<endl;

}

Single_Process::~Single_Process()
{
  if (fl)       { delete [] fl;    fl       = 0; }
  if (flin)     { delete [] flin;  flin     = 0; }
  if (flout)    { delete [] flout; flout    = 0; }
  if (pl)       { delete [] pl;    pl       = 0; }
  if (plin)     { delete [] plin;  plin     = 0; }
  if (plout)    { delete [] plout; plout    = 0; }
  if (b)        { delete [] b;     b        = 0; }
  if (moms)     { delete [] moms;  moms     = 0; }
 
  if (sel)      { delete sel;      sel      = 0; }
  if (cuts)     { delete cuts;     cuts     = 0; }
  if (ps)       { delete ps;       ps       = 0; }
  if (analysis) { delete analysis; analysis = 0; }

  if (hel)      delete hel;
  if (BS)       delete BS;
  if (shand)    delete shand;
  if (ampl)     delete ampl;
}

/*------------------------------------------------------------------------------
  
  Generic stuff for initialization of Single_Processes
      
  ------------------------------------------------------------------------------*/

void Single_Process::Initialize(APHYTOOLS::Selector_Data * _seldata) {
  cuts = 0;
  InitCuts();
  ps   = new Phase_Space_Handler(this,isr,beam);

  if (_seldata) sel = new Combined_Selector(nin,nout,fl,_seldata);
  else {
    msg.Error()<<"Potential Error in Single_Process "<<name<<endl
	       <<"   No selection cuts specified. Init No_Selector !"<<endl;
    sel = new No_Selector();
  }

  moms = 0;

  totalxs  = totalerr = totalsum = totalsumsqr = 0.;
  last     = lastdxs  = max      = 0.;
  lastlumi = 1.;

  partner  = this;
  selected = this;
  atoms    = 1;
  analyse  = 0;
}

void Single_Process::PolarizationNorm() {
  int                   is_massless_pol  = pol.Massless_Vectors(nin,flin);
  if (!is_massless_pol) is_massless_pol  = pol.Massless_Vectors(nout,flout);
  int                   nmassive_pols    = pol.Massive_Vectors(nin,flin);
  nmassive_pols                         += pol.Massive_Vectors(nout,flout);

  // arrange Flavours
  nvec = nin+nout+is_massless_pol+nmassive_pols;
  fl   = new Flavour[nvec];
  pl   = new Pol_Info[nvec];
  b    = new int[nvec];
  for (short int i=0;i<nin;i++)         { fl[i] = flin[i]         ; pl[i] = plin[i]     ; b[i] = -1; }
  for (short int i=nin;i<nin+nout;i++)  { fl[i] = flout[i-nin]    ; pl[i] = plout[i-nin]; b[i] = 1; } 
  for (short int i=nin+nout;i<nvec;i++) { fl[i] = Flavour(kf::pol); b[i] = 1; }

  ModNorm = pol.Spin_Average(nin,flin);
  Norm = SymmetryFactors() * ModNorm;
#ifndef Explicit_Pols
  pol.Attach(nin+nout,fl);
  Norm *= pol.Massive_Norm();
  ModNorm *= pol.Massive_Norm();
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
    for (int j=0;j<nout;j++) {
      if (flout[j]==hflav)                                      cp++;
      else {
	if ((flout[j]==hflav.Bar()) && (hflav != hflav.Bar()))  cap++;
      }
    }
    if (cp>1)  sym *= double(fak(cp));
    if (cap>1) sym *= double(fak(cap));
  } 
  return 1./sym;
}

void Single_Process::InitCuts() 
{
  if (cuts == 0) {
    cuts = new Cut_Data();
    cuts->Init(nin+nout,fl);
  }
}

/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/

int Single_Process::InitAmplitude(Topology* top,Vec4D *& _testmoms,
				   vector<double> & results,vector<Single_Process *> & links)
{
  if (xsflag && (runmode != AMPLITUDE_MODE)) {
    results.push_back(1.);
    links.push_back(this);
    return 1;
  }
  if (!xsflag && (runmode == XS_MODE)) return 0;
  if (_testmoms==0) {
    msg.Debugging()<<"Init moms : "<<_testmoms<<" : "<<nin+nout<<endl;
    _testmoms = new Vec4D[nvec];
    ps->TestPoint(_testmoms);
    for (int i=0;i<nin+nout;i++) msg.Debugging()<<i<<" th mom : "<<_testmoms[i]<<" ("<<_testmoms[i].Abs2()<<")"<<endl;
  }
  if (moms) { delete [] moms; }
  moms = new Vec4D[nvec]; 
  for (int i=0;i<nin+nout;i++) moms[i] = _testmoms[i];

  double result = 0.;
  hel    = new Helicity(nin,nout,fl,pl);
  BS     = new Basic_Sfuncs(nin+nout,nvec,fl,b);  
  shand  = new String_Handler(gen_str,BS);

  ampl   = new Amplitude_Handler(nin+nout,fl,b,&pol,top,BS,shand,ptypename+string("/")+name);
  if (ampl->GetGraphNumber()==0) {
    msg.Error()<<"Single_Process::InitAmplitude : No diagrams for "<<name<<"."<<endl;
    return -1;
  }
  pol.Add_Extern_Polarisations(BS,fl,hel);
  BS->Initialize();
  shand->Initialize(ampl->GetGraphNumber(),hel->Max_Hel());

  switch (Tests(result)) {
  case 2 : 
    for (int j=0;j<results.size();j++) {
      if (AMATOOLS::IsZero((results[j]-result)/(results[j]+result))) {
	msg.Tracking()<<"Test : 2.  Can map "<<name<<" on "<<links[j]->Name()<<endl;
	partner = links[j];
      }
    }
    if (partner==this) {
      results.push_back(result);
      links.push_back(this);
    }
    return 1;
  case 1 :
    for (int j=0;j<results.size();j++) {
      if (AMATOOLS::IsZero((results[j]-result)/(results[j]+result))) {
	msg.Tracking()<<"Test : 1.  Can map "<<name<<" on "<<links[j]->Name()<<endl;
	partner = links[j];
      }
    }
    if (partner==this) {
      results.push_back(result);
      links.push_back(this);
    }
    return InitLibrary(result);
  default :
    msg.Error()<<"Error in Single_Process::InitAmplitude : Failed for "<<name<<"."<<endl;
    return -2;
  }
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
  

  double M2 = 0.;
  if (gauge_test) {
    pol.Set_Gauge_Vectors(nin+nout,moms,Vec4D(sqrt(3.),1.,1.,-1.));
    BS->set_k0(0);
    BS->setS(moms);  
    BS->ResetS_GT(.9);

    ampl->SetStringOff();

    msg.Debugging()<<number<<" :";AORGTOOLS::msg.Debugging().flush();
    for (short int i=0;i<hel->Max_Hel();i++) { 
      if (hel->On(i)) {
	M2 +=  ampl->Differential(i,(*hel)[i])*hel->PolarizationFactor(i); 
	msg.Debugging()<<"*";msg.Debugging().flush();
      } 
      else {
	msg.Debugging()<<"0";msg.Debugging().flush();
      }
    }
    msg.Debugging()<<endl;
    M2     *= sqr(pol.Massless_Norm(nin+nout,fl,BS));
    result  = M2;
  }

  ampl->ClearCalcList();
  // To prepare for the string test.
  ampl->SetStringOn();
  (shand->Get_Generator())->Reset();
  ampl->FillCoupling(shand);

  /* ---------------------------------------------------
     
     First test : gauge test

     --------------------------------------------------- */
#ifndef Explicit_Pols 
  Vec4D gauge(sqrt(3.),1.,1.,1.);
  if (nout==4) gauge = moms[4];
  if (nout==5) gauge = moms[4];
  //?????????
  if (nout==6) gauge = moms[4];  
  pol.Reset_Gauge_Vectors(nin+nout,moms,gauge);
#else
  BS->set_k0(1);
#endif

  BS->setS(moms);
  number++;
  msg.Debugging()<<number<<" :";

  double M2g = 0.;
  double * M_doub = new double[hel->Max_Hel()];

  for (short int i=0;i<hel->Max_Hel();i++) { 
    if (hel->On(i)) {
      M_doub[i]  = ampl->Differential(i,(*hel)[i])*hel->PolarizationFactor(i);  
      M2g       += M_doub[i];
      msg.Debugging()<<"*";msg.Debugging().flush();
    }
  }
  msg.Debugging()<<endl;

  //shorten helicities
  for (short int i=0;i<hel->Max_Hel();i++) {
    if (M_doub[i]/M2g<1.e-30) {
      hel->switch_off(i);
      msg.Debugging()<<"Switch off zero helicity "<<i<<" : "
		     <<ampl->Differential(i,(*hel)[i])<<"/"<<M_doub[i]/M2g<<endl;
    }
    if (hel->On(i)) {
      for (short int j=i+1;j<hel->Max_Hel();j++) {
	if (hel->On(j)) {
	  if (AMATOOLS::IsEqual(M_doub[i]/M2g,M_doub[j]/M2g)) {
	    msg.Debugging()<<"Switch off helicity "<<j<<" -> "<<i<<" : "
			   <<ampl->Differential(i,(*hel)[i])<<"/"
			  <<M_doub[i]/M2g<<" -> "<<M_doub[j]/M2g<<endl;
	    hel->switch_off(j);
	    hel->SetPartner(i,j);
	    hel->Inc_Mult(i);
	  }
	}
      }
    }
  }
  M2g    *= sqr(pol.Massless_Norm(nin+nout,fl,BS));
  result  = M2g;

  delete[] M_doub;
  ampl->ClearCalcList();  
  

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
    if (shand->SearchValues(gen_str,testname,BS)) {
      pol.Set_Gauge_Vectors(nin+nout,moms,Vec4D(sqrt(3.),1.,1.,-1.));
      BS->setS(moms);  
      shand->Initialize(ampl->GetGraphNumber(),hel->Max_Hel());
      (shand->Get_Generator())->Reset();
      ampl->FillCoupling(shand);
      shand->Complete(hel);
      M2 = operator()(moms);
      gauge_test = string_test = 0;
    }
    msg.Debugging()<<"Mapping file(1) : "<<abs(M2)<<endl
		   <<"Original    (2) : "<<abs(M2g)<<endl
		   <<"Cross check     : "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
    if (!AMATOOLS::IsZero(abs(M2/M2g-1.))) {
      msg.Tracking()<<"Cross check not satisfied: "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
      return 0;
    }
    libname    = testname;
    return 2;
  }

  /* ---------------------------------------------------
     
     Second test : string test

     --------------------------------------------------- */

  {
    PROFILE_LOCAL("Shand.Complete()");
    shand->Complete(hel);
  }
  if (string_test) {
    //String-Test
    if (shand->Is_String()) {
      double  M2S = 0.;
      shand->Calculate();
      msg.Debugging()<<"3:";msg.Debugging().flush();
      for (short int i=0;i<hel->Max_Hel();i++) {
	if (hel->On(i)) {
	  msg.Debugging()<<"*";msg.Debugging().flush();
	  M2S      += ampl->Differential(i)*hel->PolarizationFactor(i)*hel->Multiplicity(i);
	}
      }
      msg.Debugging()<<endl;
      M2S *= sqr(pol.Massless_Norm(nin+nout,fl,BS));
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
  if (shand->IsLibrary()) {
    msg.Tracking()<<"No library needs to be initialized. Already done."<<endl;
    return 1;
  }

  char help[20];
  sprintf(help,"%i",ampl->GetGraphNumber());
  libname += string("_");
  libname += string(help);
  sprintf(help,"%i",shand->NumberOfCouplings());
  libname += string("_");
  libname += string(help);
  sprintf(help,"%i",shand->NumberOfZfuncs());
  libname += string("_");
  libname += string(help);
  msg.Debugging()<<"In Init library for "<<name<<endl;


  String_Handler * shand1;
  shand1      = new String_Handler(shand->Get_Generator());
  
  int number  = 0;
  string proc = string("Process/")+ptypename+string("/V");
  string testname;
  double M2s;

  for (;;) {
    ++number;
    sprintf(help,"%i",number);
    testname  = libname+string("_")+string(help);
    if (shand1->SearchValues(gen_str,testname,BS)) {
      shand1->Initialize(ampl->GetGraphNumber(),hel->Max_Hel());
      (shand1->Get_Generator())->Reset();
      ampl->FillCoupling(shand1);
      shand1->Calculate();
      
      M2s = 0.;
      AORGTOOLS::msg.Debugging()<<"Check "<<number<<" :";AORGTOOLS::msg.Debugging().flush();
      for (short int i=0;i<hel->Max_Hel();i++) {
	M2s     += ampl->Differential(shand1,i) * hel->PolarizationFactor(i) *
	  hel->Multiplicity(i);
	msg.Debugging()<<"*";AORGTOOLS::msg.Debugging().flush();
      }
      msg.Debugging()<<endl;
      M2s *= sqr(pol.Massless_Norm(nin+nout,fl,BS));
      msg.Debugging()<<"Cross check: "<<abs(M2s/result-1.)*100.<<"%"<<"  : "
		     <<M2s<<"/"<<result<<endl;
      if (AMATOOLS::IsZero(abs((M2s-result)/(M2s+result)))) {
	msg.Tracking()<<"Found a suitable string."<<endl;
	libname = testname;
	if (shand->SearchValues(gen_str,testname,BS)) {
	  shand->Initialize(ampl->GetGraphNumber(),hel->Max_Hel());
	  (shand->Get_Generator())->Reset();
	  ampl->FillCoupling(shand);
	  
	  shand->Calculate();
	  M2s = 0.;
	  AORGTOOLS::msg.Debugging()<<number<<" :";AORGTOOLS::msg.Debugging().flush();
	  for (short int i=0;i<hel->Max_Hel();i++) {
	    M2s     += ampl->Differential(shand,i) * hel->PolarizationFactor(i) * 
	      hel->Multiplicity(i);
	    msg.Debugging()<<"*";AORGTOOLS::msg.Debugging().flush();
	  }
	  msg.Debugging()<<endl;
	  M2s *= sqr(pol.Massless_Norm(nin+nout,fl,BS));
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
    testname    = libname+string("_")+string(help);
    if (!(IsFile(string("Process/")+ptypename+string("/")+testname+string("/V.H")))) break;
    ++number;
  }
  libname = testname;
  if (partner==this) {
    msg.Debugging()<<"Write Library for "<<name<<" = "<<libname<<endl;
    int  mode_dir = 448;
    msg.Debugging()<<" ptypename = "<<ptypename<<endl<<" libname = "<<libname<<endl;
    mkdir((string("Process/")+ptypename+string("/")+libname).c_str(),mode_dir); 
    shand->Output(hel,ptypename+string("/")+libname);
    CreateMappingFile();
    if (shand1) { delete shand1; shand1 = 0; }
    return 0;
  }
  else {
    if (partner->shand->IsLibrary()) {
      if (partner->shand->SearchValues(gen_str,partner->libname,BS)) {
	partner->shand->Initialize(ampl->GetGraphNumber(),hel->Max_Hel());
	(partner->shand->Get_Generator())->Reset();
	ampl->FillCoupling(partner->shand);
	partner->shand->Calculate();
	
	M2s = 0.;
	AORGTOOLS::msg.Debugging()<<"Check "<<number<<" :";AORGTOOLS::msg.Debugging().flush();
	for (short int i=0;i<hel->Max_Hel();i++) {
	  M2s     += ampl->Differential(partner->shand,i) * hel->PolarizationFactor(i) *
	    hel->Multiplicity(i);
	  msg.Debugging()<<"*";AORGTOOLS::msg.Debugging().flush();
	}
	msg.Debugging()<<endl;
	M2s *= sqr(pol.Massless_Norm(nin+nout,fl,BS));
	msg.Debugging()<<"Cross check: "<<abs(M2s/result-1.)*100.<<"%"<<"  : "
		       <<M2s<<"/"<<result<<endl;
	if (AMATOOLS::IsZero(abs((M2s-result)/(M2s+result)))) {
	  libname = partner->libname;
	  CreateMappingFile();
	}
	else {
	  libname = testname;
	  msg.Debugging()<<"Write Library for "<<name<<" = "<<libname<<endl;
	  int  mode_dir = 448;
	  msg.Debugging()<<" ptypename = "<<ptypename<<endl<<" libname = "<<libname<<endl;
	  mkdir((string("Process/")+ptypename+string("/")+libname).c_str(),mode_dir); 
	  shand->Output(hel,ptypename+string("/")+libname);
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
  sprintf(outname,"%s.map",(string("Process/")+ptypename+string("/")+name).c_str());
  if (IsFile(outname)) {
    msg.Debugging()<<"Mapping file already exists. Check identity."<<endl;
    ifstream from;
    from.open(outname);
    string tempname;
    from>>tempname;
    if (tempname != libname) {
      msg.Error()<<"Files do not coincide. Maybe changed input data ?"<<endl;
      abort();
    }
    else return;
  }
  msg.Debugging()<<" Write Mapping Information to : "<<outname<<endl;

  ofstream to;
  to.open(outname,ios::out);
  to<<libname<<endl;
  to.close();
  msg.Debugging()<<"File "<<outname<<" saved."<<endl;  

}

bool Single_Process::FoundMappingFile(std::string & tempname) {
  char outname[100];
  sprintf(outname,"%s.map",(string("Process/")+ptypename+string("/")+name).c_str());
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

bool Single_Process::IsFile(string filename)
{
  msg.Debugging()<<"Check for "<<filename<<endl;
  ifstream from;
  bool     hit = 0;
  from.open(filename.c_str());
  if (from) hit = 1;
  from.close();
  return hit;
}

void Single_Process::InitAnalysis(std::vector<APHYTOOLS::Primitive_Observable_Base *> _obs) {
  analysis = new APHYTOOLS::Primitive_Analysis(this->Name());
  for (int i=0;i<_obs.size();i++) {
    analysis->AddObservable(_obs[i]->GetCopy());
  }
  analyse  = 1;
}

bool Single_Process::SetUpIntegrator() {  
  sel->BuildCuts(cuts);
  if (nin==2) {
    if ( (fl[0].Mass() != rpa.gen.Beam1().Mass()) ||
	 (fl[1].Mass() != rpa.gen.Beam2().Mass()) ) isr->SetPartonMasses(fl);
  }

  if (ps->CreateChannelLibrary(ptypename,libname)) {
    if (ps->CreateIntegrators()) return 1;
  }
  return 0;
}

/*------------------------------------------------------------------------------
  
  Process management
  
  ------------------------------------------------------------------------------*/

void Single_Process::Empty() {
  if (cuts)        { delete cuts; cuts = 0; }
  if (ps)          { delete ps; ps = 0; } 
  /*
    cout<<" sel="<<sel<<endl;
    cout<<" hel="<<hel<<endl;
    cout<<" BS ="<<BS<<endl;
    cout<<"shand="<<shand<<endl;
    cout<<"sgen="<<shand->Get_Generator()<<endl;
    cout<<"moms="<<moms<<endl;
  */
  if (partner != this) {
    /*
    if (sel)         { delete sel; sel = 0; } 
    if (hel)         { delete hel; hel = 0; }
    if (BS)          { delete BS; BS = 0; }
    if (shand)       { delete shand; shand = 0; }
    if (moms)        { delete [] moms; moms = 0; }
    */
    msg.Debugging()<<"Emptied "<<name<<" completely"<<endl;
    return;
  }
  msg.Debugging()<<"Emptied "<<name<<" partially"<<endl; 
}

void Single_Process::SetTotalXS(int tables)  { 
  if (analyse) analysis->FinishAnalysis(resdir+string("/Tab")+name,tables);
  totalxs  = totalsum/n; 
  totalerr = sqrt( (totalsumsqr/n - 
		    (AMATOOLS::sqr(totalsum)-totalsumsqr)/n/(n-1) )  / n); 
  AORGTOOLS::msg.Events()<<"      xs for "<<name<<" : "
			 <<totalxs*AORGTOOLS::rpa.Picobarn()<<" pb"
			 <<" +/- "<<totalerr/totalxs*100.<<"%,"<<endl
			 <<"       max : "<<max<<endl;
}

/*------------------------------------------------------------------------------

  Calculating total cross sections
  
  ------------------------------------------------------------------------------*/

bool Single_Process::CalculateTotalXSec() { 
  msg.Events()<<"In Single_Process::CalculateTotalXSec() for "<<name<<endl; 

  totalxs = ps->Integrate()/AORGTOOLS::rpa.Picobarn();
  if (!(AMATOOLS::IsZero((n*totalxs-totalsum)/(n*totalxs+totalsum)))) {
    msg.Error()<<"Result of PS-Integrator and internal summation to not coincide!"<<endl
	       <<"  "<<name<<" : "<<totalxs<<" vs. "<<totalsum/n<<endl;
  }
  SetTotalXS(0);
  if (totalxs>0.) return 1;
  return 0;
}

bool Single_Process::LookUpXSec(double ycut,bool calc,string obs) { 
  string filename = (resdir+string("/Tab")+name+string("/")+obs).c_str();
  if (IsFile(filename)) {
    Histogram * histo = new Histogram(filename);
    double * res      = new double[histo->Depth()];
    histo->Extrapolate(ycut,res,1);
    totalxs = res[0];
    max     = res[1];
    msg.Events()<<name<<" : Set total xsec and max at ycut = "<<ycut
		<<" : "<<endl<<"   "<<totalxs<<" / "<<max<<endl;
    delete histo;
    delete res;

    if (calc) {
      ps->ReadIn(resdir+string("/MC_")+name);
      if (ps->BeamIntegrator() != 0) ps->BeamIntegrator()->Print();
      if (ps->ISRIntegrator()  != 0) ps->ISRIntegrator()->Print();
      if (ps->FSRIntegrator()  != 0) ps->FSRIntegrator()->Print();
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
  msg.Events()<<"In Single_Process::PrepareXSecTables() for "<<name<<endl; 

  string filename = (resdir+string("/Tab")+name+string("dY_cut")).c_str();
  if (IsFile(filename)) {
    msg.Events()<<"Found "<<filename<<endl;
  }

  totalxs = ps->Integrate()/AORGTOOLS::rpa.Picobarn();
  if (!(AMATOOLS::IsZero((n*totalxs-totalsum)/(n*totalxs+totalsum)))) {
    msg.Error()<<"Result of PS-Integrator and internal summation to not coincide!"<<endl;
    msg.Error()<<"  "<<name<<" : "<<totalxs<<" vs. "<<totalsum/n<<endl;
  }
  SetTotalXS(1);
  ps->WriteOut(resdir+string("/MC_")+name);
  if (totalxs>0.) return 1;
  return 0;
}



void Single_Process::AddPoint(const double value) {
  if (rpa.gen.Debugging()) {
    if (n<=10 || (n>=500000 && n<=500010)) {
      msg.Debugging()<<"In Process_Base::AddPoint("<<value<<")"<<endl;  
      msg.Out()<<" n        = "<<n<<endl;
      msg.Out()<<" totalsum = "<<totalsum<<endl;
    }
  }
  n++;
  totalsum    += value;
  totalsumsqr += value*value;
  if (value>max) max = value;
  if (analyse) analysis->DoAnalysis(value*rpa.Picobarn());
}

double Single_Process::Differential(AMATOOLS::Vec4D* _moms) { return DSigma(_moms,0); }

double Single_Process::Differential2() { 
  if (isr->On()==0) return 0.;
  return DSigma2(); 
}

double Single_Process::DSigma(AMATOOLS::Vec4D* _moms,bool lookup)
{
  last = lastdxs = 0.;
  if (partner == this) {
    lastdxs = operator()(_moms);
    //    msg.Out()<<"Check this : "<<name<<" : "<<lastdxs * Norm<<" @ "<<isr->Weight(flin)<<endl;
  }
  else {
    if (lookup)        lastdxs = partner->LastXS();
           else        lastdxs = partner->operator()(_moms);
  }
  for (int i=0;i<nin+nout;i++) {
    if (_moms[i][0] < fl[i].PSMass()) return last = 0.;
  }
  if (lastdxs <= 0.)                  return lastdxs = last = 0.;
  if (nin==2)          lastlumi = isr->Weight(flin);
         else          lastlumi = 1.;

  return last = Norm * lastdxs * lastlumi;
}

double Single_Process::DSigma2() { 
  if ((flin[0]==flin[1]) || (isr->On()==0) ) return 0.;
  if (partner == this) {
    //    msg.Out()<<"Check this : "<<name<<" : "<<lastdxs * Norm<<" @ "<<isr->Weight2(flin)<<endl;
  }
  double tmp = Norm * lastdxs * isr->Weight2(flin); 
  last      += tmp;
  return tmp;
}

double Single_Process::operator()(AMATOOLS::Vec4D * mom)
{
  if (xsflag) return xsec->operator()(mom) / ModNorm;
  double M2 = 0.;

#ifndef Explicit_Pols   
  Vec4D gauge(sqrt(3.),1.,1.,1.);
  
  if (nout==4) gauge = mom[4];
  if (nout==5) gauge = mom[4];
  if (nout==6) gauge = mom[4];
  
  pol.Set_Gauge_Vectors(nin+nout,mom,gauge);
#endif
  BS->setS(mom);
  
  if (shand->Is_String()) {
    shand->Calculate();
    for (short int i=0;i<hel->Max_Hel();i++) {
      if (hel->On(i)) 
	M2 += ampl->Differential(i) * hel->Multiplicity(i) * hel->PolarizationFactor(i);
    }
  }
  else {
    for (short int i=0;i<hel->Max_Hel();i++)
      if (hel->On(i)) 
	M2 += ampl->Differential(i,(*hel)[i]) * hel->PolarizationFactor(i);
    shand->Complete(hel);
    ampl->ClearCalcList();
  }

  return M2 * sqr(pol.Massless_Norm(nin+nout,fl,BS));
}


bool Single_Process::OneEvent()        { return (ps->OneEvent()); }
bool Single_Process::SameEvent()       { return (ps->SameEvent()); }
double Single_Process::WeightedEvent() { return (ps->WeightedEvent()); }




int Single_Process::NumberOfDiagrams() { 
  if (xsflag) return IS_XS_FLAG;
  if (partner==this) return ampl->GetGraphNumber(); 
  return             partner->NumberOfDiagrams();
}

Point * Single_Process::Diagram(int i) { return ampl->GetPointlist(i); } 



/*------------------------------------------------------------------------------
  
  Helpers
  
  ------------------------------------------------------------------------------*/

void Single_Process::PrintDifferential()
{
  if (!(AORGTOOLS::rpa.gen.Debugging())) return;
  AORGTOOLS::msg.Out()<<name<<" : "<<last<<" -> "
		      <<lastdxs<<" @ "<<lastlumi<<", "<<endl;
}

int Single_Process::FindXS()
{
  xsec = 0;
  XS_Selector * xsselector = new XS_Selector(); 
  if ((xsec = xsselector->GetXS(nin, nout, fl))) {
    return 1; 
  }
  return 0;
}

void Single_Process::SetXS(EXTRAXS::Single_XS * _xsec)
{
  if (xsec) delete xsec;
  xsec = _xsec;
}



/*
    for (int j=0;j<results.size();j++) {
      if (AMATOOLS::IsZero((results[j]-result)/(results[j]+result))) {
	msg.Tracking()<<"Test : 1. Can map "<<name<<" on "<<links[j]->Name()<<endl;
	partner = links[j];

	string testname = string("");
	if (!FoundMappingFile(testname)) {
	  ampl->ClearCalcList();
	  String_Handler * shand1;
	  shand1 = new String_Handler(partner->shand->Get_Generator());
	  //if (shand1->SearchValues(gen_str,partner->libname,BS)) {
	  shand1->Initialize(ampl->GetGraphNumber(),hel->Max_Hel());
	  (shand1->Get_Generator())->Reset();
	  ampl->FillCoupling(shand1);
	  shand1->Calculate();
	  
	  double M1 = 0.;
	  double M2 = 0.;
	  for (short int i=0;i<hel->Max_Hel();i++) {
	    M2        += ampl->Differential(i) * hel->PolarizationFactor(i) *
	                 hel->Multiplicity(i);
	    M1        += ampl->Differential(i) * partner->hel->PolarizationFactor(i) *
	                 partner->hel->Multiplicity(i);
	    msg.Tracking()<<"check : helicity "<<i
			  <<" : "<<ampl->Differential(i)
			  <<" : "<<ampl->Differential(shand1,i)
			  <<" :: "<<partner->hel->Multiplicity(i)
			  <<" * "<<partner->hel->PolarizationFactor(i)<<endl;
	  }
	  M2    *= sqr(pol.Massless_Norm(nin+nout,fl,BS));
	  M2    *= sqr(partner->pol.Massless_Norm(nin+nout,fl,BS));
	  delete shand1;
	  
	  if (AMATOOLS::IsZero(abs(result/M2-1.))) {
	    msg.Tracking()<<"Cross check satisfied: "<<abs(result/M2-1.)*100.<<"%"<<endl;
	    msg.Tracking()<<"Cross check satisfied: "<<abs(result/M1-1.)*100.<<"%"<<endl;
	    libname = partner->LibName();
	    CreateMappingFile();
	    return 1;
	  }
	  else {
	    msg.Tracking()<<"Cross check not satisfied: "<<abs(result/M2-1.)*100.<<"%"<<endl;
	  }
	  //}
	}
	return 0;
	//return InitLibrary(result);


	return 1;
      }
    }
    results.push_back(result);
    links.push_back(this);
    return 0;
    // return InitLibrary(result);
*/
	  //string tmp = string("");
	  //if (!FoundMappingFile(tmp)) CreateMappingFile();
