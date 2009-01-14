#include "Single_Process.H"
#include "Run_Parameter.H"
#include "ISR_Base.H"

#include "Running_AlphaS.H"
#include "Combined_Selector.H"

#include "Histogram.H"

#include "Random.H"

#include "Shell_Tools.H"
#include "MyStrStream.H"
#include "Data_Reader.H"

#include <unistd.h>
#include <sys/stat.h>

using namespace AMEGIC;
using namespace MODEL;
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
			       int _kfactorscheme, PHASIC::scl::scheme _scalescheme,double _scale,
			       Pol_Info * _pl,int _nex,Flavour * _ex_fl, std::string cuttag,
			       double error,std::string e_func) :
  Process_Base(NULL,_nin,_nout,_fl,_isr,_beam,_gen_str,_orderQCD,_orderEW,
	       _scalescheme,_kfactorscheme,_scale,_pl,_nex,_ex_fl,cuttag,error),
  m_sfactor(1.), p_hel(0), p_BS(0), p_ampl(0), p_shand(0), p_partner(this), 
  m_helsample(false), m_inithelsample(false), m_throws(0), m_helresult(0.), m_helresult2(0.)
{
  string newpath=rpa.gen.Variable("SHERPA_CPP_PATH");
  ATOOLS::MakeDir(newpath);
  struct stat fst;
  if (stat((newpath+"/makelibs").c_str(),&fst)==-1 || 
      (fst.st_mode&S_IFMT)!=S_IFREG) {
    CopyFile(rpa.gen.Variable("SHERPA_SHARE_PATH")+"/makelibs",
	     newpath+"/makelibs");
    CopyFile(rpa.gen.Variable("SHERPA_SHARE_PATH")+"/grid_makelibs",
	     newpath+"/grid_makelibs");
  }
  m_newlib   = false;
  m_libnumb  = 0;
  m_save_max = 0.;
  m_efunc=e_func;
  m_pslibname = m_libname = ToString(m_nin)+"_"+ToString(m_nout);
  if (m_gen_str>1) m_ptypename = "P"+m_libname;
  else m_ptypename = "N"+m_libname;

  PolarizationNorm();
  if (_seldata) p_selector = new Combined_Selector(m_nin,m_nout,p_flavours,_seldata,m_cuttag);
  else {
    if (m_nout>2)
      msg_Out()<<"WARNING in Single_Process "<<m_name<<endl
	       <<"   No selection cuts specified. Init No_Selector !"<<endl;
    p_selector = new No_Selector();
  }
  p_selector->SetProcessName(Name());
  SetScaleScheme(m_scalescheme);

  double sum_massin = 0.,sum_massout = 0.;
  for (size_t i=0;i<m_nin;i++)  sum_massin  += p_flin[i].Mass();
  for (size_t i=0;i<m_nout;i++) sum_massout += p_flout[i].Mass();
  m_threshold = ATOOLS::Max(sum_massin,sum_massout);

  p_pshandler = new Phase_Space_Handler(this,p_isrhandler,p_beamhandler,m_maxerror);
  SetPSHandler(p_pshandler);

  // making directory
  if (m_gen_str>1) {
    ATOOLS::MakeDir(rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename); 
  }
  msg_Tracking()<<"Initialized Single_Process : "<<m_name<<"."<<std::endl;
}


Single_Process::Single_Process(Process_Info* pinfo,int _nin,int _nout,Flavour * _fl,
			       ISR_Handler * _isr,Beam_Spectra_Handler * _beam,Selector_Data * _seldata,
			       int _gen_str,int _orderQCD, int _orderEW,
			       int _kfactorscheme, PHASIC::scl::scheme _scalescheme,double _scale,
			       Pol_Info * _pl,int _nex,Flavour * _ex_fl, std::string cuttag,double error,std::string e_func) :   
  Process_Base(pinfo,_nin,_nout,_fl,_isr,_beam,_gen_str,_orderQCD,_orderEW,
	       _scalescheme,_kfactorscheme,_scale,_pl,_nex,_ex_fl,cuttag,error),
  m_sfactor(1.), p_hel(0), p_BS(0), p_ampl(0), p_shand(0), p_partner(this), 
  m_helsample(false), m_inithelsample(false), m_throws(0), m_helresult(0.), m_helresult2(0.)
{
  string newpath=rpa.gen.Variable("SHERPA_CPP_PATH");
  ATOOLS::MakeDir(newpath);
  struct stat fst;
  if (stat((newpath+"/makelibs").c_str(),&fst)==-1 || 
      (fst.st_mode&S_IFMT)!=S_IFREG) {
    CopyFile(rpa.gen.Variable("SHERPA_SHARE_PATH")+"/makelibs",
	     newpath+"/makelibs");
    CopyFile(rpa.gen.Variable("SHERPA_SHARE_PATH")+"/grid_makelibs",
	     newpath+"/grid_makelibs");
  }
  m_newlib   = false;
  m_libnumb  = 0;
  m_save_max = 0.;
  m_efunc=e_func;
  m_pslibname = m_libname = ToString(m_nin)+"_"+ToString(m_nout);
  if (m_gen_str>1) m_ptypename = "P"+m_libname;
  else m_ptypename = "N"+m_libname;

  PolarizationNorm();
  if (_seldata) p_selector = new Combined_Selector(m_nin,m_nout,p_flavours,_seldata,cuttag);
  else {
    if (m_nout>2)
      msg_Out()<<"WARNING in Single_Process "<<m_name<<endl
	       <<"   No selection cuts specified. Init No_Selector !"<<endl;
    p_selector = new No_Selector();
  }
  p_selector->SetProcessName(Name());
  SetScaleScheme(m_scalescheme);

  double sum_massin = 0.,sum_massout = 0.;
  for (size_t i=0;i<m_nin;i++)  sum_massin  += p_flin[i].Mass();
  for (size_t i=0;i<m_nout;i++) sum_massout += p_flout[i].Mass();
  m_threshold = ATOOLS::Max(sum_massin,sum_massout);

  p_pshandler = new Phase_Space_Handler(this,p_isrhandler,p_beamhandler,m_maxerror);
  SetPSHandler(p_pshandler);

  // making directory
  if (m_gen_str>1) {
    ATOOLS::MakeDir(rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename); 
  }
  msg_Tracking()<<"Initialized Single_Process : "<<m_name<<"."<<std::endl;
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
  for (size_t i=0;i<m_nin;i++) { 
    p_flavours[i] = p_flin[i]; 
    p_pl[i]       = p_plin[i]; 
    p_b[i]        = -1; 
  }
  for (size_t i=m_nin;i<m_nin+m_nout;i++)  { 
    p_flavours[i] = p_flout[i-m_nin]; 
    p_pl[i]       = p_plout[i-m_nin]; 
    p_b[i]        = 1; 
  } 
  for (size_t i=m_nin+m_nout;i<m_nvector;i++) { 
    p_flavours[i] = Flavour(kf_none); 
    p_b[i]        = 1; 
  }

  m_Norm = SymmetryFactors() * m_pol.Spin_Average(m_nin,p_flin);
}

double Single_Process::SymmetryFactors()
{
  Flavour* flav  = new Flavour[m_nvector];
  double sym = 1.;
  int ndecays=p_pinfo->Ndecays();
  for (int i=0;i<=ndecays;i++) {
    int j=i;
    Process_Info* pi=p_pinfo->GetDecay(j);
    size_t nout=pi->GetStableFlavList(flav);
    sym*=SBSymmetryFactor(flav,nout);
  }
  delete[] flav;
  return sym;
}

double Single_Process::SBSymmetryFactor(Flavour* flout,size_t nout)
{
  double sym = 1.;
  for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
      kfit!=s_kftable.end();++kfit) {
    Flavour hflav(kfit->first);
    if (hflav.IsHadron()) continue; 
    int cp  = 0;
    int cap = 0;
    for (size_t j=0;j<nout;j++) {
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

void Single_Process::FixISRThreshold()
{
  double m_mass_in  = 0.;
  double m_mass_out = 0.;
  
  for (size_t i = 0;i<m_nin;i++)  m_mass_in  += p_flin[i].Mass(); 
  for (size_t i = 0;i<m_nout;i++) m_mass_out += p_flout[i].Mass(); 
  
  double isrth = ATOOLS::Max(m_mass_in,m_mass_out);
  
  SetISRThreshold(isrth);

}

/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/



int Single_Process::InitAmplitude(Model_Base * model,Topology* top,Vec4D *& _testmoms,
				  vector<Process_Base *> & links,vector<Process_Base *> & errs,
				  int & totalsize, int & procs, int & current_atom)
{
  if (CheckAlternatives(links,current_atom,Name())) return 1;

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
	if (p_pinfo->OSDecays()) TestPoint(dummys);
	else p_pshandler->TestPoint(dummys);
	delete [] dummys;
      }
      else {
	_testmoms = new Vec4D[m_nvector];
	if (p_pinfo->OSDecays()) TestPoint(_testmoms);
	else p_pshandler->TestPoint(_testmoms);
      }
    }
    else {
      _testmoms = new Vec4D[m_nvector];
      if (p_pinfo->OSDecays()) TestPoint(_testmoms);
      else p_pshandler->TestPoint(_testmoms);
    }
  }
  if (p_momenta) { delete [] p_momenta; }
  p_momenta = new Vec4D[m_nvector]; 
  for (size_t i=0;i<m_nin+m_nout;i++) p_momenta[i] = _testmoms[i];

  p_hel    = new Helicity(m_nin,m_nout,p_flavours,p_pl);

  bool directload = true;
  int libchk=0; 
  Data_Reader reader(" ",";","!","=");
  if (reader.ReadFromFile(libchk,"ME_LIBCHECK")) {
    if (libchk==1) {
      msg_Info()<<"Enforce full library check. This may take some time"<<std::endl;
      directload = false;
    }
  }  
  if (directload) directload = FoundMappingFile(m_libname,m_pslibname);
  if (directload) {
    string hstr=rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename+"/"+m_libname;
    string hstr2=rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename+"/"+m_name+".map";
    p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nvector,p_flavours,p_b,hstr,hstr2);  
  }
  else p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nvector,p_flavours,p_b);  
  p_shand  = new String_Handler(m_gen_str,p_BS,model->GetVertex()->GetCouplings());
  
  p_ampl   = new Amplitude_Handler(m_nin+m_nout,p_flavours,p_b,p_pinfo,model,top,m_orderQCD,m_orderEW,
				   p_BS,p_shand,m_print_graphs,!directload);
  if (p_ampl->GetGraphNumber()==0) {
    msg_Tracking()<<"Single_Process::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return -1;
  }
  procs++;
  map<string,Complex> cplmap;
  for (size_t j=current_atom;j<links.size();j++) {
    cplmap.clear();
    if (p_ampl->CompareAmplitudes(links[j]->GetAmplitudeHandler(),m_sfactor,cplmap)) {
      if (p_hel->Compare(links[j]->GetHelicity(),m_nin+m_nout)) {
	m_sfactor = sqr(m_sfactor);
	msg_Tracking()<<"Single_Process::InitAmplitude : Found compatible process for "<<m_name<<" : "<<links[j]->Name()<<endl;
	  
	if (!FoundMappingFile(m_libname,m_pslibname)) {
	  string mlname = rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename+"/"+links[j]->Name();
	  string mnname = rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename+"/"+Name();
	  if (IsFile(mlname+string(".map"))) { 
	    if (m_sfactor==1.) CopyFile(mlname+".map",mnname+".map");
	    else {
	      UpdateMappingFile(mlname,cplmap);
	      CreateMappingFile((Single_Process*)links[j]);
	    }
	    CopyFile(mlname+".col",mnname+".col");
	  }
	}

	p_partner = (Single_Process*)links[j];
	WriteAlternativeName(p_partner->Name());
	
	Minimize();
	return 1;
      }
    }
  }
  if (directload) {
    p_ampl->CompleteLibAmplitudes(m_nin+m_nout,m_ptypename+string("/")+m_name,
				  m_ptypename+string("/")+m_libname);    
    if (!p_shand->SearchValues(m_gen_str,m_libname,p_BS)) return 0;
    if (!TestLib()) return 0;
    for (size_t j=current_atom;j<links.size();j++) {
      if (ATOOLS::IsEqual(links[j]->Result()*m_sfactor,Result())) {
	if (CheckMapping(links[j])) {
	  msg_Tracking()<<"Single_Process::InitAmplitude : "<<std::endl
			<<"   Found an equivalent partner process for "<<m_name<<" : "<<links[j]->Name()<<std::endl
			<<"   Map processes."<<std::endl;
	  p_partner = (Single_Process*)links[j];
	  break;
	}
      } 
    }
    if (p_partner==this) {
      links.push_back(this);
      totalsize++;
    }
    Minimize();
    return 1;
  }

  p_ampl->CompleteAmplitudes(m_nin+m_nout,p_flavours,p_b,&m_pol,
			     top,p_BS,m_ptypename+string("/")+m_name);
  m_pol.Add_Extern_Polarisations(p_BS,p_flavours,p_hel);
  p_BS->Initialize();

  switch (Tests()) {
  case 2 : 
    for (size_t j=current_atom;j<links.size();j++) {
      if (ATOOLS::IsEqual(links[j]->Result(),Result())) {
	if (CheckMapping(links[j])) {
	  msg_Tracking()<<"Single_Process::InitAmplitude : "<<std::endl
			<<"   Found an equivalent partner process for "<<m_name<<" : "<<links[j]->Name()<<std::endl
			<<"   Map processes."<<std::endl;
	  p_partner = (Single_Process*)links[j];
	  break;
	}
      } 
    }
    if (p_partner==this) {
      links.push_back(this);
      totalsize++;
    }
    Minimize();
    WriteAlternativeName(p_partner->Name());
    return 1;
  case 1 :
    for (size_t j=current_atom;j<links.size();j++) {
      if (ATOOLS::IsEqual(links[j]->Result(),Result())) {
	msg_Tracking()<<"Single_Process::InitAmplitude : "<<std::endl
		      <<"   Found a partner for process "<<m_name<<" : "<<links[j]->Name()<<std::endl;
	p_partner   = (Single_Process*)links[j];
	m_pslibname = links[j]->PSLibName();
	WriteAlternativeName(p_partner->Name());
	break;
      } 
    }
    if (Result()==0.) return -3;
    if (p_partner==this) links.push_back(this);
    
    if (CheckLibraries()) return 1;
    for (size_t j=0;j<links.size();j++) {
      if (links[j]->NewLibs()) {
	if (CheckStrings((Single_Process*)links[j])) return 1;	
      }      
    }
    if (p_partner!=this) links.push_back(this);
    
    if (m_gen_str<2) return 1;
    totalsize++;
    if (p_partner!=this) {
      msg_Tracking()<<"Single_Process::InitAmplitude : "<<std::endl
		    <<"   Strings of process "<<m_name<<" and partner "
		    <<p_partner->Name()<<" did not fit."<<std::endl
		    <<"   Have to write new library."<<std::endl;
    }
    WriteLibrary();
    if (p_partner==this && Result()>0.) SetUpIntegrator();
    return 0;
  case -3: return -3;
  default :
    msg_Error()<<"ERROR in Single_Process::InitAmplitude : "<<std::endl
	       <<"   Failed for "<<m_name<<"."<<endl;
    errs.push_back(this);
    return 1;
  }
}


int Single_Process::InitAmplitude(Model_Base * model,Topology * top)
{
  if (p_momenta) { delete [] p_momenta; }
  p_momenta   = new Vec4D[m_nvector]; 
  if (p_pinfo->OSDecays()) TestPoint(p_momenta);
  else p_pshandler->TestPoint(p_momenta);

  p_hel    = new Helicity(m_nin,m_nout,p_flavours,p_pl);

  bool directload = true;
  int libchk=0; 
  Data_Reader reader(" ",";","!","=");
  if (reader.ReadFromFile(libchk,"ME_LIBCHECK")) {
    if (libchk==1) {
      msg_Info()<<"Enforce full library check. This may take some time"<<std::endl;
      directload = false;
    }
  }  
  if (directload) directload = FoundMappingFile(m_libname,m_pslibname);
  if (directload) {
    string hstr=rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename+"/"+m_libname;
    string hstr2=rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename+"/"+m_name+".map";
    p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nvector,p_flavours,p_b,hstr,hstr2);  
  }
  else p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nvector,p_flavours,p_b);  
  p_shand  = new String_Handler(m_gen_str,p_BS,model->GetVertex()->GetCouplings());
  p_ampl   = new Amplitude_Handler(m_nin+m_nout,p_flavours,p_b,p_pinfo,model,top,m_orderQCD,m_orderEW,
				   p_BS,p_shand,!directload);
  if (p_ampl->GetGraphNumber()==0) {
    msg_Tracking()<<"Single_Process::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return -1;
  }
  if (directload) {
    p_ampl->CompleteLibAmplitudes(m_nin+m_nout,m_ptypename+string("/")+m_name,
				  m_ptypename+string("/")+m_libname);    
    if (!p_shand->SearchValues(m_gen_str,m_libname,p_BS)) return 0;
    if (!TestLib()) return 0;
    return 1;
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
    msg_Error()<<"Error in Single_Process::InitAmplitude : Failed for "<<m_name<<"."<<endl;
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

    msg_Info()<<"Single_Process::Tests for "<<m_name<<std::endl
	      <<"   Prepare gauge test and init helicity amplitudes. This may take some time."
	      <<std::endl;
    for (size_t i=0;i<p_hel->MaxHel();i++) { 
      if (p_hel->On(i)) {
	helvalue = p_ampl->Differential(i,(*p_hel)[i])*p_hel->PolarizationFactor(i); 
	M2      +=  helvalue;
      } 
    }
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
  p_BS->Setk0(1);
  p_BS->CalcEtaMu(p_momenta);
  number++;

  if (!gauge_test) p_ampl->SetStringOff();  //second test without string production 

  double M2g = 0.;
  double * M_doub = new double[p_hel->MaxHel()];

  /* Calculate the squared amplitude of the polarisation states. If a certain external
     polarisation combination is found not to contribute for the point in phase space
     tested, it is assumed that is doesnt contribute at all and is switched off.      */
  for (size_t i=0;i<p_hel->MaxHel();i++) { 
    if (p_hel->On(i)) {
      M_doub[i]  = p_ampl->Differential(i,(*p_hel)[i])*p_hel->PolarizationFactor(i);  
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
  msg_Tracking()<<"Single_Process::Tests for "<<m_name<<std::endl
		<<"   Switched off or mapped "<<switchhit<<" helicities."<<std::endl;

  M2g    *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS));
  m_iresult  = M2g;

  p_ampl->ClearCalcList();  
  p_ampl->FillCoupling(p_shand);
  p_ampl->KillZList();  
  p_BS->StartPrecalc();

  if (gauge_test) {
    if (!ATOOLS::IsEqual(M2,M2g)) {
      msg_Out()<<"WARNING:  Gauge test not satisfied: "
	       <<M2<<" vs. "<<M2g<<" : "<<dabs(M2/M2g-1.)*100.<<"%"<<endl
	       <<"Gauge(1): "<<abs(M2)<<endl
	       <<"Gauge(2): "<<abs(M2g)<<endl;
    }
    /*
      else {
      msg_Debugging()<<"Gauge(1): "<<abs(M2)<<endl
      <<"Gauge(2): "<<abs(M2g)<<endl;
      if (M2g!=0.)
      msg_Debugging()<<"Gauge test: "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
      else
      msg_Debugging()<<"Gauge test: "<<0.<<"%"<<endl;
      }
    */
  }
  else {
    delete[] M_doub;
    number++;
    if (p_shand->SearchValues(m_gen_str,testname,p_BS)) {
      p_shand->Initialize(p_ampl->GetRealGraphNumber(),p_hel->MaxHel());
      (p_shand->Get_Generator())->Reset();

      // Get a cross section from the operator() method to compare with M2g later on.
      p_hel->ForceNoTransformation();
      M2 = operator()(p_momenta);
      p_hel->AllowTransformation();
      gauge_test = string_test = 0;
    }
    else {
      string searchfilename = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+m_ptypename+string("/")+testname+string("/V.H");
      if (IsFile(searchfilename)) {
      	msg_Error()<<"ERROR in Single_Process::Tests()"<<std::endl
			   <<"   No compiled & linked library found for process "<<testname<<std::endl
			   <<"   but files already written out !"<<std::endl
			   <<om::bold<<"   Interrupt run and execute \"makelibs\" in '"
			   <<rpa.gen.Variable("SHERPA_CPP_PATH")<<"'."
			   <<om::reset<<std::endl;
	CopyFile(rpa.gen.Variable("SHERPA_SHARE_PATH")+"/makelibs",
		 rpa.gen.Variable("SHERPA_CPP_PATH")+"/makelibs");
	THROW(normal_exit,"Failed to load library.");
      }
      else {
      	msg_Error()<<"ERROR in Single_Process::Tests()"<<std::endl
			   <<"   Mapping file exists, but no compiled & linked library found for process "
			   <<testname<<std::endl
			   <<"   and no files written out !"<<std::endl
			   <<om::bold<<"   Interrupt run, execute \"makeclean\" in Run-directory and re-start."
			   <<om::reset<<std::endl;
	THROW(critical_error,"Failed to load library.");
      }
    }
    if (!ATOOLS::IsEqual(M2,M2g)) {
      if (abs(M2/M2g-1.)>rpa.gen.Accu()) {
	msg_Out()<<"WARNING: Library cross check not satisfied: "
		 <<M2<<" vs. "<<M2g<<"  difference:"<<abs(M2/M2g-1.)*100.<<"%"<<endl
		 <<"   Mapping file(1) : "<<abs(M2)<<endl
		 <<"   Original    (2) : "<<abs(M2g)<<endl
		 <<"   Cross check (T) : "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
	return 0;
      }
      else {
	msg_Out()<<"WARNING: Library cross check not satisfied: "
		 <<M2<<" vs. "<<M2g<<"  difference:"<<abs(M2/M2g-1.)*100.<<"%"<<endl
		 <<"   assuming numerical reasons with small numbers, continuing "<<endl;
      }
    }
    else {
      if (M2g==0.) {
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
    if (!rpa.gen.SpinCorrelation()) {
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
    }
    delete[] M_doub;
    p_shand->Complete(p_hel);

    if (p_shand->Is_String()) {
      double  M2S = 0.;
      p_shand->Calculate();
      
      for (size_t i=0;i<p_hel->MaxHel();i++) {
	if (p_hel->On(i)) {
	  M2S += p_ampl->Differential(i)*p_hel->PolarizationFactor(i)*p_hel->Multiplicity(i);
	}
      }
      M2S *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS));
      if (!ATOOLS::IsEqual(M2g,M2S)) {
	msg_Out()<<"WARNING: String test not satisfied: "
		 <<M2g<<" vs. "<<M2S<<"  difference:"<<abs(M2g/M2S-1.)*100.<<"%"<<endl;
	if (abs(M2g/M2S-1.)>rpa.gen.Accu()) {
	  return 0;
	}
	msg_Out()<<"         assuming numerical reasons, continuing "<<endl;
      }
      return 1;
    }
    return 1;
  }
  delete[] M_doub;

  return 0;
}

int Single_Process::TestLib()
{
  double M2(0.);
  double * M_doub = new double[p_hel->MaxHel()];

  p_BS->CalcEtaMu((ATOOLS::Vec4D*)p_momenta);
  p_hel->InitializeSpinorTransformation(p_BS);
  p_shand->Calculate();
  
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    M2 += M_doub[i] = p_ampl->Differential(i)*p_hel->Multiplicity(i)*p_hel->PolarizationFactor(i);
  } 
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (M_doub[i]==0. || dabs(M_doub[i]/M2)<(ATOOLS::Accu()*1.e-2)) {
      p_hel->SwitchOff(i);
    }
  }
  if (!rpa.gen.SpinCorrelation()) {
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
  }
  delete[] M_doub;

  m_iresult = M2 * sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS));
  if (m_iresult>0.) return 1;
  return 0;
}

int Single_Process::CheckLibraries() {
  if (m_gen_str==0) return 1;
  if (p_shand->IsLibrary()) return 1;

  msg_Info()<<"Single_Process::CheckLibraries : Looking for a suitable library. This may take some time."<<std::endl;
  String_Handler * shand1;
  shand1      = new String_Handler(p_shand->Get_Generator());
  
  m_libnumb  = 0;
  string proc = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+m_ptypename+string("/V");
  string testname;
  double M2s, helvalue;

  for (;;) {
    testname  = CreateLibName()+string("_")+ToString(m_libnumb);
    if (shand1->SearchValues(m_gen_str,testname,p_BS)) {
      shand1->Calculate();
      
      M2s = 0.;
      for (size_t i=0;i<p_hel->MaxHel();i++) {
	helvalue = p_ampl->Differential(shand1,i) * p_hel->PolarizationFactor(i) *
	  p_hel->Multiplicity(i);
	M2s     += helvalue;
	} 
      M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS));
      if (ATOOLS::IsEqual(M2s,Result())) {
	m_libname = testname;
	m_pslibname = testname;
	if (shand1) { delete shand1; shand1 = 0; }
	//Clean p_shand!!!!
	CreateMappingFile(this);
	Minimize();
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
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    helvalue = p_ampl->Differential(shand1,i) * p_hel->PolarizationFactor(i) *
      p_hel->Multiplicity(i);
    M2s     += helvalue;
  }
  M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,p_flavours,p_BS));
  (shand1->Get_Generator())->ReStore();
  delete shand1;

  if (ATOOLS::IsEqual(M2s,Result())) {
    m_libname = tproc->LibName();
    m_pslibname = tproc->PSLibName();
    CreateMappingFile(this);
    Minimize();
    return 1;
  }
  return 0;
}
  
void Single_Process::WriteAlternativeName(string aname) 
{
  if (aname==m_name) return;
  std::string altname = rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename+"/"+m_name+".alt";
  if (IsFile(altname)) return;
  std::ofstream to;
  to.open(altname.c_str(),ios::out);
  to<<aname<<" "<<m_sfactor<<endl;
  to.close();
}

bool Single_Process::CheckAlternatives(vector<Process_Base *> & links,int current_atom,string procname)
{
  std::string altname = rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename+"/"+procname+".alt";
  if (IsFile(altname)) {
    double factor;
    string name; 
    ifstream from;
    from.open(altname.c_str(),ios::in);
    from>>name>>factor;
    from.close();
    m_sfactor *= factor;
    for (size_t j=current_atom;j<links.size();j++) {
      if (links[j]->Name()==name) {
	p_partner = (Single_Process*)links[j];
	m_orderQCD=p_partner->OrderStrong();
	m_orderEW=p_partner->OrderEWeak();
	msg_Tracking()<<"Found Alternative process: "<<m_name<<" "<<name<<endl;
	return true;
      }
    }
    if (CheckAlternatives(links,current_atom,name)) return true;
  }
  m_sfactor = 1.;
  return false;
}

void Single_Process::WriteLibrary() 
{
  if (m_gen_str<2) return;
  string testname;
  string newpath=rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/");
  for (;;) {
    testname    = CreateLibName()+string("_")+ToString(m_libnumb);
    if (!(IsFile(newpath+m_ptypename+string("/")+testname+string("/V.H")))) break;
    ++m_libnumb;
  }
  m_libname = testname;
  if (p_partner==this) m_pslibname = m_libname;
                  else m_pslibname = p_partner->PSLibName();
  ATOOLS::MakeDir(newpath+m_ptypename+"/"+m_libname); 
  p_shand->Output(p_hel,m_ptypename+string("/")+m_libname);
  CreateMappingFile(this);
  p_BS->Output(newpath+m_ptypename+string("/")+m_libname);
  p_ampl->StoreAmplitudeConfiguration(newpath+m_ptypename+string("/")+m_libname);
  m_newlib=true;
  msg_Info()<<"Single_Process::WriteLibrary : "<<std::endl
	    <<"   Library for "<<m_name<<" has been written, name is "<<m_libname<<std::endl;
  sync();
}

std::string  Single_Process::CreateLibName()
{
  string name=m_ptypename;
  name+="_"+ToString(p_ampl->GetGraphNumber());
  name+="_"+ToString(p_shand->NumberOfCouplings());
  name+="_"+ToString(p_shand->NumberOfZfuncs());
  name+="_"+ToString(p_hel->MaxHel());
  name+="_"+ToString(p_BS->MomlistSize());
  return name;
}

void Single_Process::CreateMappingFile(Single_Process* partner) {
  if (m_gen_str<2) return;
  std::string outname = rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename+"/"+m_name+".map";
  if (IsFile(outname)) {
    string MEname,PSname;
    FoundMappingFile(MEname,PSname);
    if (MEname != m_libname || PSname != m_pslibname) {
      msg_Error()<<"ERROR in Single_Process::CreateMappingFile() :"<<std::endl
		 <<"   Files do not coincide. Maybe changed input data ? Abort the run."<<std::endl;
      abort();
    }
    return;
  }

  std::ofstream to;
  to.open(outname.c_str(),ios::out);
  to<<"ME: "<<m_libname<<endl
    <<"PS: "<<m_pslibname<<endl;
  p_shand->Get_Generator()->WriteCouplings(to);
  partner->WriteMomFlavs(to);
  to.close();
}

bool Single_Process::FoundMappingFile(std::string & MEname, std::string & PSname) 
{
  std::string buf;
  int pos;
  std::string outname = rpa.gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_ptypename+"/"+m_name+".map";
  if (IsFile(outname)) {
    ifstream from;
    from.open(outname.c_str());
    getline(from,buf);
    pos = buf.find(string("ME:"));
    if (pos==-1) MEname = PSname = buf;
    else {
      MEname = buf.substr(pos+4);
      getline(from,buf);
      pos = buf.find(string("PS:"));
      if (pos==-1) PSname = MEname;
      else PSname = buf.substr(pos+4);
      if (PSname==string("")) PSname = MEname;
    }
    return 1;
  }
  return 0;
}

void Single_Process::UpdateMappingFile(std::string name, map<string,Complex> & cmap) 
{
  std::string buf;
  int pos;
  name+=".map";
  ifstream from;
  from.open(name.c_str());
  getline(from,buf);
  pos = buf.find(string("ME:"));
  if (pos==-1) m_libname = m_pslibname = buf;
  else {
    m_libname = buf.substr(pos+4);
    getline(from,buf);
    pos = buf.find(string("PS:"));
    if (pos==-1) m_pslibname = m_libname;
    else m_pslibname = buf.substr(pos+4);
    if (m_pslibname==string("")) m_pslibname = m_libname;
  }
  p_shand->Get_Generator()->ReadCouplings(from);
  from.close();
  p_shand->Get_Generator()->UpdateCouplings(cmap);
}

/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/


bool Single_Process::SetUpIntegrator() 
{  
  if (m_nin==2) {
    if ( (p_flavours[0].Mass() != p_isrhandler->Flav(0).Mass()) ||
	 (p_flavours[1].Mass() != p_isrhandler->Flav(1).Mass()) ) p_isrhandler->SetPartonMasses(p_flavours);
    if (CreateChannelLibrary()) return 1;
  }
  if (m_nin==1) if (CreateChannelLibrary()) return 1;
  return 0;
}

bool Single_Process::CreateChannelLibrary()
{
  p_psgen     = new Phase_Space_Generator(m_nin,m_nout);
  bool newch  = 0;
  if (m_nin>=1)  newch = p_psgen->Construct(p_pshandler->GetChannelLibNames(),m_ptypename,m_pslibname,p_flavours,this); 
  if (newch>0) {
    msg_Tracking()<<"Single_Process::CreateChannelLibrary() :"<<std::endl
		  <<"   "<<p_pshandler->NumberOfFSRIntegrators()<<" new channels produced for "
		  <<m_pslibname<<". After program termination please enter \"makelibs\" and rerun."<<endl;
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

  if (p_selector && m_ownselector) { 
    delete p_selector;      
  }
  m_orderQCD      = p_partner->OrderStrong();
  m_orderEW       = p_partner->OrderEWeak();
  p_selector      = p_partner->Selector(); 
  p_jf            = p_partner->p_jf;
  m_ownselector=false;
  if (p_pshandler) { delete p_pshandler; p_pshandler = 0; }
}

void Single_Process::Empty() {
  if (p_pshandler)          { delete p_pshandler; p_pshandler = 0; } 
  if (p_partner != this) {
    return;
  }
}

void Single_Process::SetTotal(int flag, int depth)  { 
  if (flag!=2) {
    m_totalxs  = TotalResult(); 
    m_totalerr = TotalVar();
//     m_totalerr = sqrt( (m_totalsumsqr/m_n - 
// 			(ATOOLS::sqr(m_totalsum)-m_totalsumsqr)/(m_n*(m_n-1.)) )  / m_n); 
    if ((m_nin==1 && m_nout==2) || m_n==1) m_totalerr = 0.;
  }
  else {
    //   flag==2  means  check xs with sum of subprocesses
    //   nothing to do for a Single_Process
  }
  if (m_nin==2 && flag==0) {
    if ( (depth<=0 && msg_LevelIsInfo()) || msg_LevelIsTracking()) {
      for (int i=0;i<depth;++i) msg_Out()<<"  ";
      msg_Out()<<om::bold<<m_name<<om::reset<<" : "
	       <<om::blue<<om::bold<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"<<om::reset
	       <<" +/- "<<om::reset<<om::blue<<m_totalerr/m_totalxs*100.<<" %,"<<om::reset
	       <<om::bold<<" exp. eff: "<<om::red<<(100.*m_totalxs/m_max)<<" %."<<om::reset<<endl;
    }
  }
  if (m_nin==1) {
    msg_Info()<<"Total width for "<<m_name<<" : "
	      <<m_totalxs<<" GeV"
	      <<" +/- "<<m_totalerr/m_totalxs*100.<<"%, max : "<<m_max<<endl;
  }
}

/*------------------------------------------------------------------------------

  Calculating total cross sections
  
  ------------------------------------------------------------------------------*/

bool Single_Process::CalculateTotalXSec(std::string _resdir) { 
  msg_Info()<<"In Single_Process::CalculateTotalXSec("<<_resdir<<") for "<<m_name<<endl; 
  
  string _name;
  double _totalxs,_totalerr,_max,sum,sqrsum,ssum,ssqrsum,ss2,wmin;
  long int n,sn,son;
  std::string filename = _resdir+"/"+m_name+".xs_tot";
  std::string histofile =_resdir+string("/WD_")+m_name+"/";
  if (_resdir!=string("")) {
    if (IsFile(filename)) {
      ifstream from;
      from.open(filename.c_str(),ios::in);
      from>>_name>>_totalxs>>_max>>_totalerr>>sum>>sqrsum>>n>>ssum>>ssqrsum>>ss2>>sn>>wmin>>son;
      if (_name==m_name) {
	m_totalxs  = _totalxs;
	m_totalerr = _totalerr; 
	m_max      = _max;
	m_n        = n;
	m_totalsum = sum;
	m_totalsumsqr = sqrsum;
	m_vsmax.clear(); 
	m_vsmax.push_back(_max);
	m_vsn.clear();   
	m_vsn.push_back(n);
	m_sn = sn; m_smax=0.;
	m_ssum     = ssum;
	m_ssumsqr  = ssqrsum;
	m_ssigma2  = ss2;
	m_wmin     = wmin;
	m_son      = son;
      }
      msg_Tracking()<<"Found result : xs for "<<m_name<<" : "
		    <<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
		    <<" +/- "<<m_totalerr/m_totalxs*100.<<"%, max : "<<m_max<<endl;
      from.close();
      p_pshandler->ReadIn(_resdir+string("/MC_")+m_name);
      if (p_pshandler->BeamIntegrator() != 0) p_pshandler->BeamIntegrator()->Print();
      if (p_pshandler->ISRIntegrator()  != 0) p_pshandler->ISRIntegrator()->Print();
      if (p_pshandler->FSRIntegrator()  != 0) p_pshandler->FSRIntegrator()->Print();
      p_pshandler->InitIncoming();
      ReadInHistogram(histofile);
    }
  }
  m_resultpath=_resdir;
  m_resultfile=filename;
  m_histofile=histofile;
  ATOOLS::exh->AddTerminatorObject(this);
  long unsigned int points=m_n;
  m_totalxs = p_pshandler->Integrate();
  if (m_nin==2) m_totalxs /= ATOOLS::rpa.Picobarn();
  if (!(ATOOLS::IsZero((m_totalxs-TotalResult())/(m_totalxs+TotalResult())))) {
    msg_Error()<<"ERROR in Single_Process::CalculateTotalXSec :"<<std::endl
	       <<"   Result of PS-Integrator and internal summation do not coincide for "<<endl
	       <<m_name<<" : "<<m_totalxs<<" vs. "<<TotalResult()<<endl;
  }
  SetTotal(0);
  if (m_totalxs>=0.) {
    if (points==m_n) {
      ATOOLS::exh->RemoveTerminatorObject(this);
      return 1;
    }
    if (_resdir!=string("")) {
      msg_Info()<<"Store result : xs for "<<m_name<<" : "
		<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
		<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
		<<"       max : "<<m_max<<endl;
      ATOOLS::MakeDir(histofile,0); 
      WriteOutHistogram(histofile);
      p_pshandler->WriteOut(_resdir+string("/MC_")+m_name);
      std::ofstream to;
      to.open(filename.c_str(),ios::out);
      WriteOutXSecs(to);
      to.close();
    }
    ATOOLS::exh->RemoveTerminatorObject(this);
    return 1;
  }
  ATOOLS::exh->RemoveTerminatorObject(this);
  return 0;      
}

void Single_Process::PrepareTerminate()
{
  if (rpa.gen.BatchMode()) return;
  if (m_resultpath.length()==0 && m_resultfile.length()==0) return;
  SetTotal(0);
  ATOOLS::MakeDir(m_histofile,0); 
  WriteOutHistogram(m_histofile);
  msg_Info()<<"Store result : xs for "<<m_name<<" : "
	    <<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb"
	    <<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
	    <<"       max : "<<m_max<<endl;
  p_pshandler->WriteOut(m_resultpath+string("/MC_")+m_name);
  std::ofstream to;
  to.open(m_resultfile.c_str(),ios::out);
  WriteOutXSecs(to);
  to.close();
}

void Single_Process::WriteOutXSecs(std::ofstream & _to)    
{ 
//   _to<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<" "
//      <<m_totalsum<<" "<<m_totalsumsqr<<" "<<m_n<<endl; 
  _to.precision(12);
  _to<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<" "
     <<m_totalsum<<" "<<m_totalsumsqr<<" "<<m_n<<" "
     <<m_ssum<<" "<<m_ssumsqr<<" "<<m_ssigma2<<" "<<m_sn<<" "<<m_wmin<<" "<<m_son<<endl; 
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
    msg_Info()<<m_name<<" : Set total xsec and max at ycut = "<<ycut
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
  string filename = (m_resdir+string("/Tab")+m_name+string("dY_cut")).c_str();

  m_totalxs = p_pshandler->Integrate();
  if (m_nin==2) m_totalxs /= ATOOLS::rpa.Picobarn();

  if (!(ATOOLS::IsZero((m_totalxs-TotalResult())/(m_totalxs+TotalResult())))) {
    msg_Error()<<"ERROR in Single_Process::PrepareXSecTables :"<<std::endl
	       <<"   Result of PS-Integrator and internal summation do not coincide for "
	       <<m_name<<" : "<<m_totalxs<<" vs. "<<TotalResult()<<endl;
  }
  SetTotal(1);
  p_pshandler->WriteOut(m_resdir+string("/MC_")+m_name);
  if (m_totalxs>0.) return 1;
  return 0;
}

void Single_Process::AddPoint(const double value) {
  Integrable_Base::AddPoint(value);

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

void Single_Process::OptimizeResult()
{
//   double ssigma2 = (m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1);
//   m_ssigma2  += 1./ssigma2; 
//   m_totalsum += m_ssum/ssigma2/m_sn;
//   m_totalsumsqr+= m_ssumsqr/ssigma2/m_sn;
  double ssigma2 = ATOOLS::sqr(m_ssum/m_sn)/((m_ssumsqr/m_sn - ATOOLS::sqr(m_ssum/m_sn))/(m_sn-1));
  if (ssigma2>m_wmin) {
    m_ssigma2  += ssigma2; 
    m_totalsum += m_ssum*ssigma2/m_sn;
    m_totalsumsqr+= m_ssumsqr*ssigma2/m_sn;
    m_ssum     = 0.;
    m_ssumsqr  = 0.;
    m_sn       = 0;
    if (ssigma2/m_son>m_wmin) m_wmin = ssigma2/m_son;
    m_son      = 0;
  }
  m_son++;
}

void Single_Process::ResetMax(int flag)
{
  if (flag==0) {
    if (m_vsmax.size()>1) {
      m_vsmax.erase(m_vsmax.begin());
      m_vsn.erase(m_vsn.begin());
    }
    if (m_vsmax.empty()) {
      m_vsmax.push_back(m_max);
      m_vsn.push_back(m_n);
    }
    m_vsmax.back() = ATOOLS::Max(m_smax,m_vsmax.back());
    m_vsn.back()   = m_n;
  }
  else {
    if (flag==2 && m_vsmax.size()==4) {
      m_vsmax.erase(m_vsmax.begin());
      m_vsn.erase(m_vsn.begin());
    }
    m_vsmax.push_back(m_smax);
    m_vsn.push_back(m_n);
    if (flag==2) m_smax = 0.;
  }
  m_max  = 0.;
  for (size_t i=0;i<m_vsmax.size();i++) m_max=ATOOLS::Max(m_max,m_vsmax[i]);
}

double Single_Process::Differential(const ATOOLS::Vec4D* _moms) 
{ 
  return DSigma(_moms,0); 
}

double Single_Process::Differential2() { 
  if (p_isrhandler->On()==0) return 0.;
  return DSigma2(); 
}


double Single_Process::DSigma(const ATOOLS::Vec4D* _moms,bool lookup)
{
  m_last = m_lastdxs = 0.;
  if (m_nin==2) {
    for (size_t i=0;i<m_nin+m_nout;i++) {
      if (_moms[i][0]<p_flavours[i].PSMass()) return m_last = 0.;
    }
  }
  if (m_nin==1) {
    for (size_t i=m_nin;i<m_nin+m_nout;i++) {
      if (_moms[i][0]<p_flavours[i].PSMass()) return m_last = 0.;
    }
  }
  if (p_partner == this) {
    if (m_helsample) {
      if (!m_inithelsample) InitializeHelicityWeights();
      m_lastdxs = operator()(_moms,SelectedHelicity());
    }
    else m_lastdxs = operator()(_moms);
  }
  else {
    if (lookup) m_lastdxs = p_partner->LastXS()*m_sfactor;
           else m_lastdxs = p_partner->operator()(_moms)*m_sfactor;
  }
  if (m_lastdxs <= 0.) return m_lastdxs = m_last = 0.;
  if (m_nin==2) {
#ifdef USING__Threading
    static pthread_mutex_t s_pdf_mtx;
    pthread_mutex_lock(&s_pdf_mtx);
#endif
    m_lastlumi = p_isrhandler->Weight(p_flin);
#ifdef USING__Threading
    pthread_mutex_unlock(&s_pdf_mtx);
#endif
    int    pols[2] = {p_pl[0].type[0],p_pl[1].type[0]};
    double dofs[2] = {p_pl[0].factor[0],p_pl[1].factor[0]};
    if (p_pl[0].num>1) pols[0] = 99;
    if (p_pl[1].num>1) pols[1] = 99;
    m_lastlumi *= p_beamhandler->Weight(pols,dofs);
  }
  else  m_lastlumi = 1.;

  return m_last = m_Norm * m_lastdxs * m_lastlumi*p_partner->KFactor();
}

double Single_Process::DSigma2() { 
  if ((p_flin[0]==p_flin[1]) || (p_isrhandler->On()==0) ) return 0.;
  if (p_partner == this) {
  }
#ifdef USING__Threading
  static pthread_mutex_t s_pdf_mtx;
  pthread_mutex_lock(&s_pdf_mtx);
#endif
  double tmp = m_Norm * m_lastdxs * p_isrhandler->Weight2(p_flin); 
#ifdef USING__Threading
  pthread_mutex_unlock(&s_pdf_mtx);
#endif
  m_last    += tmp*=p_partner->KFactor();
  return tmp;
}

double Single_Process::operator()(const ATOOLS::Vec4D * mom)
{
  double M2(0.);

  p_BS->CalcEtaMu((ATOOLS::Vec4D*)mom);
  p_hel->InitializeSpinorTransformation(p_BS);

  double helvalue;
  if (p_shand->Is_String()) {
    p_shand->Calculate();

    if (p_hel->UseTransformation()) {
      M2 = p_ampl->Zvalue(p_hel);
    } else {
      for (size_t i=0;i<p_hel->MaxHel();i++) {
	if (p_hel->On(i)) {
	  M2 += p_ampl->Differential(i) * p_hel->Multiplicity(i) 
	             * p_hel->PolarizationFactor(i);
	  //	  M2     += helvalue;     
	}
      } 
    }
  }
  else {
    // *** is this ever called ?
    for (size_t i=0;i<p_hel->MaxHel();i++) {
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

// ATOOLS::Spin_Correlation_Tensor* Single_Process::GetSpinCorrelations()
// {
//   Spin_Correlation_Tensor* SCT = p_ampl->GetSpinCorrelations(p_hel, m_nin);
//   if (SCT != NULL) SCT->Set_k0(p_BS->Getk0_n());
//   return SCT;
// }

void Single_Process::FillAmplitudes(HELICITIES::Amplitude_Tensor* atensor,double sfactor)
{
  if (p_partner==this) p_ampl->FillAmplitudes(atensor,p_hel,sfactor);
  else p_partner->FillAmplitudes(atensor,sfactor*sqrt(m_sfactor));
}

ATOOLS::Blob_Data_Base *Single_Process::OneEvent(double _mass) { 
  if (p_partner==this) return p_pshandler->OneEvent(_mass,1); 
  return p_partner->OneEvent(_mass);
}

ATOOLS::Blob_Data_Base *Single_Process::SameEvent() { 
  if (p_partner==this) return p_pshandler->SameEvent(); 
  return p_partner->SameEvent(); 
}

ATOOLS::Blob_Data_Base * Single_Process::WeightedEvent(const int mode)     
{ 
  if (p_partner==this) return p_pshandler->WeightedEvent(mode);
  return p_partner->WeightedEvent(mode); 
}

ATOOLS::Blob_Data_Base * Single_Process::SameWeightedEvent() 
{ 
  if (p_partner==this) return p_pshandler->SameWeightedEvent(); 
  return p_partner->SameWeightedEvent(); 
}




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
  if (!(msg_LevelIsDebugging())) return;
  msg_Out()<<m_name<<" : "<<m_last<<" -> "
		      <<m_lastdxs<<" @ "<<m_lastlumi<<", "<<endl;
}


/*------------------------------------------------------------------------------
  
  Stuff to sample over helicities

  ------------------------------------------------------------------------------*/

void Single_Process::InitializeHelicityWeights()
{
  int activehels = 0, active = 0;
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (p_hel->On(i)) {
      active++,
      activehels += p_hel->Multiplicity(i);
    }
  }

  double alpha_start = 1./activehels;
  double alpha;
  for (size_t i=0;i<p_hel->MaxHel();i++) {
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
    msg_Out()<<"WARNING in Single_Process::SelectedHelicity() after "<<m_throws<<std::endl;
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

  msg_Tracking()<<"After Optimize Helicity Weights for "<<m_name<<" after "<<m_throws<<" throws: "<<std::endl; 
  for (i=0;i<m_helnumber;i++) {
    msg_Tracking()<<i<<" th helicity: "<<m_helnumbers[i]<<"("<<p_hel->Multiplicity(m_helnumbers[i])
		  <<"), alpha = "<<m_helalphas[i]<<std::endl;
  }

  double variance = (m_throws*m_helresult2)/((m_throws-1)*ATOOLS::sqr(m_helresult)) - 1./(m_throws-1);
  if (variance>0.) variance = m_helresult/m_throws * sqrt(variance);
  else {
    msg_Error()<<"Negative variance."<<std::endl;
    variance = m_helresult/m_throws * sqrt(-variance);
  }

  msg_Tracking()<<"S1X: "<<s1x<<", variance : "<<variance<<std::endl
		<<"result,result2,n: "<<m_helresult<<", "
		<<m_helresult2<<", "<<m_throws<<std::endl
		<<"-----------------------------------------------"<<endl;
}


double Single_Process::operator()(const ATOOLS::Vec4D * mom,const int hel)
{
  if (!p_shand->Is_String()) {
    msg_Error()<<"Error in Single_Process::operator()(ATOOLS::Vec4D * p,int hel)"<<std::endl
	       <<"   Sampling over helicities in processes implemented only for libs."<<std::endl
	       <<"   Will abort the run. Check for libraries."<<std::endl;
    abort();
  }

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

bool Single_Process::CheckMapping(const Process_Base * proc)
{
  const Flavour * flavs=Flavours();
  const Flavour * partner_flavs=proc->Flavours();
  // create map
  std::map<ATOOLS::Flavour,ATOOLS::Flavour> flmap;
  for (size_t i=0;i<NIn()+NOut();++i) {
    if (flmap.find(partner_flavs[i])==flmap.end()) {
      flmap[partner_flavs[i]]=flavs[i];
      if (partner_flavs[i]!=(Flavour(partner_flavs[i])).Bar()) {
	flmap[(Flavour(partner_flavs[i])).Bar()]=(Flavour(flavs[i])).Bar();
      }
    }
  }
  // check map
  for (size_t i=0;i<NIn()+NOut();++i) {
    if (flmap[partner_flavs[i]]!=flavs[i]) {
      msg_Tracking()<<" mapping test failed "<<std::endl;
      return false;
    }
  }
  return true;
}

bool             Single_Process::SelectOne()                        { return true;          }
bool             Single_Process::SelectOneFromList()                { return true;          }
void             Single_Process::DeSelect()                         {                       }

bool             Single_Process::ReSelect(int i)                    
{   
  if (i==1) return true;
  if (p_parent==Parent()) return false;
  return Parent()->SelectOne();
}

size_t           Single_Process::Size()                             { return 1;             }
Process_Base   * Single_Process::operator[] (int idx)               { return this;          }
