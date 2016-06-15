#include "AMEGIC++/Main/Single_Process.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "AMEGIC++/Phasespace/Phase_Space_Generator.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_File.H"

#include <unistd.h>

using namespace AMEGIC;
using namespace MODEL;
using namespace PHASIC;
using namespace PDF;
using namespace BEAM;
using namespace ATOOLS;
using namespace std;

/*-------------------------------------------------------------------------------

  Constructors

  ------------------------------------------------------------------------------- */

AMEGIC::Single_Process::Single_Process():
  m_gen_str(2), p_hel(0), p_BS(0), p_ampl(0), p_shand(0), p_psgen(0), p_partner(this)
{
  m_lastk=1.0;
}

AMEGIC::Single_Process::~Single_Process()
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

void AMEGIC::Single_Process::PolarizationNorm() {
  m_Norm = SymmetryFactors() * m_pol.Spin_Average(m_nin,&m_flavs.front());
}


void AMEGIC::Single_Process::WriteAlternativeName(string aname) 
{
  if (aname==m_name) return;
  std::string altname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_name+".alt";
  if (FileExists(altname,1)) return;
  My_Out_File to(altname);
  to.Open();
  *to<<aname<<" "<<m_sfactor<<endl;
  for (map<string,ATOOLS::Flavour>::const_iterator fit=p_ampl->GetFlavourmap().begin();fit!=p_ampl->GetFlavourmap().end();fit++)
    *to<<fit->first<<" "<<(long int)fit->second<<endl;
  to.Close();
}

bool AMEGIC::Single_Process::CheckAlternatives(vector<Process_Base *>& links,string procname)
{
  std::string altname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+procname+".alt";
  if (FileExists(altname,1)) {
    double factor;
    string name,dummy; 
    My_In_File from(altname);
    from.Open();
    *from>>name>>factor;
    m_sfactor *= factor;
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (links[j]->Name()==name) {
	p_mapproc = p_partner = (Single_Process*)links[j];
	m_iresult = p_partner->Result()*m_sfactor;
	m_maxcpl=p_partner->MaxOrders();
	m_mincpl=p_partner->MinOrders();
	m_ntchanmin=p_partner->NTchanMin();
	msg_Tracking()<<"Found Alternative process: "<<m_name<<" "<<name<<endl;

	while (*from) {
	  string f1;
	  long int f2;
	  getline(*from,dummy);
	  if (dummy!="") {
	    MyStrStream str;
	    str<<dummy;
	    str>>f1>>f2;
	    AddtoFlavmap(f1,Flavour(abs(f2),f2<0));
	  }
	}
	from.Close();
	InitFlavmap(p_partner);
	FillCombinations();
	return true;
      }
    }
    from.Close();
    if (CheckAlternatives(links,name)) return true;
  }
  m_sfactor = 1.;
  return false;
}

int AMEGIC::Single_Process::InitAmplitude(Amegic_Model *model,Topology* top,
					  vector<Process_Base *> & links,
					  vector<Process_Base *> & errs)
{
  Init();
  model->p_model->GetCouplings(m_cpls);
  if (!model->p_model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;
  m_newlib   = false;
  m_pslibname = m_libname = ToString(m_nin)+"_"+ToString(m_nout);
  if (m_gen_str>1) m_ptypename = "P"+m_libname;
  else m_ptypename = "N"+m_libname;
  PolarizationNorm();

  if (CheckAlternatives(links,Name())) return 1;

  p_hel    = new Helicity(m_nin,m_nout,&m_flavs.front(),p_pl);

  bool directload = true;
  int libchk=ToType<int>(rpa->gen.Variable("ME_LIBCHECK")); 
    if (libchk==1) {
      msg_Info()<<"Enforce full library check. This may take some time"<<std::endl;
      directload = false;
    }
  if (directload) directload = FoundMappingFile(m_libname,m_pslibname);
  if (directload) {
    string hstr=rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_libname;
    string hstr2=rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_name+".map";
    p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nin+m_nout,&m_flavs.front(),p_b,hstr,hstr2);  
  }
  else p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nin+m_nout,&m_flavs.front(),p_b);  
  p_BS->Setk0(s_gauge);
  p_shand  = new String_Handler(m_gen_str,p_BS,model->p_model->GetCouplings());
  int ntchanmin(m_ntchanmin);
  bool cvp(ToType<int>(rpa->gen.Variable("AMEGIC_CUT_MASSIVE_VECTOR_PROPAGATORS")));
  p_ampl   = new Amplitude_Handler(m_nin+m_nout,&m_flavs.front(),p_b,p_pinfo,model,top,m_maxcpl,m_mincpl,ntchanmin,
                                   &m_cpls,p_BS,p_shand,m_print_graphs,!directload,cvp);
  if (p_ampl->GetGraphNumber()==0) {
    msg_Tracking()<<"AMEGIC::Single_Process::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return 0;
  }
  if (!directload) {
  map<string,Complex> cplmap;
  for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
    cplmap.clear();
    if (m_pinfo.m_special.find("MapOff")!=std::string::npos) continue;
    if (m_allowmap && FlavCompare(links[j]) && p_ampl->CompareAmplitudes(links[j]->GetAmplitudeHandler(),m_sfactor,cplmap)) {
      if (p_hel->Compare(links[j]->GetHelicity(),m_nin+m_nout)) {
	m_sfactor = sqr(m_sfactor);
	msg_Tracking()<<"AMEGIC::Single_Process::InitAmplitude : Found compatible process for "<<Name()<<" : "<<links[j]->Name()<<endl;
	  
	bool found(true);
	if (!FoundMappingFile(m_libname,m_pslibname)) {
	  string mlname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+links[j]->Name();
	  string mnname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+Name();
	  if (FileExists(mlname+string(".map"),1)) { 
	    if (m_sfactor==1.) My_In_File::CopyInDB(mlname+".map",mnname+".map");
	    else {
	      UpdateMappingFile(mlname,cplmap);
	      CreateMappingFile((Single_Process*)links[j]);
	    }
            My_In_File::CopyInDB(mlname+".col",mnname+".col");
	  }
	  found=false;
	}

	p_mapproc = p_partner = (Single_Process*)links[j];
	m_iresult = p_partner->Result()*m_sfactor;
	for (std::map<string,Flavour>::const_iterator fit=p_ampl->GetFlavourmap().begin();
	     fit!=p_ampl->GetFlavourmap().end();fit++) AddtoFlavmap(fit->first,fit->second);
	InitFlavmap(p_partner);
	FillCombinations();
	if (!found) WriteAlternativeName(p_partner->Name());

	Minimize();
	return 1;
      }
    } 
  }
  }
  if (directload) {
    p_ampl->CompleteLibAmplitudes(m_nin+m_nout,m_ptypename+string("/")+m_name,
				  m_ptypename+string("/")+m_libname);    
    if (!p_shand->SearchValues(m_gen_str,m_libname,p_BS)) {
      m_newlib=true;
      return 1;
    }
    if (!TestLib()) return 0;
    links.push_back(this);
    FillCombinations();
    msg_Info()<<"."<<std::flush;
    Minimize();
    return 1;
  }

  p_ampl->CompleteAmplitudes(m_nin+m_nout,&m_flavs.front(),p_b,&m_pol,
			     top,p_BS,m_ptypename+string("/")+m_name);
  m_pol.Add_Extern_Polarisations(p_BS,&m_flavs.front(),p_hel);
  p_BS->Initialize();
  FillCombinations();

  int result(Tests());
  switch (result) {
    case 2 : 
    if (p_partner==this) links.push_back(this);
    Minimize();
    WriteAlternativeName(p_partner->Name());
    return 1;
  case 1 :
    if (Result()==0.) return -3;
    if (p_partner==this) links.push_back(this);

    if (CheckLibraries()) return 1;
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (links[j]->NewLibs()) {
	if (CheckStrings((Single_Process*)links[j])) return 1;	
      }      
    }
    if (p_partner!=this) links.push_back(this);

    if (m_gen_str<2) return 1;
    if (p_partner!=this) {
      msg_Tracking()<<"AMEGIC::Single_Process::InitAmplitude : "<<std::endl
		    <<"   Strings of process "<<m_name<<" and partner "
		    <<p_partner->Name()<<" did not fit."<<std::endl
		    <<"   Have to write new library."<<std::endl;
    }
    WriteLibrary();
    if (p_partner==this && Result()>0.) SetUpIntegrator();
    return 1;
  case -3: return 0;
  default :
    msg_Error()<<"ERROR in AMEGIC::Single_Process::InitAmplitude : "<<std::endl
	       <<"   Failed for "<<m_name<<" with result "<<result<<"."<<endl;
    errs.push_back(this);
    return 0;
  }

  return 1;
}



int AMEGIC::Single_Process::Tests() 
{
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
    m_pol.Set_Gauge_Vectors(m_nin+m_nout,p_testmoms,Vec4D(sqrt(3.),1.,1.,-1.));
    p_BS->Setk0(0);
    p_BS->CalcEtaMu(p_testmoms);  
    p_BS->InitGaugeTest(.9);

    msg_Info()<<"AMEGIC::Single_Process::Tests for "<<m_name<<std::endl
	      <<"   Prepare gauge test and init helicity amplitudes. This may take some time."
	      <<std::endl;
    for (size_t i=0;i<p_hel->MaxHel();i++) { 
      if (p_hel->On(i)) {
	helvalue = p_ampl->Differential(i,(*p_hel)[i])*p_hel->PolarizationFactor(i); 
	M2      +=  helvalue;
      } 
    }
    M2     *= sqr(m_pol.Massless_Norm(m_nin+m_nout,&m_flavs.front(),p_BS));
    m_iresult  = M2;
  }
  p_ampl->ClearCalcList();
  // To prepare for the string test.
  p_ampl->SetStringOn();
  (p_shand->Get_Generator())->Reset(1);
  /* ---------------------------------------------------
     
  First test : gauge test
  
  --------------------------------------------------- */
  p_BS->Setk0(s_gauge);
  p_BS->CalcEtaMu(p_testmoms);
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
    if (M_doub[i]==0.) {
#ifdef FuckUp_Helicity_Mapping
      p_hel->SwitchOff(i);
      switchhit++;
#endif
    }
  }
  msg_Tracking()<<"AMEGIC::Single_Process::Tests for "<<m_name<<std::endl
		<<"   Switched off or mapped "<<switchhit<<" helicities."<<std::endl;

  M2g    *= sqr(m_pol.Massless_Norm(m_nin+m_nout,&m_flavs.front(),p_BS));
  m_iresult  = M2g;
  p_ampl->ClearCalcList();  
  p_ampl->FillCoupling(p_shand);
  p_ampl->KillZList();  
  p_BS->StartPrecalc();

  if (gauge_test) {
    if (!ATOOLS::IsEqual(M2,M2g)) {
      msg_Info()<<"WARNING:  Gauge test not satisfied: "
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
      M2 = operator()(p_testmoms);
      p_hel->AllowTransformation();
      string_test = 0;
    }
    else {
      string searchfilename = rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+m_ptypename+string("/")+testname+string("/V.H");
      if (FileExists(searchfilename)) {
      	msg_Error()<<"ERROR in AMEGIC::Single_Process::Tests()"<<std::endl
			   <<"   No compiled & linked library found for process "<<testname<<std::endl
			   <<"   but files already written out !"<<std::endl
			   <<om::bold<<"   Interrupt run and execute \"makelibs\" in '"
			   <<rpa->gen.Variable("SHERPA_CPP_PATH")<<"'."
			   <<om::reset<<std::endl;
	Copy(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/makelibs",
		 rpa->gen.Variable("SHERPA_CPP_PATH")+"/makelibs");
	THROW(normal_exit,"Failed to load library.");
      }
      else {
      	msg_Error()<<"ERROR in AMEGIC::Single_Process::Tests()"<<std::endl
			   <<"   Mapping file exists, but no compiled & linked library found for process "
			   <<testname<<std::endl
			   <<"   and no files written out !"<<std::endl
			   <<om::bold<<"   Interrupt run, execute \"makeclean\" in Run-directory and re-start."
			   <<om::reset<<std::endl;
	THROW(critical_error,"Failed to load library.");
      }
    }
    if (!ATOOLS::IsEqual(M2,M2g)) {
      if (abs(M2/M2g-1.)>rpa->gen.Accu()) {
	msg_Info()<<"WARNING: Library cross check not satisfied: "
		 <<M2<<" vs. "<<M2g<<"  difference:"<<abs(M2/M2g-1.)*100.<<"%"<<endl
		 <<"   Mapping file(1) : "<<abs(M2)<<endl
		 <<"   Original    (2) : "<<abs(M2g)<<endl
		 <<"   Cross check (T) : "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
	THROW(critical_error,"Check output above. Increase NUM_ACCURACY if you wish to skip this test.");
      }
      else {
	msg_Info()<<"WARNING: Library cross check not satisfied: "
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
    if (!(rpa->gen.SoftSC()||rpa->gen.HardSC())) {
      for (size_t i=0;i<p_hel->MaxHel();i++) {
	if (p_hel->On(i)) {
	  for (size_t j=i+1;j<p_hel->MaxHel();j++) {
	    if (p_hel->On(j)) {
#ifdef FuckUp_Helicity_Mapping
	      if (ATOOLS::IsEqual(M_doub[i],M_doub[j])) {
		p_hel->SwitchOff(j);
		p_hel->SetPartner(i,j);
		p_hel->IncMultiplicity(i);
	      }
#endif
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
      M2S *= sqr(m_pol.Massless_Norm(m_nin+m_nout,&m_flavs.front(),p_BS));
      if (!ATOOLS::IsEqual(M2g,M2S)) {
	msg_Info()<<"WARNING: String test not satisfied: "
		 <<M2g<<" vs. "<<M2S<<"  difference:"<<abs(M2g/M2S-1.)*100.<<"%"<<endl;
	if (abs(M2g/M2S-1.)>rpa->gen.Accu()) {
	  THROW(critical_error,"Check output above. Increase NUM_ACCURACY if you wish to skip this test.");
	}
	msg_Info()<<"         assuming numerical reasons, continuing "<<endl;
      }
      return 1;
    }
    return 1;
  }
  delete[] M_doub;

  return 0;
}

int AMEGIC::Single_Process::TestLib()
{
  double M2(0.);
  double * M_doub = new double[p_hel->MaxHel()];
  p_BS->CalcEtaMu(p_testmoms);
  p_hel->InitializeSpinorTransformation(p_BS);
  p_shand->Calculate();
  
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    M2 += M_doub[i] = p_ampl->Differential(i)*p_hel->Multiplicity(i)*p_hel->PolarizationFactor(i);
  } 
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (M_doub[i]==0.) {
#ifdef FuckUp_Helicity_Mapping
     p_hel->SwitchOff(i);
#endif
    }
  }
  if (!(rpa->gen.SoftSC()||rpa->gen.HardSC())) {
    for (size_t i=0;i<p_hel->MaxHel();i++) {
      if (p_hel->On(i)) {
	for (size_t j=i+1;j<p_hel->MaxHel();j++) {
	  if (p_hel->On(j)) {
#ifdef FuckUp_Helicity_Mapping
	    if (ATOOLS::IsEqual(M_doub[i],M_doub[j])) {
	      p_hel->SwitchOff(j);
	      p_hel->SetPartner(i,j);
	      p_hel->IncMultiplicity(i);
	    }
#endif
	  }
	}
      }
    }
  }
  delete[] M_doub;
  m_iresult = M2 * sqr(m_pol.Massless_Norm(m_nin+m_nout,&m_flavs.front(),p_BS));
  if (m_iresult>0. || m_iresult<0.) return 1;
  return 0;
}

int AMEGIC::Single_Process::CheckLibraries() {
  if (m_gen_str==0) return 1;
  if (p_shand->IsLibrary()) return 1;

  msg_Info()<<"AMEGIC::Single_Process::CheckLibraries : Looking for a suitable library. This may take some time."<<std::endl;
  String_Handler * shand1;
  shand1      = new String_Handler(p_shand->Get_Generator());
  
  string proc = rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+m_ptypename+string("/V");
  string testname;
  double M2s, helvalue;

  for (;;) {
    testname  = CreateLibName();
    if (shand1->SearchValues(m_gen_str,testname,p_BS)) {
      shand1->Calculate();
      
      M2s = 0.;
      for (size_t i=0;i<p_hel->MaxHel();i++) {
	helvalue = p_ampl->Differential(shand1,i) * p_hel->PolarizationFactor(i) *
	  p_hel->Multiplicity(i);
	M2s     += helvalue;
	} 
      M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,&m_flavs.front(),p_BS));
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
  }
  if (shand1) { delete shand1; shand1 = 0; }
  return 0;
}

int AMEGIC::Single_Process::CheckStrings(Single_Process* tproc)
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
  M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,&m_flavs.front(),p_BS));
  (shand1->Get_Generator())->ReStore();
  delete shand1;

  if (ATOOLS::IsEqual(M2s,Result())) {
    m_newlib = true;
    m_libname = tproc->LibName();
    m_pslibname = tproc->PSLibName();
    CreateMappingFile(this);
    Minimize();
    return 1;
  }
  return 0;
}
  
void AMEGIC::Single_Process::WriteLibrary() 
{
  if (m_gen_str<2) return;
  string newpath=rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/");
  m_libname = CreateLibName();
  if (p_partner==this) m_pslibname = m_libname;
                  else m_pslibname = p_partner->PSLibName();
  if (!FileExists(newpath+m_ptypename+string("/")+m_libname+string("/V.H"))) {
  ATOOLS::MakeDir(newpath+m_ptypename+"/"+m_libname); 
  p_shand->Output(p_hel,m_ptypename+string("/")+m_libname);
  }
  CreateMappingFile(this);
  p_BS->Output(newpath+m_ptypename+string("/")+m_libname);
  p_ampl->StoreAmplitudeConfiguration(newpath+m_ptypename+string("/")+m_libname);
  m_newlib=true;
  if (!FileExists(rpa->gen.Variable("SHERPA_CPP_PATH")+"/makelibs"))
    Copy(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/makelibs",
	 rpa->gen.Variable("SHERPA_CPP_PATH")+"/makelibs");
  msg_Info()<<"AMEGIC::Single_Process::WriteLibrary : "<<std::endl
	    <<"   Library for "<<m_name<<" has been written, name is "<<m_libname<<std::endl;
  sync();
}

void AMEGIC::Single_Process::CreateMappingFile(Single_Process* partner) {
  if (m_gen_str<2) return;
  std::string outname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_name+".map";
  if (FileExists(outname,1)) {
    string MEname,PSname;
    FoundMappingFile(MEname,PSname);
    if (MEname != m_libname || PSname != m_pslibname) {
      msg_Error()<<"ERROR in AMEGIC::Single_Process::CreateMappingFile() :"<<std::endl
		 <<"   Files do not coincide. Maybe changed input data ? Abort the run."<<std::endl;
      abort();
    }
    return;
  }

  My_Out_File to(outname);
  to.Open();
  *to<<"ME: "<<m_libname<<endl
    <<"PS: "<<m_pslibname<<endl;
  p_shand->Get_Generator()->WriteCouplings(*to);
  to.Close();
}

bool AMEGIC::Single_Process::FoundMappingFile(std::string & MEname, std::string & PSname) 
{
  std::string buf;
  int pos;
  std::string outname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_name+".map";
  if (FileExists(outname,1)) {
    My_In_File from(outname);
    from.Open();
    getline(*from,buf);
    pos = buf.find(string("ME:"));
    if (pos==-1) MEname = PSname = buf;
    else {
      MEname = buf.substr(pos+4);
      getline(*from,buf);
      pos = buf.find(string("PS:"));
      if (pos==-1) PSname = MEname;
      else PSname = buf.substr(pos+4);
      if (PSname==string("")) PSname = MEname;
    }
    return 1;
  }
  return 0;
}

void AMEGIC::Single_Process::UpdateMappingFile(std::string name, map<string,Complex> & cmap) 
{
  std::string buf;
  int pos;
  name+=".map";
  My_In_File from(name);
  from.Open();
  getline(*from,buf);
  pos = buf.find(string("ME:"));
  if (pos==-1) m_libname = m_pslibname = buf;
  else {
    m_libname = buf.substr(pos+4);
    getline(*from,buf);
    pos = buf.find(string("PS:"));
    if (pos==-1) m_pslibname = m_libname;
    else m_pslibname = buf.substr(pos+4);
    if (m_pslibname==string("")) m_pslibname = m_libname;
  }
  p_shand->Get_Generator()->ReadCouplings(*from);
  from.Close();
  p_shand->Get_Generator()->UpdateCouplings(cmap);
}

int AMEGIC::Single_Process::PerformTests()
{
  return 1;
}

/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/

bool AMEGIC::Single_Process::FillIntegrator
(PHASIC::Phase_Space_Handler *const psh)
{
  if (p_partner!=this) return true;
  if (!SetUpIntegrator()) THROW(fatal_error,"No integrator");
  return Process_Base::FillIntegrator(psh);
}

bool AMEGIC::Single_Process::SetUpIntegrator() 
{  
  if (m_nin==2) {
    if ( (m_flavs[0].Mass() != p_int->ISR()->Flav(0).Mass()) ||
	 (m_flavs[1].Mass() != p_int->ISR()->Flav(1).Mass()) ) p_int->ISR()->SetPartonMasses(m_flavs);
    if (CreateChannelLibrary()) return 1;
  }
  if (m_nin==1) if (CreateChannelLibrary()) return 1;
  return 0;
}

bool AMEGIC::Single_Process::CreateChannelLibrary()
{
  if (p_partner!=this || p_psgen) return true;
  p_psgen     = new Phase_Space_Generator(m_nin,m_nout);
  bool newch  = 0;
  if (m_nin>=1)  newch = p_psgen->Construct(p_channellibnames,m_ptypename,m_pslibname,&m_flavs.front(),this); 
  if (newch>0) return 0;
  return 1;
}

/*------------------------------------------------------------------------------
  
  Process management
  
  ------------------------------------------------------------------------------*/
void AMEGIC::Single_Process::Minimize()
{
  if (p_partner==this) return;
  if (p_hel)      {delete p_hel; p_hel=0;}
  if (p_BS)       {delete p_BS;   p_BS=0;}
  if (p_shand)    {delete p_shand;p_shand=0;}
  if (p_ampl)     {delete p_ampl; p_ampl=0;}
  if (p_psgen)    {delete p_psgen; p_psgen=0;}

  m_maxcpl    = p_partner->MaxOrders();
  m_mincpl    = p_partner->MinOrders();
  m_ntchanmin = p_partner->NTchanMin();
}

double AMEGIC::Single_Process::Partonic(const Vec4D_Vector &_moms,const int mode) 
{ 
  if (mode==1) return m_mewgtinfo.m_B=m_lastxs;
  if (!Selector()->Result()) return m_mewgtinfo.m_B=m_lastxs=0.0;
  if (!(IsMapped() && LookUp())) {
    p_partner->ScaleSetter()->CalculateScale(_moms,m_cmode);
  }
  m_mewgtinfo.m_K=p_partner->LastK();
  return m_mewgtinfo.m_B=DSigma(_moms,m_lookup);
}

double AMEGIC::Single_Process::DSigma(const ATOOLS::Vec4D_Vector &mom,bool lookup)
{
  m_lastxs = 0.;
  if (p_partner == this) {
    m_lastxs = m_Norm * operator()((ATOOLS::Vec4D*)&mom.front());
  }
  
  else {
    if (lookup && p_partner->m_lookup)
      m_lastxs = p_partner->LastXS()*m_sfactor;
    else m_lastxs = m_Norm * p_partner->operator()((ATOOLS::Vec4D*)&mom.front())*m_sfactor;
  }
  return m_lastxs;
}

double AMEGIC::Single_Process::operator()(const ATOOLS::Vec4D* mom)
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
  m_lastk=KFactor();
  return M2 * sqr(m_pol.Massless_Norm(m_nin+m_nout,&m_flavs.front(),p_BS)) * m_lastk;
}

void AMEGIC::Single_Process::FillAmplitudes(vector<METOOLS::Spin_Amplitudes>& amps,
                                            std::vector<std::vector<Complex> >& cols)
{
  if (p_partner==this) p_ampl->FillAmplitudes(amps, cols, p_hel, 1.0);
  else p_partner->FillAmplitudes(amps, cols, sqrt(m_sfactor));
}

void AMEGIC::Single_Process::FillAmplitudes(vector<METOOLS::Spin_Amplitudes>& amps,
                                            std::vector<std::vector<Complex> >& cols,
                                            double sfactor)
{
  if (p_partner==this) p_ampl->FillAmplitudes(amps, cols, p_hel, sfactor);
  else p_partner->FillAmplitudes(amps, cols, sfactor*sqrt(m_sfactor));
}

int AMEGIC::Single_Process::NumberOfDiagrams() { 
  if (p_partner==this) return p_ampl->GetGraphNumber(); 
  return p_partner->NumberOfDiagrams();
}

Point * AMEGIC::Single_Process::Diagram(int i) { 
  if (p_partner==this) return p_ampl->GetPointlist(i); 
  return p_partner->Diagram(i);
} 


void AMEGIC::Single_Process::AddChannels(std::list<std::string>* tlist) 
{
  if (p_partner==this) {    
    list<string>* clist = p_channellibnames;
    for (list<string>::iterator it=clist->begin();it!=clist->end();++it) {
      bool hit = 0;
      for (list<string>::iterator jt=tlist->begin();jt!=tlist->end();++jt) {
	if ((*it)==(*jt)) {
	  hit = 1;
	  break;
	}
      }
      if (!hit) tlist->push_back((*it));
    }
  }
}

void AMEGIC::Single_Process::FillCombinations
(Point *const p,size_t &id)
{
  if (p->middle) return;
  if (p->left==NULL || p->right==NULL) {
    id=1<<p->number;
    return;
  }
  size_t ida(id), idb(id);
  FillCombinations(p->left,ida);
  FillCombinations(p->right,idb);
  id=ida+idb;
  size_t idc((1<<(m_nin+m_nout))-1-id);
#ifdef DEBUG__Fill_Combinations
  msg_Debugging()<<"  comb "<<ID(ida)
		 <<" "<<ID(idb)<<" "<<ID(idc)<<"\n";
#endif
  m_ccombs.insert(std::pair<size_t,size_t>(ida,idb));
  m_ccombs.insert(std::pair<size_t,size_t>(idb,ida));
  m_ccombs.insert(std::pair<size_t,size_t>(idb,idc));
  m_ccombs.insert(std::pair<size_t,size_t>(idc,idb));
  m_ccombs.insert(std::pair<size_t,size_t>(idc,ida));
  m_ccombs.insert(std::pair<size_t,size_t>(ida,idc));
  if (idc!=1) {
    bool in(false);
    Flavour fl(ReMap(p->fl,p->GetPropID()));
    Flavour_Vector cf(m_cflavs[id]);
    for (size_t i(0);i<cf.size();++i)
      if (cf[i]==fl) {
	in=true;
	break;
      }
    if (!in) {
      m_cflavs[idc].push_back(fl.Bar());
      m_cflavs[id].push_back(fl);
#ifdef DEBUG__Fill_Combinations
      msg_Debugging()<<"  flav "<<ID(idc)<<" / "
		     <<ID(id)<<" -> "<<fl<<"\n";
#endif
    }
  }
}

void AMEGIC::Single_Process::FillCombinations()
{
#ifdef DEBUG__Fill_Combinations
  msg_Debugging()<<METHOD<<"(): '"<<m_name<<"' {\n";
#endif
  size_t nd(p_partner->NumberOfDiagrams());
  for (size_t i(0);i<nd;++i) {
    Point *p(p_partner->Diagram(i));
    size_t id(1<<p->number);
    FillCombinations(p,id);
  }
#ifdef DEBUG__Fill_Combinations
  msg_Debugging()<<"  } -> "<<m_cflavs.size()
		 <<" flavours, "<<m_ccombs.size()
		 <<" combinations\n";
  msg_Debugging()<<"}\n";
#endif
}

bool AMEGIC::Single_Process::Combinable
(const size_t &idi,const size_t &idj)
{
  Combination_Set::const_iterator 
    cit(m_ccombs.find(std::pair<size_t,size_t>(idi,idj)));
  return cit!=m_ccombs.end();
}

const Flavour_Vector &AMEGIC::Single_Process::
CombinedFlavour(const size_t &idij)
{
  CFlavVector_Map::const_iterator fit(m_cflavs.find(idij));
  if (fit==m_cflavs.end()) THROW(fatal_error,"Invalid request");
  return fit->second;
}
