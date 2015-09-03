#include "AMEGIC++/DipoleSubtraction/Single_LOProcess.H"

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


Single_LOProcess::Single_LOProcess(const Process_Info &pi,
                                   BEAM::Beam_Spectra_Handler *const beam,
                                   PDF::ISR_Handler *const isr) :   
  m_gen_str(2), m_emit(-1), m_spect(-1),
  p_hel(0), p_BS(0), p_ampl(0), p_shand(0), p_partner(this), p_sub(NULL)
{
  m_nin=pi.m_ii.NExternal();
  m_nout=pi.m_fi.NExternal();

  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(rpa->GetPath());
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  int ord=reader.GetValue<int>("AMEGIC_SORT_LOPROCESS",1);
  static bool print(false);
  if (!print && !ord) {
    print=true;
    msg_Info()<<METHOD<<"(): "<<om::red
	      <<"Sorting flavors!\n"<<om::reset;    
  }
  PHASIC::Process_Base::Init(pi, beam, isr, ord?0:1);
  AMEGIC::Process_Base::Init();

  m_rsmap.resize(m_nin+m_nout);
  m_srmap.resize(m_nin+m_nout+1,-1);
  for (size_t i(0);i<m_nin;++i) {
    m_rsmap[i]=m_pinfo.m_ii.m_ps[i].m_tag;
    if (m_rsmap[i]>=0) m_srmap[m_rsmap[i]]=i;
  }
  vector<int> fi_tags;
  m_pinfo.m_fi.GetTags(fi_tags);
  if (fi_tags.size()!=m_nout) THROW(fatal_error, "Internal error.");
  for (size_t i(0);i<m_nout;++i) {
    m_rsmap[m_nin+i]=fi_tags[i];
    if (m_rsmap[m_nin+i]>=0) m_srmap[m_rsmap[m_nin+i]]=m_nin+i;
  }

  m_newlib   = false;
  m_pslibname = m_libname = ToString(m_nin)+"_"+ToString(m_nout);
  if (m_gen_str>1) m_ptypename = "P"+m_libname;
  else m_ptypename = "N"+m_libname;

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
  for (size_t i(0);i<m_nout;++i) {
    if (fi_tags[i]==-1) {
      m_emit=i+NIn();
      cnt++;
    }
    if (fi_tags[i]==-2) {
      m_spect=i+NIn();
      cnt+=10;
    }
  }
  if (cnt!=0&&cnt!=11) THROW(critical_error,"mistagged process "+m_name);
}


Single_LOProcess::~Single_LOProcess()
{
  if (p_hel)      {delete p_hel; p_hel=0;}
  if (p_BS)       {delete p_BS;   p_BS=0;}
  if (p_shand)    {delete p_shand;p_shand=0;}
  if (p_ampl)     {delete p_ampl; p_ampl=0;}
}


/*------------------------------------------------------------------------------

  Initializing libraries, amplitudes, etc.

  ------------------------------------------------------------------------------*/


void AMEGIC::Single_LOProcess::WriteAlternativeName(string aname) 
{
  if (aname==m_name) return;
  std::string altname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_name+".alt";
  if (FileExists(altname)) return;
  My_Out_File to(altname);
  to.Open();
  *to<<aname<<" "<<m_sfactor<<endl;
  for (map<string,ATOOLS::Flavour>::const_iterator fit=p_ampl->GetFlavourmap().begin();fit!=p_ampl->GetFlavourmap().end();fit++)
    *to<<fit->first<<" "<<(long int)fit->second<<endl;
  to.Close();
}

bool AMEGIC::Single_LOProcess::CheckAlternatives(vector<Process_Base *>& links,string procname)
{
  std::string altname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+procname+".alt";
  if (FileExists(altname)) {
    double factor;
    string name,dummy; 
    My_In_File from(altname);
    from.Open();
    *from>>name>>factor;
    m_sfactor *= factor;
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (links[j]->Name()==name) {
	p_mapproc = p_partner = (Single_LOProcess*)links[j];
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



int AMEGIC::Single_LOProcess::InitAmplitude(Amegic_Model * model,Topology* top,
					    vector<Process_Base *> & links,
					    vector<Process_Base *> & errs)
{
  THROW(fatal_error,"Invalid function call");
}

int AMEGIC::Single_LOProcess::InitAmplitude(Amegic_Model * model,Topology* top,
					    vector<Process_Base *> & links,
					    vector<Process_Base *> & errs,int checkloopmap)
{
  m_type = 20;
  if (!model->p_model->CheckFlavours(m_nin,m_nout,&m_flavs.front())) return 0;
  model->p_model->GetCouplings(m_cpls);

  m_partonlist.clear();
  for (size_t i=0;i<m_nin;i++) if (m_flavs[i].Strong()) m_partonlist.push_back(i);
  for (size_t i=m_nin;i<m_nin+m_nout;i++) if (m_flavs[i].Strong()) m_partonlist.push_back(i);

  if (m_gen_str>1) {
    ATOOLS::MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename);
  }
  string newpath=rpa->gen.Variable("SHERPA_CPP_PATH");
  ATOOLS::MakeDir(newpath);

  if (CheckAlternatives(links,Name())) return 1;

  p_hel    = new Helicity(m_nin,m_nout,&m_flavs.front(),p_pl);

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
    string hstr=rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_libname;
    string hstr2=rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_name+".map";
    p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nin+m_nout,&m_flavs.front(),p_b,hstr,hstr2);  
  }
  else p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nin+m_nout,&m_flavs.front(),p_b);  
  p_BS->Setk0(s_gauge);
  p_shand  = new String_Handler(m_gen_str,p_BS,model->p_model->GetCouplings());
  int ntchanmin(m_ntchanmin);
  bool cvp(reader.GetValue<int>("AMEGIC_CUT_MASSIVE_VECTOR_PROPAGATORS",1));
  p_ampl   = new Amplitude_Handler(m_nin+m_nout,&m_flavs.front(),p_b,p_pinfo,model,top,m_maxcpl,m_mincpl,ntchanmin,
                                   &m_cpls,p_BS,p_shand,m_print_graphs,!directload,cvp);
  m_ntchanmin=ntchanmin;
  if (p_ampl->GetGraphNumber()==0) {
    msg_Tracking()<<"AMEGIC::Single_LOProcess::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return 0;
  }
  map<string,Complex> cplmap;
  for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
    cplmap.clear();
    if (checkloopmap && !NaiveMapping(links[j])) continue;
    if (checkloopmap==2 || m_pinfo.m_special.find("MapOff")!=std::string::npos) continue;
    if (FlavCompare(links[j]) && p_ampl->CompareAmplitudes(links[j]->GetAmplitudeHandler(),m_sfactor,cplmap)) {
      if (p_hel->Compare(links[j]->GetHelicity(),m_nin+m_nout)) {
	m_sfactor = sqr(m_sfactor);
	msg_Tracking()<<"AMEGIC::Single_LOProcess::InitAmplitude : Found compatible process for "<<Name()<<" : "<<links[j]->Name()<<endl;
	  
	bool found(true);
	if (!FoundMappingFile(m_libname,m_pslibname)) {
	  string mlname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+links[j]->Name();
	  string mnname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+Name();
	  if (FileExists(mlname+string(".map"))) { 
	    if (m_sfactor==1.) My_In_File::CopyInDB(mlname+".map",mnname+".map");
	    else {
	      UpdateMappingFile(mlname,cplmap);
	      CreateMappingFile((Single_LOProcess*)links[j]);
	    }
	    My_In_File::CopyInDB(mlname+".col",mnname+".col");
	    for (size_t i=0;i<m_nin+m_nout-1;i++) if (m_flavs[i].Strong()) {
	      for (size_t j=i+1;j<m_nin+m_nout;j++) if (m_flavs[j].Strong()) {
		string sij=string("__S")+ToString(i)+string("_")+ToString(j);
		My_In_File::CopyInDB(mlname+sij+".col",mnname+sij+".col");
	      }
	    }
	  }
	  found=false;
	}

	p_mapproc = p_partner = (Single_LOProcess*)links[j];
	for (std::map<string,Flavour>::const_iterator fit=p_ampl->GetFlavourmap().begin();
	     fit!=p_ampl->GetFlavourmap().end();fit++) AddtoFlavmap(fit->first,fit->second);
	InitFlavmap(p_partner);
	FillCombinations();
	if (!found) WriteAlternativeName(p_partner->Name());
	m_iresult = p_partner->Result()*m_sfactor;

	Minimize();
	return 1;
      }
    }
  }
  if (directload) {
    p_ampl->CompleteLibAmplitudes(m_nin+m_nout,m_ptypename+string("/")+m_name,
				  m_ptypename+string("/")+m_libname,127,127,&m_flavs.front());    
    if (!p_shand->SearchValues(m_gen_str,m_libname,p_BS)) return 0;
    if (!TestLib()) return 0;
    if (p_partner==this) links.push_back(this);
    FillCombinations();
    msg_Info()<<".";
    Minimize();
    return 1;
  }

  p_ampl->CompleteAmplitudes(m_nin+m_nout,&m_flavs.front(),p_b,&m_pol,
			     top,p_BS,m_ptypename+string("/")+m_name,127,127);
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
    if (p_partner==this) links.push_back(this);
    if (CheckLibraries()) return 1;
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (links[j]->NewLibs()) {
	if (CheckStrings((Single_LOProcess*)links[j])) return 1;	
      }      
    }
    if (p_partner!=this) links.push_back(this);
    
    if (m_gen_str<2) return 1;
    if (p_partner!=this) {
      msg_Tracking()<<"AMEGIC::Single_LOProcess::InitAmplitude : "<<std::endl
		    <<"   Strings of process "<<m_name<<" and partner "
		    <<p_partner->Name()<<" did not fit."<<std::endl
		    <<"   Have to write new library."<<std::endl;
    }
    WriteLibrary();
    if (p_partner==this && Result()>0.) SetUpIntegrator();
    return 1;
  case -3: return -1;
  default :
    msg_Error()<<"ERROR in AMEGIC::Single_LOProcess::InitAmplitude : "<<std::endl
	       <<"   Failed for "<<m_name<<" with result "<<result<<"."<<endl;
    errs.push_back(this);
    return 0;
  }

  return 1;
}

int Single_LOProcess::InitAmplitude(Amegic_Model * model,Topology* top,
				    vector<Process_Base *> & links,
				    vector<Process_Base *> & errs,
				    std::vector<ATOOLS::Vec4D>* epol,std::vector<double> * pfactors)
{
  m_type = 10;
  model->p_model->GetCouplings(m_cpls);
  
  if (m_gen_str>1) {
    ATOOLS::MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename);
  }
  string newpath=rpa->gen.Variable("SHERPA_CPP_PATH");
  ATOOLS::MakeDir(newpath);

  m_name+= "__S"+ToString((int)m_emit)+"_"+ToString((int)m_spect);
  if (m_flavs[m_emit].IsGluon()) {
    p_pl[m_emit]=Pol_Info();
    p_pl[m_emit].Init(2);
    p_pl[m_emit].pol_type = 'e'; 
    p_pl[m_emit].type[0] = 90;
    p_pl[m_emit].type[1] = 91;
    p_pl[m_emit].factor[0] = 1.;
    p_pl[m_emit].factor[1] = 1.;
  }
  m_epol.resize(epol->size());
  for (size_t i=0;i<epol->size();i++) m_epol[i]=(*epol)[i];

  p_hel    = new Helicity(m_nin,m_nout,&m_flavs.front(),p_pl);

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
  if (m_libname=="0") {
    return 0;
  }
  if (directload) {
    string hstr=rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_libname;
    string hstr2=rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_name+".map";
    p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nin+m_nout,&m_flavs.front(),p_b,hstr,hstr2);
  }
  else p_BS     = new Basic_Sfuncs(m_nin+m_nout,m_nin+m_nout,&m_flavs.front(),p_b);  
  p_BS->Setk0(s_gauge);
  p_BS->SetEPol(&m_epol); 
  p_shand  = new String_Handler(m_gen_str,p_BS,model->p_model->GetCouplings());

 
  int ntchanmin(m_ntchanmin);
  bool cvp(reader.GetValue<int>("AMEGIC_CUT_MASSIVE_VECTOR_PROPAGATORS",1));
  p_ampl   = new Amplitude_Handler(m_nin+m_nout,&m_flavs.front(),p_b,p_pinfo,model,top,m_maxcpl,m_mincpl,ntchanmin,
                                   &m_cpls,p_BS,p_shand,m_print_graphs,!directload,cvp);
  m_ntchanmin=ntchanmin;
  if (p_ampl->GetGraphNumber()==0) {
    msg_Tracking()<<"Single_LOProcess::InitAmplitude : No diagrams for "<<m_name<<"."<<endl;
    return 0;
  }

  map<string,Complex> cplmap;
  for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
    cplmap.clear();
    if (m_pinfo.m_special.find("MapOff")!=std::string::npos) continue;
    if (FlavCompare(links[j]) && p_ampl->CompareAmplitudes(links[j]->GetAmplitudeHandler(),m_sfactor,cplmap)) {
      if (p_hel->Compare(links[j]->GetHelicity(),m_nin+m_nout)) {
	Single_LOProcess *pp=dynamic_cast<Single_LOProcess*>(links[j]);
	if (m_emit!=pp->m_emit || m_spect!=pp->m_spect ||
	    p_sub->m_ijt!=pp->p_sub->m_ijt || p_sub->m_kt!=pp->p_sub->m_kt ||
	    p_sub->m_i!=pp->p_sub->m_i || p_sub->m_j!=pp->p_sub->m_j || p_sub->m_k!=pp->p_sub->m_k) continue;
	m_sfactor = sqr(m_sfactor);
	msg_Tracking()<<"AMEGIC::Single_LOProcess::InitAmplitude : Found compatible process for "<<Name()<<" : "<<links[j]->Name()<<endl;
	  
	if (!FoundMappingFile(m_libname,m_pslibname)) {
	  string mlname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+links[j]->Name();
	  string mnname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+Name();
	  if (FileExists(mlname+string(".map"))) { 
	    if (m_sfactor==1.) My_In_File::CopyInDB(mlname+".map",mnname+".map");
	    else {
	      UpdateMappingFile(mlname,cplmap);
	      CreateMappingFile((Single_LOProcess*)links[j]);
	    }
	    My_In_File::CopyInDB(mlname+".col",mnname+".col");
	    for (size_t i=0;i<m_nin+m_nout-1;i++) if (m_flavs[i].Strong()) {
	      for (size_t j=i+1;j<m_nin+m_nout;j++) if (m_flavs[j].Strong()) {
		string sij=string("__S")+ToString(i)+string("_")+ToString(j);
		My_In_File::CopyInDB(mlname+sij+".col",mnname+sij+".col");
	      }
	    }
	  }
	  WriteAlternativeName(p_partner->Name());
	}

	p_mapproc = p_partner = dynamic_cast<Single_LOProcess*>(links[j]);
	for (std::map<string,Flavour>::const_iterator fit=p_ampl->GetFlavourmap().begin();
	     fit!=p_ampl->GetFlavourmap().end();fit++) AddtoFlavmap(fit->first,fit->second);
	InitFlavmap(p_partner);
	FillCombinations();
	m_iresult = p_partner->Result()*m_sfactor;

	Minimize();
	return 1;
      }
    }
  }
  if (directload) {
    p_ampl->CompleteLibAmplitudes(m_nin+m_nout,m_ptypename+string("/")+m_name,
				  m_ptypename+string("/")+m_libname,
				  m_emit,m_spect,&m_flavs.front());
    if (!p_shand->SearchValues(m_gen_str,m_libname,p_BS)) return 1;
    if (!TestLib(pfactors)) return 0;
    if (p_partner==this) links.push_back(this);
    FillCombinations();
    Minimize();
    return 1;
  }

  p_ampl->CompleteAmplitudes(m_nin+m_nout,&m_flavs.front(),p_b,&m_pol,
			     top,p_BS,m_ptypename+string("/")+m_name,
			     m_emit,m_spect);
  m_pol.Add_Extern_Polarisations(p_BS,&m_flavs.front(),p_hel);
  p_BS->Initialize();
  FillCombinations();

  int tr=Tests(pfactors);
  switch (tr) {
  case 2 : 
    if (p_partner==this) links.push_back(this);
    return 1;
  case 1 :
  case 100 :
    if (Result()==0.) {
      CreateMappingFile(this);
      return 0;
    }
    if (p_partner==this) links.push_back(this);
    
    if (CheckLibraries(pfactors)) return 1;
    for (size_t j=0;j<links.size();j++) if (Type()==links[j]->Type()) {
      if (links[j]->NewLibs()) {
	if (CheckStrings((Single_LOProcess*)links[j],pfactors)) return 1;	
      }      
    }
    if (p_partner!=this) links.push_back(this);
    
    if (m_gen_str<2) return 1;
    if (p_partner!=this) {
      msg_Tracking()<<"Single_LOProcess::InitAmplitude : "<<std::endl
		    <<"   Strings of process "<<m_name<<" and partner "
		    <<p_partner->Name()<<" did not fit."<<std::endl
		    <<"   Have to write new library."<<std::endl;
    }
    if (tr==1) WriteLibrary();
    return 1;
  case -3: return -3;
  default :
    msg_Error()<<"ERROR in Single_LOProcess::InitAmplitude : "<<std::endl
	       <<"   Failed for "<<m_name<<" with result "<<tr<<"."<<endl;
//     errs.push_back(this);
    return -1;
  }
  return 1;
}



int Single_LOProcess::Tests(std::vector<double> * pfactors) {

  int number      = 1;
  int gauge_test  = 1;
  int string_test = 1;

  /* ---------------------------------------------------
     
     The reference result for momenta moms

     --------------------------------------------------- */

  string testname = string("");
  int fmfbnl=0;
  if (FoundMappingFile(testname,m_pslibname)) {
    if (testname != string("")) {
      if (FoundLib(testname)) gauge_test = string_test = 0;
      else fmfbnl=99;
    }
  }
  
  if (gauge_test) p_shand->Initialize(p_ampl->GetRealGraphNumber(),p_hel->MaxHel());

  p_ampl->SetStringOff();

  double M2 = 0.;
  double helvalue;
  std::vector<ATOOLS::Vec4D> epol;
  for (size_t i=0;i<m_epol.size();i++) epol.push_back(m_epol[i]);

  if (gauge_test) {
    m_pol.Set_Gauge_Vectors(m_nin+m_nout,p_testmoms,Vec4D(sqrt(3.),1.,1.,-1.));
    if (m_epol.size()>0)
      for (size_t i=0;i<m_epol.size();i++) m_epol[i]+=p_testmoms[m_emit]/p_testmoms[m_emit][0];
    p_BS->Setk0(0);
    p_BS->CalcEtaMu(p_testmoms);  
    if (m_epol.size()==0) p_BS->InitGaugeTest(.9);

    msg_Info()<<"Single_LOProcess::Tests for "<<m_name<<std::endl
	      <<"   Prepare gauge test and init helicity amplitudes. This may take some time."
	      <<std::endl;
    for (size_t i=0;i<p_hel->MaxHel();i++) { 
      if (p_hel->On(i)) {
	helvalue = p_ampl->Differential(i,(*p_hel)[i])*p_hel->PolarizationFactor(i);
	if (pfactors) helvalue*=(double)(*pfactors)[p_hel->GetEPol(i)-90]; 
	M2      +=  helvalue;
      } 
    }
    M2     *= sqr(m_pol.Massless_Norm(m_nin+m_nout,&m_flavs.front(),p_BS));
    m_iresult  = M2;
  }
  for (size_t i=0;i<m_epol.size();i++) m_epol[i]=epol[i];

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
     tested, it is assumed that is doesn�t contribute at all and is switched off.      */
  for (size_t i=0;i<p_hel->MaxHel();i++) { 
    if (p_hel->On(i)) {
      M_doub[i]  = p_ampl->Differential(i,(*p_hel)[i])*p_hel->PolarizationFactor(i);  
      if (pfactors) M2g+= M_doub[i]*(double)(*pfactors)[p_hel->GetEPol(i)-90];
      else M2g+= M_doub[i];
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
  msg_Tracking()<<"Single_LOProcess::Tests for "<<m_name<<std::endl
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

      p_BS->CalcEtaMu((ATOOLS::Vec4D*)p_testmoms);
      p_hel->InitializeSpinorTransformation(p_BS);
      p_shand->Calculate();
      
      for (size_t i=0;i<p_hel->MaxHel();i++) {
	if (p_hel->On(i)) {
	  double help = p_ampl->Differential(i)*p_hel->Multiplicity(i)*p_hel->PolarizationFactor(i);
	  if (pfactors) M2 += help*(*pfactors)[p_hel->GetEPol(i)-90];
	  else M2 += help;
	}
      } 

      p_hel->AllowTransformation();
      gauge_test = string_test = 0;
    }
    else {
      string searchfilename = rpa->gen.Variable("SHERPA_CPP_PATH")+string("/Process/Amegic/")+m_ptypename+string("/")+testname+string("/V.H");
      if (FileExists(searchfilename)) {
      	msg_Error()<<"ERROR in Single_LOProcess::Tests()"<<std::endl
		   <<"   No compiled & linked library found for process "<<testname<<std::endl
		   <<"   but files already written out !"<<std::endl
		   <<om::bold<<"   Interrupt run and execute \"makelibs\" in '"
		   <<rpa->gen.Variable("SHERPA_CPP_PATH")<<"'."
		   <<om::reset<<std::endl;
	int stat;
	stat=system((string("cp ")+rpa->gen.Variable("SHERPA_SHARE_PATH")+
		string("/makelibs ")+rpa->gen.Variable("SHERPA_CPP_PATH")).c_str());
	THROW(normal_exit,"Failed to load library.");
      }
      else {
      	msg_Error()<<"ERROR in Single_LOProcess::Tests()"<<std::endl
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
	msg_Out()<<"XX: Library cross check: "
		 <<M2<<" vs. "<<M2g<<"  difference:"<<abs(M2/M2g-1.)*100.<<"%"<<endl
		 <<"   Mapping file(1) : "<<abs(M2)<<endl
		 <<"   Original    (2) : "<<abs(M2g)<<endl
		 <<"   Cross check (T) : "<<abs(M2/M2g-1.)*100.<<"%"<<endl;
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
    
    if (m_emit==m_spect) {
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
    else {
      for (size_t i=0;i<p_hel->MaxHel();i++) {
	if (p_hel->On(i)) {
	  for (size_t j=i+1;j<p_hel->MaxHel();j++) {
	    if (p_hel->On(j)) {
#ifdef FuckUp_Helicity_Mapping
	      if (ATOOLS::IsEqual(M_doub[i]*(*pfactors)[p_hel->GetEPol(i)-90],M_doub[j]*(*pfactors)[p_hel->GetEPol(j)-90])) {
		p_hel->SwitchOff(j);
		p_hel->SetPartner(i,j);
		p_hel->IncMultiplicity(i);
	      }
#endif
	    }
	  }
	}
      }
      for (size_t i=0;i<p_hel->MaxHel();i++) {
	if (p_hel->On(i)) {
	  for (size_t j=i+1;j<p_hel->MaxHel();j++) {
	    if (p_hel->On(j)) {
#ifdef FuckUp_Helicity_Mapping
	      if (ATOOLS::IsEqual(M_doub[i],M_doub[j]) && p_hel->Multiplicity(i)==p_hel->Multiplicity(j)) {
		p_hel->SwitchOff(j);
		p_hel->SetPartner(i,j);
		p_hel->IncMultiplicity(i,p_hel->GetEPol(j)*1024);
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
	  double multi = (p_hel->Multiplicity(i)%1024);
	  if (pfactors) {
	    if (p_hel->Multiplicity(i)<1024) multi*=(*pfactors)[p_hel->GetEPol(i)-90];
	    else multi*=(*pfactors)[p_hel->GetEPol(i)-90]+(*pfactors)[p_hel->Multiplicity(i)/1024-90];
	  }
	  M2S += p_ampl->Differential(i)*p_hel->PolarizationFactor(i)*multi;
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
      return 1+fmfbnl;
    }
    return 1+fmfbnl;
  }      

  delete[] M_doub;

  return 0;
}

int Single_LOProcess::TestLib(std::vector<double> * pfactors)
{
  double M2(0.);
  double * M_doub = new double[p_hel->MaxHel()];
  p_BS->CalcEtaMu((ATOOLS::Vec4D*)p_testmoms);
  p_hel->InitializeSpinorTransformation(p_BS);
  p_shand->Calculate();
  
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    M_doub[i] = p_ampl->Differential(i)*p_hel->Multiplicity(i)*p_hel->PolarizationFactor(i);
    if (IsNan(M_doub[i])) {
      msg_Error()<<METHOD<<"("<<m_name<<"): Helicity "<<i<<" yields "<<M_doub[i]<<". Continuing."<<std::endl;
      continue;
    }
    if (pfactors) M2 += M_doub[i]*(*pfactors)[p_hel->GetEPol(i)-90];
    else M2 += M_doub[i];
  } 
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (M_doub[i]==0.) {
#ifdef FuckUp_Helicity_Mapping
      p_hel->SwitchOff(i);
#endif
    }
  }

  if (!(rpa->gen.SoftSC()||rpa->gen.HardSC())) {
    if (m_emit==m_spect) {
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
    else {
      for (size_t i=0;i<p_hel->MaxHel();i++) {
	if (p_hel->On(i)) {
	  for (size_t j=i+1;j<p_hel->MaxHel();j++) {
	    if (p_hel->On(j)) {
#ifdef FuckUp_Helicity_Mapping
	      if (ATOOLS::IsEqual(M_doub[i]*(*pfactors)[p_hel->GetEPol(i)-90],M_doub[j]*(*pfactors)[p_hel->GetEPol(j)-90])) {
		p_hel->SwitchOff(j);
		p_hel->SetPartner(i,j);
		p_hel->IncMultiplicity(i);
	      }
#endif
	    }
	  }
	}
      }
      for (size_t i=0;i<p_hel->MaxHel();i++) {
	if (p_hel->On(i)) {
	  for (size_t j=i+1;j<p_hel->MaxHel();j++) {
	    if (p_hel->On(j)) {
#ifdef FuckUp_Helicity_Mapping
	      if (ATOOLS::IsEqual(M_doub[i],M_doub[j]) && p_hel->Multiplicity(i)==p_hel->Multiplicity(j)) {
		p_hel->SwitchOff(j);
		p_hel->SetPartner(i,j);
		p_hel->IncMultiplicity(i,p_hel->GetEPol(j)*1024);
	      }
#endif
	    }
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

int Single_LOProcess::CheckLibraries(std::vector<double> * pfactors) {
  if (m_gen_str==0) return 1;
  if (p_shand->IsLibrary()) return 1;

  msg_Info()<<"Single_LOProcess::CheckLibraries : Looking for a suitable library. This may take some time."<<std::endl;
  String_Handler * shand1;
  shand1      = new String_Handler(p_shand->Get_Generator());
  
  string testname;
  double M2s, helvalue;

  for (;;) {
    testname  = CreateLibName();
    if (shand1->SearchValues(m_gen_str,testname,p_BS)) {
      shand1->Calculate();
      
      M2s = 0.;
      for (size_t i=0;i<p_hel->MaxHel();i++) {
	double multi = (p_hel->Multiplicity(i)%1024);
	if (pfactors) {
	  if (p_hel->Multiplicity(i)<1024) multi*=(*pfactors)[p_hel->GetEPol(i)-90];
	  else multi*=(*pfactors)[p_hel->GetEPol(i)-90]+(*pfactors)[p_hel->Multiplicity(i)/1024-90];
	}
	helvalue = p_ampl->Differential(shand1,i) * p_hel->PolarizationFactor(i);
	M2s     += helvalue * multi;
      } 
      M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,&m_flavs.front(),p_BS));
      if (ATOOLS::IsEqual(M2s,Result())) {
	m_libname = testname;
	m_pslibname = testname;
	if (shand1) { delete shand1; shand1 = 0; }
	//Clean p_shand!!!!
// 	Minimize();
	CreateMappingFile(this);
	return 1;
      }
    } 
    else break;
  }
  if (shand1) { delete shand1; shand1 = 0; }
  return 0;
}

int Single_LOProcess::CheckStrings(Single_LOProcess* tproc,std::vector<double> * pfactors)
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
    double multi = (p_hel->Multiplicity(i)%1024);
    if (pfactors) {
      if (p_hel->Multiplicity(i)<1024) multi*=(*pfactors)[p_hel->GetEPol(i)-90];
      else multi*=(*pfactors)[p_hel->GetEPol(i)-90]+(*pfactors)[p_hel->Multiplicity(i)/1024-90];
    }
    helvalue = p_ampl->Differential(shand1,i) * p_hel->PolarizationFactor(i);
    M2s     += helvalue * multi;
  }
  M2s *= sqr(m_pol.Massless_Norm(m_nin+m_nout,&m_flavs.front(),p_BS));
  (shand1->Get_Generator())->ReStore();
  delete shand1;

  if (ATOOLS::IsEqual(M2s,Result())) {
    m_libname = tproc->LibName();
    m_pslibname = tproc->PSLibName();
//     Minimize();
    CreateMappingFile(this);
    return 1;
  }
  return 0;
}
  
void Single_LOProcess::WriteLibrary() 
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

void Single_LOProcess::CreateMappingFile(Single_LOProcess* partner) {
  if (m_gen_str<2) return;
  std::string outname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_name+".map";
  if (FileExists(outname)) {
    string MEname,PSname;
    FoundMappingFile(MEname,PSname);
    if (MEname != m_libname || PSname != m_pslibname) {
      msg_Error()<<"ERROR in Single_LOProcess::CreateMappingFile() :"<<std::endl
		 <<"   Files do not coincide. Maybe changed input data ? Abort the run."<<std::endl
		 <<MEname<<" v "<<m_libname<<" || "<<PSname<<" v "<<m_pslibname<<endl;
       abort();
    }
    return;
  }

  My_Out_File to(outname);
  to.Open();
  if (Result()!=0.) {
    *to<<"ME: "<<m_libname<<endl
      <<"PS: "<<m_pslibname<<endl;
    p_shand->Get_Generator()->WriteCouplings(*to);
  }
  else {
    *to<<"ME: 0"<<endl
      <<"PS: 0"<<endl;
  }
  to.Close();
}

bool Single_LOProcess::FoundMappingFile(std::string & MEname, std::string & PSname) {
  
  std::string buf;
  int pos;
  std::string outname = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_ptypename+"/"+m_name+".map";
  if (FileExists(outname)) {
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

bool Single_LOProcess::FoundLib(std::string& pID)
{
  std::string libname=ATOOLS::rpa->gen.Variable("SHERPA_LIB_PATH")+
    std::string("/libProc_P")+pID.substr(1)+std::string(LIB_SUFFIX);
  if (FileExists(libname)) return 1;
  return 0;
}

void AMEGIC::Single_LOProcess::UpdateMappingFile(std::string name, map<string,Complex> & cmap) 
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

bool AMEGIC::Single_LOProcess::CompareTestMoms(const ATOOLS::Vec4D* p)
{
  for (size_t i=0;i<m_nin+m_nout-1;i++) if (!(p[i]==p_testmoms[i])) return 0;
  return 1;
}


int AMEGIC::Single_LOProcess::PerformTests()
{
  return 1;
}


/*------------------------------------------------------------------------------
  
  Phase space initialization
  
  ------------------------------------------------------------------------------*/


bool Single_LOProcess::SetUpIntegrator() 
{  
  return 0;
}

/*------------------------------------------------------------------------------
  
  Process management
  
  ------------------------------------------------------------------------------*/
void Single_LOProcess::Minimize()
{
  if (p_partner==this) return;
  if (p_hel)      {delete p_hel; p_hel=0;}
  if (p_BS)       {delete p_BS;   p_BS=0;}
  if (p_shand)    {delete p_shand;p_shand=0;}
  if (p_ampl)     {delete p_ampl; p_ampl=0;}

  m_maxcpl    = p_partner->MaxOrders();
  m_mincpl    = p_partner->MinOrders();
  m_ntchanmin = p_partner->NTchanMin();
}


/*------------------------------------------------------------------------------

  Calculating total cross sections
  
  ------------------------------------------------------------------------------*/

double Single_LOProcess::Partonic(const ATOOLS::Vec4D_Vector& _moms,const int mode) { return 0.; }

double Single_LOProcess::operator()(const ATOOLS::Vec4D_Vector &labmom,const ATOOLS::Vec4D *mom,
				    std::vector<double> * pfactors,std::vector<ATOOLS::Vec4D>* epol,const int mode)
{
  if (p_partner!=this) {
    return m_lastxs = p_partner->operator()(labmom,mom,pfactors,epol,mode)*m_sfactor;
  }

  double M2(0.);
  p_int->SetMomenta(labmom);
  p_scale->CalculateScale(labmom,m_cmode);
 
  for (size_t i=0;i<m_epol.size();i++) m_epol[i]=(*epol)[i];
  p_BS->CalcEtaMu((ATOOLS::Vec4D*)mom);
  p_hel->InitializeSpinorTransformation(p_BS);

  if (p_shand->Is_String()) {
    p_shand->Calculate();

    for (size_t i=0;i<p_hel->MaxHel();i++) {
      if (p_hel->On(i)) {
	double multi = (p_hel->Multiplicity(i)%1024);
	if (p_hel->Multiplicity(i)<1024) multi*=(*pfactors)[p_hel->GetEPol(i)-90];
	else multi*=(*pfactors)[p_hel->GetEPol(i)-90]+(*pfactors)[p_hel->Multiplicity(i)/1024-90];
	double mh=p_ampl->Differential(i);
	M2 += mh * multi * p_hel->PolarizationFactor(i);
      }
    }
  }
  m_lastxs = M2;
  return M2;
}

void Single_LOProcess::FillAmplitudes(vector<METOOLS::Spin_Amplitudes>& amps,
                                      std::vector<std::vector<Complex> >& cols)
{
  if (p_partner==this) p_ampl->FillAmplitudes(amps, cols, p_hel, 1.0);
  else p_partner->FillAmplitudes(amps, cols, sqrt(m_sfactor));
}

void Single_LOProcess::FillAmplitudes(vector<METOOLS::Spin_Amplitudes>& amps,
                                      std::vector<std::vector<Complex> >& cols,
                                      double sfactor)
{
  if (p_partner==this) p_ampl->FillAmplitudes(amps, cols, p_hel, sfactor);
  else p_partner->FillAmplitudes(amps, cols, sfactor*sqrt(m_sfactor));
}

double Single_LOProcess::Calc_M2ik(int ci, int ck) 
{
  double M2=0.;
  for (size_t i=0;i<p_hel->MaxHel();i++) {
    if (p_hel->On(i)) {
      M2 += p_ampl->Differential(i,ci,ck) * p_hel->Multiplicity(i) 
	* p_hel->PolarizationFactor(i);
      //	  M2     += helvalue;     
    }
  } 
  return M2;
}

void Single_LOProcess::Calc_AllXS(const ATOOLS::Vec4D_Vector &labmom,
				  const ATOOLS::Vec4D *mom,std::vector<std::vector<double> > &dsij,const int mode) 
{
  if (p_partner!=this) {
    p_partner->Calc_AllXS(labmom,mom,dsij,mode);
    dsij[0][0]*=m_sfactor;
    for (size_t i=0;i<m_partonlist.size();i++) {
      for (size_t k=i+1;k<m_partonlist.size();k++) {
	dsij[i][k] = dsij[k][i]*=m_sfactor;
      }
    }
    return;
  }
  p_int->SetMomenta(labmom);
  p_scale->CalculateScale(labmom,m_cmode);

  p_BS->CalcEtaMu((ATOOLS::Vec4D*)mom);
  p_hel->InitializeSpinorTransformation(p_BS);

  if (p_shand->Is_String()) {
    p_shand->Calculate();

    dsij[0][0] = Calc_M2ik(0,0);
    for (size_t i=0;i<m_partonlist.size();i++) {
      for (size_t k=i+1;k<m_partonlist.size();k++) {
	dsij[i][k] = dsij[k][i] = Calc_M2ik(m_partonlist[i],m_partonlist[k]);
      }
    }
  }
}


String_Handler *AMEGIC::Single_LOProcess::GetStringHandler()
{ 
  if (p_partner==this) return p_shand;
  return p_partner->GetStringHandler();
}

Amplitude_Handler *AMEGIC::Single_LOProcess::GetAmplitudeHandler()
{ 
  if (p_partner==this) return p_ampl;
  return p_partner->GetAmplitudeHandler();
}

Helicity *AMEGIC::Single_LOProcess::GetHelicity() 
{ 
  if (p_partner==this) return p_hel; 
  return p_partner->GetHelicity();
}    

int AMEGIC::Single_LOProcess::NumberOfDiagrams() { 
  if (p_partner==this) return p_ampl->GetGraphNumber(); 
  return p_partner->NumberOfDiagrams();
}

Point * AMEGIC::Single_LOProcess::Diagram(int i) { 
  if (p_partner==this) return p_ampl->GetPointlist(i); 
  return p_partner->Diagram(i);
} 

void Single_LOProcess::AddChannels(std::list<std::string>* tlist) 
{ }

void AMEGIC::Single_LOProcess::FillCombinations
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
    Flavour_Vector cf(m_cflavs[idc]);
    for (size_t i(0);i<cf.size();++i)
      if (cf[i]==p->fl) {
	in=true;
	break;
      }
    if (!in) {
      m_cflavs[idc].push_back(p->fl.Bar());
      m_cflavs[id].push_back(p->fl);
#ifdef DEBUG__Fill_Combinations
      msg_Debugging()<<"  flav "<<ID(idc)<<" / "
		     <<ID(id)<<" -> "<<p->fl<<"\n";
#endif
    }
  }
}

void AMEGIC::Single_LOProcess::FillCombinations()
{
#ifdef DEBUG__Fill_Combinations
  msg_Debugging()<<METHOD<<"(): '"<<m_name<<"' {\n";
#endif
  std::string fname=rpa->gen.Variable("SHERPA_CPP_PATH")
    +"/Process/Amegic/"+m_ptypename+"/"+m_name+".clu";
  if (FileExists(fname)) {
    My_In_File from(fname);
    if (!from.Open()) THROW(critical_error,"No clu for "+m_name);
    size_t size(0);
    *from>>size;
    for (size_t i(0);i<size;++i) {
      size_t ida(0), idb(0);
      *from>>ida>>idb;
      m_ccombs.insert(std::pair<size_t,size_t>(ida,idb));
    }
    *from>>size;
    for (size_t i(0);i<size;++i) {
      size_t id(0), fsize(0);
      *from>>id>>fsize;
      m_cflavs[id].resize(fsize);
      for (size_t j(0);j<fsize;++j) {
	long int fl(0);
	*from>>fl;
	m_cflavs[id][j]=Flavour(abs(fl),fl<0);
      }
    }
    std::string eof;
    *from>>eof;
    if (eof!="eof") THROW(critical_error,"Corrupted clu for "+m_name);
    from.Close();
    return;
  }
  size_t nd(NumberOfDiagrams());
  for (size_t i(0);i<nd;++i) {
    Point *p(Diagram(i));
    size_t id(1<<p->number);
    FillCombinations(p,id);
  }
#ifdef DEBUG__Fill_Combinations
  msg_Debugging()<<"  } -> "<<m_cflavs.size()
		 <<" flavours, "<<m_ccombs.size()
		 <<" combinations\n";
  msg_Debugging()<<"}\n";
#endif
  My_Out_File to(fname);
  to.Open();
  *to<<m_ccombs.size()<<"\n";
  for (Combination_Set::const_iterator
	 cit(m_ccombs.begin());cit!=m_ccombs.end();++cit)
    *to<<cit->first<<" "<<cit->second<<"\n";
  *to<<m_cflavs.size()<<"\n";
  for (CFlavVector_Map::const_iterator
	 cit(m_cflavs.begin());cit!=m_cflavs.end();++cit) {
    *to<<cit->first<<" "<<cit->second.size();
    for (Flavour_Vector::const_iterator fit(cit->second.begin());
	 fit!=cit->second.end();++fit) *to<<" "<<(long int)*fit;
    *to<<"\n";
  }
  *to<<"eof\n";
  to.Close();
}

bool AMEGIC::Single_LOProcess::Combinable
(const size_t &idi,const size_t &idj)
{
  Combination_Set::const_iterator 
    cit(m_ccombs.find(std::pair<size_t,size_t>(idi,idj)));
  return cit!=m_ccombs.end();
}

const Flavour_Vector &AMEGIC::Single_LOProcess::
CombinedFlavour(const size_t &idij)
{
  CFlavVector_Map::const_iterator fit(m_cflavs.find(idij));
  if (fit==m_cflavs.end()) THROW(fatal_error,"Invalid request");
  return fit->second;
}

