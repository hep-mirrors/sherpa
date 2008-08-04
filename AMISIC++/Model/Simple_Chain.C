#include "Simple_Chain.H"

#include "Phase_Space_Handler.H"
#include "Running_AlphaS.H"
#include "Single_XS.H"
#include "XS_Selector.H"
#include "Semihard_QCD.H"
#include "Particle.H"
#include "Random.H"
#include "My_Limits.H"
#include "Vegas.H"
#include "Run_Parameter.H"
#include "Shell_Tools.H"
#include "Remnant_Base.H"
#include "Data_Reader.H"
#include "MyStrStream.H"
#include "ISR_Handler.H"

#ifdef PROFILE__all
#define PROFILE__Simple_Chain
#endif
#ifdef PROFILE__Simple_Chain
#include "prof.hh"
#else
#define PROFILE_HERE
#endif

// #define DEBUG__Simple_Chain

#ifdef DEBUG__Simple_Chain
const std::string differentialfile=std::string("differential.dat");
const std::string integralfile=std::string("integral.dat");
const std::string normalizedfile=std::string("normalized.dat");
#endif

static double s_epsilon=1.0e-3, s_xsnd, s_xstot;

using namespace AMISIC;
using namespace ATOOLS;

Simple_Chain::Simple_Chain():
  MI_Base("Simple Chain",MI_Base::HardEvent,5,4,1),
  p_differential(NULL), p_total(NULL), m_norm(1.0), m_enhance(1.0), 
  m_maxreduction(1.0), m_xsextension("_xs.dat"), m_mcextension("MC"), 
  p_processes(NULL), p_fsrinterface(NULL), p_model(NULL),
  p_beam(NULL), p_isr(NULL), p_profile(NULL), m_nflavour(5), m_maxtrials(1000), 
  m_scalescheme(PHASIC::scl::gmeanpt+PHASIC::scl::div_by_2), m_kfactorscheme(1), 
  m_ecms(rpa.gen.Ecms()), m_external(false), m_regulate(false)
{
  Init();
}

Simple_Chain::Simple_Chain(MODEL::Model_Base *const model,
			   BEAM::Beam_Spectra_Handler *const beam,
			   PDF::ISR_Handler *const isr):
  MI_Base("Simple Chain",MI_Base::HardEvent,5,4,1),
  p_differential(NULL), p_total(NULL), m_norm(1.0), m_enhance(1.0),
  m_maxreduction(1.0), m_xsextension("_xs.dat"), m_mcextension("MC"), 
  p_processes(NULL), p_fsrinterface(NULL), p_model(model), 
  p_beam(beam), p_isr(isr), p_profile(NULL), m_nflavour(5), m_maxtrials(1000), 
  m_scalescheme(PHASIC::scl::gmeanpt+PHASIC::scl::div_by_2), m_kfactorscheme(1), 
  m_ecms(rpa.gen.Ecms()), m_external(true), m_regulate(false)
{
  Init();
  p_remnants[0]=p_isr->GetRemnant(0);
  p_remnants[1]=p_isr->GetRemnant(1);
}

void Simple_Chain::Init()
{
  SetInputFile("MI.dat");
  SetInputFile("XS.dat",1);
  SetInputFile("Run.dat",2);
  SetInputFile("Model.dat",3);
  SetOutputFile("SC.log");
  m_start[4]=m_start[0]=m_ecms/2;
  m_stop[4]=m_stop[0]=0.0;
  m_start[3]=m_start[2]=m_ecms/2;
  m_stop[3]=m_stop[2]=0.0;
  m_spkey.Assign("s' isr",4,0,PHASIC::Phase_Space_Handler::GetInfo());
  m_ykey.Assign("y isr",3,0,PHASIC::Phase_Space_Handler::GetInfo());
  m_isrspkey.Assign("s' isr mi",3,0,PHASIC::Phase_Space_Handler::GetInfo());
  m_isrykey.Assign("y isr mi",2,0,PHASIC::Phase_Space_Handler::GetInfo());
  p_remnants[1]=p_remnants[0]=NULL;
  p_gridcreator=NULL;
}

Simple_Chain::~Simple_Chain()
{
  CleanUp();
}

void Simple_Chain::CleanUp() 
{
  if (p_gridcreator!=NULL) {
    exh->RemoveTerminatorObject(this);
    delete p_gridcreator;
    p_gridcreator=NULL;
  }
  if (p_fsrinterface!=NULL) {
    for (size_t i=0;i<p_processes->Size();++i) {
      Semihard_QCD *group = dynamic_cast<Semihard_QCD*>((*p_processes)[i]);
      group->SetFSRMode(3);
      group->CreateFSRChannels();
    }
    delete p_fsrinterface;
    p_fsrinterface=NULL;
  }
  if (p_differential!=NULL) {
    delete p_differential;
    p_differential=NULL;
  }
  if (p_total!=NULL) {
    delete p_total;
    p_total=NULL;
  }
  while (m_differentials.size()>0) {
    delete m_differentials.begin()->second;
    m_differentials.erase(m_differentials.begin());
  }
  if (p_profile!=NULL) {
    delete p_profile;
    p_profile=NULL;
  }
}

bool Simple_Chain::GeneratePathName()
{
  std::string outputpath, help[2];
  MyStrStream converter;
  converter<<rpa.gen.Bunch(0);
  converter>>help[0];
  converter.clear();
  converter<<rpa.gen.Bunch(1);
  converter>>help[1];
  outputpath=std::string("MIG_")+help[0]+help[1]+
    std::string("_")+ToString(rpa.gen.Ecms());
  if (m_regulate) {
    outputpath+=std::string("_")+m_regulator[0];
    for (size_t i=0;i<m_regulation.size();++i) {
      outputpath+=std::string("_")+ToString(m_regulation[i]);
    }
  }
  if (p_isr->PDF(0)->Type()!=p_isr->PDF(1)->Type()) {
    outputpath+=std::string("_")+p_isr->PDF(0)->Type();
  }
  outputpath+=std::string("_")+p_isr->PDF(0)->Type()+std::string("_")+
    ToString(static_cast<MODEL::Running_AlphaS*>
		     (p_model->GetScalarFunction("alpha_S"))->Order())+
    std::string("_")+ToString(m_scalescheme)+
    std::string("_")+ToString(m_kfactorscheme)+std::string("/");
  SetOutputPath(OutputPath()+m_pathextra+outputpath);
  m_pathextra=outputpath;
  return true;
}

void Simple_Chain::OrderFlavours(ATOOLS::Flavour *flavs)
{
  if ((int)flavs[0].Kfcode()>(int)flavs[1].Kfcode()) 
    std::swap<Flavour>(flavs[0],flavs[1]);
  if ((int)flavs[2].Kfcode()>(int)flavs[3].Kfcode()) 
    std::swap<Flavour>(flavs[2],flavs[3]);
  if (flavs[0].IsAnti()) 
    std::swap<Flavour>(flavs[0],flavs[1]);
  if (flavs[2].IsAnti()) 
    std::swap<Flavour>(flavs[2],flavs[3]);
}

EXTRAXS::XS_Group *Simple_Chain::FindPDFGroup(const size_t nin,const size_t nout,
					      const ATOOLS::Flavour *flavours)
{
  if (p_remnants[0]==NULL || p_remnants[1]==NULL) return p_processes;
  for (size_t i=0;i<p_processes->Size();++i) {
    if (nin==2 && nout==(*p_processes)[i]->NOut()) {
      Flavour ref[2], test[2];
      for (size_t j=0;j<2;++j) {
	ref[j]=p_remnants[j]->ConstituentType((*p_processes)[i]->Flavours()[j]);
	test[j]=p_remnants[j]->ConstituentType(flavours[j]);
      }
      if (ref[0]==test[0] && ref[1]==test[1])
	return dynamic_cast<EXTRAXS::XS_Group*>((*p_processes)[i]);
    }
  }
  Flavour *copy = new Flavour[nin+nout];
  for (short unsigned int i=0;i<nin;++i) 
    copy[i]=p_remnants[i]->ConstituentType(flavours[i]);
  for (short unsigned int i=nin;i<nin+nout;++i) copy[i]=kf_jet;
  Semihard_QCD *newgroup = 
    new Semihard_QCD(p_beam,p_isr,p_processes->SelectorData(),
		     copy,m_scalescheme,m_kfactorscheme);
  newgroup->PSHandler(false)->SetError(m_error);
  newgroup->SetScaleScheme(m_scalescheme);
  newgroup->SetKFactorScheme(m_kfactorscheme);
  p_processes->Add(newgroup);
  p_processes->SetAtoms(1);
  delete [] copy;
  return newgroup;
}

bool Simple_Chain::AddProcess(EXTRAXS::XS_Group *const group,
			      const ATOOLS::Flavour *flavs)
{
  bool success=false;
  Flavour help[4];
  for (size_t i=0;i<(size_t)flavs[0].Size();++i) {
    for (size_t j=0;j<(size_t)flavs[1].Size();++j) {
      for (size_t k=0;k<(size_t)flavs[2].Size();++k) {
	for (size_t l=0;l<(size_t)flavs[3].Size();++l) {
	  help[0]=flavs[0][i];
	  help[1]=flavs[1][j];
	  help[2]=flavs[2][k];
	  help[3]=flavs[3][l];
	  OrderFlavours(help);
	  EXTRAXS::XS_Base *newxs = p_processes->XSSelector()->
	    GetXS(2,2,help,false,0,2,false);
	  if (newxs==NULL) continue;
	  EXTRAXS::XS_Base *testxs=NULL;
	  EXTRAXS::XS_Group *pdfgroup = FindPDFGroup(2,2,help);
          EXTRAXS::XS_Selector::FindInGroup(pdfgroup,testxs,2,2,help);
          if (testxs==NULL) {
	    if (m_regulate) newxs->AssignRegulator(m_regulator,m_regulation);
	    newxs->SetScaleScheme(m_scalescheme);
	    newxs->SetKFactorScheme(m_kfactorscheme);
	    pdfgroup->Add(newxs);
	    pdfgroup->CreateISRChannels();
	    m_processmap[newxs->Name()]=newxs;
	    success=true;
	    msg_Debugging()<<"Simple_Chain::AddProcess(..): "
			   <<"New process '"<<newxs->Name()<<"'.\n";
	  }
	}
      }
    }
  }
  return success;
}

bool Simple_Chain::ReadInData()
{
  PROFILE_HERE;
  Data_Reader *reader = new Data_Reader(" ",";","!","=");
  reader->AddWordSeparator("\t");
  reader->SetInterprete(true);
  reader->SetInputPath(InputPath());
  reader->SetInputFile(InputFile());
  int regulate=0;
  if (reader->ReadFromFile(regulate,"REGULATE_XS")) {
    m_regulate=regulate;
    if (!reader->ReadFromFile(m_regulator,"XS_REGULATOR")) 
      m_regulator="QCD_Trivial";
    if (!reader->VectorFromFile(m_regulation,"XS_REGULATION")) 
      m_regulation=std::vector<double>(1,.71);
    double exponent;
    if (reader->ReadFromFile(exponent,"RESCALE_EXPONENT")) {
      double scale;
      if (!reader->ReadFromFile(scale,"REFERENCE_SCALE")) scale=1960.0;
      m_regulation[0]*=pow(m_ecms/scale,exponent);
    }
  }
  int helpssc;
  reader->SetTags(PHASIC::Integrable_Base::ScaleTags());
  if (reader->ReadFromFile(helpssc,"MI_SCALE_SCHEME")) 
    m_scalescheme=(PHASIC::scl::scheme)helpssc;
  else m_scalescheme=PHASIC::scl::gmeanpt;
  if (!reader->ReadFromFile(m_kfactorscheme,"MI_K_FACTOR_SCHEME")) 
    m_kfactorscheme=1;
  if (!reader->ReadFromFile(m_nflavour,"N_FLAVOUR")) m_nflavour=5;
  if (!reader->ReadFromFile(m_error,"PS_ERROR")) m_error=1.e-2;
  if (!reader->ReadFromFile(m_pathextra,"PATH_EXTRA")) m_pathextra="";
  GeneratePathName();
  delete reader;
  return true;
}

bool Simple_Chain::CreateGrid()
{
  PROFILE_HERE;
  bool vegas=PHASIC::Vegas::OnExternal();
  PHASIC::Vegas::SetOnExternal(m_vegas);
  double min=Min(m_stop[0],m_stop[4]);
  p_isr->SetFixedSprimeMin(4.0*min*min);
  p_isr->SetFixedSprimeMax(4.0*m_start[0]*m_start[0]);
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader(" ",";","!","=");
  reader->AddWordSeparator("\t");
  reader->SetInputPath(InputPath());
  reader->SetInputFile(InputFile(2));
  if (!reader->ReadFromFile(m_selectorfile,"MI_SELECTOR_FILE")) 
    m_selectorfile="MICuts.dat";
  p_processes = new EXTRAXS::Simple_XS(InputPath(),InputFile(1),p_model);
  if (p_processes->Size()>0) p_processes->Clear();
  p_processes->InitializeProcesses(p_beam,p_isr,false);  
  delete p_processes->SelectorData();
  p_processes->SetSelectorData
    (new ATOOLS::Selector_Data(InputPath(),m_selectorfile));
  p_processes->SetScaleScheme(m_scalescheme);
  p_processes->SetKFactorScheme(m_kfactorscheme);
  reader->SetInputFile(InputFile());
  reader->AddIgnore("->");
  reader->AddIgnore("to");
  reader->AddIgnore("for");
  reader->RescanInFile();
  std::vector<std::vector<std::string> > temp;
  reader->MatrixFromFile(temp,"CREATE_GRID");
  bool found=false;
  for (unsigned int i=0;i<temp.size();++i) {
    if (temp[i].size()>3) {
      Flavour flavour[4];
      int current;
      bool success=true;
      for (unsigned int j=0;j<4;++j) {
	flavour[j]=Flavour(s_kftable.KFFromIDName(temp[i][j]));
	if (flavour[j].Kfcode()==kf_none) {
	  reader->SetString(temp[i][j]);
	  reader->ReadFromString(current);
	  flavour[j]=Flavour((kf_code)abs(current));
	  if (current<0) flavour[j]=flavour[j].Bar();
	  if (flavour[j].Kfcode()==kf_none) success=false;
	}
      }
      if (!success) continue;
      if (AddProcess(p_processes,flavour)) found=true;
    }
  }
  delete reader;
  if (!found) {
    msg_Error()<<"Simple_Chain::CreateGrid(): "
	       <<"Did not find any process in '"
	       <<InputFile()<<"'."<<std::endl;
    PHASIC::Vegas::SetOnExternal(vegas);
    return false;
  }
  p_gridcreator = new Grid_Creator(&m_differentials,p_processes);
  p_gridcreator->SetGridXMin(min);
  p_gridcreator->SetGridXMax(m_ecms/2.0);
  p_gridcreator->ReadInArguments(InputFile(),InputPath());
  p_gridcreator->SetXSExtension(m_xsextension);
  p_gridcreator->SetMCExtension(m_mcextension);
  p_gridcreator->SetOutputPath(OutputPath());
  if (msg_LevelIsTracking()) {
    msg_Out()<<"Simple_Chain::CreateGrid(..): Process group {\n";
    msg_Indentation(3);
    p_processes->Print();
  }
  msg_Tracking()<<"}"<<std::endl;
  if (!p_gridcreator->ReadInGrid()) {
    if (MakeDir(OutputPath())==0) {
      msg_Tracking()<<"Simple_Chain::CreateGrid(..): "
		    <<"Created output directory "
		    <<OutputPath()<<"."<<std::endl;
    }
    p_gridcreator->CreateGrid();
  }
  PHASIC::Vegas::SetOnExternal(vegas);
  exh->AddTerminatorObject(this);
  Reset();
  return true;
}

bool Simple_Chain::SetUpInterface()
{
  p_processes->Reset();
  Flavour flavour[4]={kf_jet,kf_jet,kf_jet,kf_jet};
  if (p_fsrinterface!=NULL) delete p_fsrinterface;
  p_fsrinterface = new FSR_Channel(2,2,flavour,
				   p_total->XAxis()->Variable()->Name());
  for (size_t i=0;i<p_processes->Size();++i) {
    Semihard_QCD *group = dynamic_cast<Semihard_QCD*>((*p_processes)[i]);
    group->InitIntegrators();
    group->CreateISRChannels();
    group->SetFSRInterface(p_fsrinterface);
    group->SetFSRMode(2);
    group->CreateFSRChannels();
  }
  p_xs=p_processes;
  return true;
}

bool Simple_Chain::CheckConsistency(EXTRAXS::XS_Group *const group,
				    Amisic_Histogram_Type *const grid,
				    const double min,const double max,
				    const double integral)
{  
  int helpi=0, criterion=grid->XAxis()->Variable()->SelectorID();
  std::vector<Flavour> flavours(1,(kf_jet));
  std::vector<std::pair<double,double> > bounds
    (1,std::pair<double,double>(min,max));
  group->SelectorData()->AddData(criterion,flavours,bounds,helpi);
  double emin=Min(m_stop[0],m_stop[4]);
  p_isr->SetFixedSprimeMin(4.0*emin*emin);
  p_isr->SetFixedSprimeMax(4.0*m_start[0]*m_start[0]);
  group->ResetSelector(group->SelectorData());
  int level=msg->Level();
  msg->SetLevel(0);
  double error=group->PSHandler(false)->Error();
  group->PSHandler(false)->SetError(m_error);
  group->CalculateTotalXSec("");
  group->PSHandler(false)->SetError(error);
  msg->SetLevel(level);
  double total=group->TotalXS();
  msg_Info()<<"Simple_Chain::CheckConsistency(): {\n"
	    <<"   \\sigma_{hard xs}   = "
	    <<(total*rpa.Picobarn()/1.e9)<<" mb\n"
	    <<"   \\sigma_{hard grid} = "
	    <<(integral*rpa.Picobarn()/1.e9)<<" mb\n"
	    <<"   relative error     = "
	    <<dabs((total-integral)/(total))*100.0
	    <<" %\n}"<<std::endl;
  if (dabs((total-integral)/total)>m_error) {
    msg_Error()<<"Simple_Chain::CheckConsistency(..): Warning.\n"
	       <<"   \\Delta_{rel}\\sigma / m_error = "
	       <<dabs((total-integral)/
		      (total*m_error))<<std::endl;
    if (dabs((total-integral)/total)>2.*m_error) 
      THROW(fatal_error,"Grid integral and \\sigma_{tot} do not coincide.");
  }
  return true;
}

void Simple_Chain::CalculateSigmaND()
{
  if(s_kftable.find(111)==s_kftable.end()) // if not initialized yet
    s_kftable[111]=new Particle_Info(111,0.134976,7.8486e-09,0,0,0,1,0,"pi","pi");
  double eps=0.0808, eta=-0.4525, X=21.70, Y=56.08, b=2.3;
  if (p_isr->Flav(0).IsAnti()^p_isr->Flav(1).IsAnti()) Y=98.39;
  double s=sqr(rpa.gen.Ecms());
  double mp=Flavour(kf_p_plus).Mass();
  double mpi=Flavour(kf_pi).Mass();
  double ap=0.25, s0=8.0, y0=log(s/(mp*mp));
  double M1res=2.0, M2res=2.0, cres=2.0;
  double M1min=mp+2.0*mpi, M2min=mp+2.0*mpi;
  double ymin=4.0*log(1.0+2.0*mpi/mp);
  double MmaxAX2=0.213*s;
  double Del0=3.2-9.0/log(s)+17.4/sqr(log(s));
  double MmaxXX2=s*(0.07-0.44/log(s)+1.36/sqr(log(s)));
  double BAX=-0.47+150.0/s;
  double BXX=-1.05+40.0/sqrt(s)+8000.0/(s*s);
  double JAX=0.5/ap*log((b+ap*log(s/(M2min*M2min)))/(b+ap*log(s/MmaxAX2)));
  JAX+=0.5*cres/(b+ap*log(s/(M2res*M2min))+BAX)*log(1.0+M2res*M2res/
						    (M2min*M2min));
  double s1, s2, s3, JXX=0.5/ap*((y0-ymin)*(log((y0-ymin)/Del0)-1.0)+Del0);
  JXX+=s1=0.5*cres/ap*log(log(s*s0/(M1min*M1min*M2res*M2min))/
		       log(s*s0/(MmaxXX2*M2res*M2min)))*
    log(1.0+M2res*M2res/(M2min*M2min));
  JXX+=s2=0.5*cres/ap*log(log(s*s0/(M2min*M2min*M1res*M1min))/
		       log(s*s0/(MmaxXX2*M1res*M1min)))*
    log(1.0+M1res*M1res/(M1min*M1min));
  JXX+=s3=cres*cres/(2.0*ap*log(s*s0/(M1res*M2res*M1min*M2min))+BXX)*
    log(1.0+M1res*M1res/(M1min*M1min))*log(1.0+M2res*M2res/(M2min*M2min));
  double xstot=X*pow(s,eps)+Y*pow(s,eta);
  double xsel=0.0511*xstot*xstot/(4*(b+pow(s,eps))-4.2);
  double xssd=0.0336*X*sqrt(X)*JAX;
  double xsdd=0.0084*X*JXX;
  s_xstot=xstot;
  s_xsnd=xstot-xsel-2.0*xssd-xsdd;
  msg_Tracking()<<"Simple_Chain::CalculateSigmaND(): Results are {\n"
		<<"   \\sigma_{tot} = "<<xstot<<" mb\n"
		<<"   \\sigma_{el}  = "<<xsel<<" mb\n"
		<<"   \\sigma_{sd}  = "<<2.0*xssd<<" mb\n"
		<<"   \\sigma_{dd}  = "<<xsdd<<" mb\n"
		<<"   \\sigma_{nd}  = "<<(xstot-xsel-2.0*xssd-xsdd)
		<<" mb.\n}"<<std::endl;
  SetNorm((xstot-xsel-2.0*xssd-xsdd)*1.0e9/rpa.Picobarn());
}

bool Simple_Chain::CalculateTotal()
{
  PROFILE_HERE;
  if (m_differentials.size()==0) return false;
  Amisic_Histogram_Type *ref=m_differentials.begin()->second;
  p_differential = new Amisic_Histogram_Type();
  Axis<double> *xaxis=p_differential->XAxis(), *refx=ref->XAxis();
  Axis<double> *yaxis=p_differential->YAxis(), *refy=ref->YAxis();
  xaxis->SetVariable(refx->Variable()->Name());
  yaxis->SetVariable(refy->Variable()->Name());
  xaxis->SetScaling(refx->Scaling()->Name());
  yaxis->SetScaling(refy->Scaling()->Name());
  p_differential->Initialize(ref->XMin(),ref->XMax(),ref->NBins());
  for (Process_Map::iterator pit=m_processmap.begin();
       pit!=m_processmap.end();++pit) {
    Amisic_Histogram_Map::iterator diffit;
    for (diffit=m_differentials.begin();
	 diffit!=m_differentials.end();++diffit) {
      if (m_processmap[diffit->first]==pit->second) {
	for (size_t i=1;i<diffit->second->NBins()-1;++i) {
	  p_differential->Add(diffit->second->BinXMean(i),
			      diffit->second->BinContent(i));
	}
      }
    }
  }
  p_differential->SetFinished(true);
#ifdef DEBUG__Simple_Chain
  std::vector<std::string> comments(1,"  Differential XS   "); 
  p_differential->WriteOut(OutputPath()+differentialfile,
			   "[x,w,w2,max,n] = ",comments);
#endif
  SetStart(p_differential->XMax(),0);
  SetStop(Max(p_differential->XMin(),m_stop[0]),0);
  p_total = p_differential->GetIntegral(true);
  xaxis=p_total->XAxis();
  yaxis=p_total->YAxis();
  xaxis->SetVariable(refx->Variable()->Name());
  yaxis->SetVariable(refy->Variable()->Name());
  xaxis->SetScaling(refx->Scaling()->Name());
  yaxis->SetScaling(refy->Scaling()->Name());
#ifdef DEBUG__Simple_Chain
  comments[0]="     Total XS      "; 
  p_total->WriteOut(OutputPath()+integralfile,"[x,w,w2,max,n] = ",comments);
#endif
  m_sigmahard=(*p_total)(m_stop[4]);
  msg_Info()<<"Simple_Chain::CalculateTotal(): Result is {\n   "
	    <<"\\sigma_{hard} = "<<(m_sigmahard*rpa.Picobarn()/1.e9)
	    <<" mb\n   at "<<xaxis->Variable()->Name()<<"_{min} = "
	    <<m_stop[4]<<" GeV\n}"<<std::endl;
  CalculateSigmaND();
  if (m_sigmahard<m_norm) {
    msg_Error()<<"Simple_Chain::CalculateTotal(): "<<om::red
	       <<"\\sigma_{hard} = "
	       <<(m_sigmahard*rpa.Picobarn()/1.e9)
	       <<" mb < \\sigma_{nd} = "
	       <<(m_norm*rpa.Picobarn()/1.e9)
	       <<" mb !"<<om::reset<<std::endl;
  }
  if (m_check) {
    Flavour help[4]={kf_jet,kf_jet,kf_jet,kf_jet};
    Semihard_QCD *group;
    group = new Semihard_QCD(p_beam,p_isr,p_processes->SelectorData(),help,
			     p_processes->ScaleScheme(),
			     p_processes->KFactorScheme());
    std::map<EXTRAXS::XS_Base*,EXTRAXS::XS_Group*> parents;
    for (Process_Map::iterator pit=m_processmap.begin();
	 pit!=m_processmap.end();++pit) {
      EXTRAXS::XS_Group *parent=
	dynamic_cast<EXTRAXS::XS_Group*>(pit->second->Parent());
      parents[pit->second]=parent;
      parent->Remove(pit->second);
      group->Add(pit->second);
    }
    CheckConsistency(group,p_total,m_stop[4],m_start[4],m_sigmahard);
    for (Process_Map::iterator pit=m_processmap.begin();
	 pit!=m_processmap.end();++pit) {
      group->Remove(pit->second);
      parents[pit->second]->Add(pit->second);
    }
    delete group;
  }
  msg_Tracking()<<"Simple_Chain::CalculateTotal(): Pythia mode {"
		<<"\n   \\sigma_{tot} = "
		<<(m_sigmahard*rpa.Picobarn()/1.e9)
		<<" mb @ p_\\perp = "<<m_stop[4]
		<<" GeV\n   \\sigma_{cut} = "
		<<((*p_total)(m_stop[0])*rpa.Picobarn()/1.e9)
		<<" mb p_\\perp = "<<m_stop[0]<<" GeV\n}"<<std::endl;
  p_total->Scale(1.0/m_norm);
  return true;
}

bool Simple_Chain::Initialize()
{
  PROFILE_HERE;
  if (InputPath()=="" && InputFile()=="") return false;
  if (!rpa.gen.Beam1().IsHadron() ||
      !rpa.gen.Beam2().IsHadron()) return false;
  CleanUp();
  Data_Reader *reader = new Data_Reader(" ",";","!","=");
  reader->AddWordSeparator("\t");
  reader->SetInterprete(true);
  reader->SetInputPath(InputPath());
  reader->SetInputFile(InputFile());
  if (!ReadInData()) return false;
  std::string xsfile=std::string("XS.dat");
  reader->ReadFromFile(xsfile,"XS_FILE");
  SetInputFile(xsfile,1);
  double stop, exponent;
  if (!reader->ReadFromFile(stop,"SCALE_MIN")) stop=Stop(0);
  if (reader->ReadFromFile(exponent,"RESCALE_EXPONENT")) {
    double scale;
    if (!reader->ReadFromFile(scale,"REFERENCE_SCALE")) scale=1960.0;
    stop*=pow(m_ecms/scale,exponent);
  }
  SetStop(stop,0);
  SetStop(stop,4); 
  if (m_regulate) {
    SetStop(rpa.gen.Accu()*stop,4);
    // // Uncomment for cross-check vs. PYHTIA
    // SetStop(0.08*stop,0);
  }
  if (!reader->ReadFromFile(m_check,"CHECK_CONSISTENCY")) m_check=0;
  if (!reader->ReadFromFile(m_vegas,"VEGAS_MI")) m_vegas=0;
  if (!reader->ReadFromFile(m_maxreduction,"MI_MAX_REDUCTION")) 
    m_maxreduction=10.0;
  std::string function;
  if (reader->ReadFromFile(function,"PROFILE_FUNCTION")) {
    std::vector<double> parameters;
    if (reader->VectorFromFile(parameters,"PROFILE_PARAMETERS")) {
      p_profile = Profile_Function_Base::SelectProfile(function,parameters);
    }
  }
  int jetveto(1);
  if (!reader->ReadFromFile(jetveto,"JET_VETO")) jetveto=1;
  delete reader;
  SetJetVeto(jetveto);
  if (!CreateGrid()) {
    CleanUp();
    THROW(critical_error,"Grid creation failed.");
  }
  if (!CalculateTotal()) {
    CleanUp();
    THROW(critical_error,"Determination of \\sigma_{tot} failed.");
  }
  if (!SetUpInterface()) {
    CleanUp();
    THROW(critical_error,"Phasespace setup failed.");
  }
  if (p_profile!=NULL) {
    if (!p_profile->CalculateOMean(m_sigmahard/m_norm)) {
      CleanUp();
      THROW(critical_error,"Determination of <\\tilde{O}> failed.");
    }
  }
  return true;
}

void Simple_Chain::SetISRRange()
{
  m_isrspkey[0]=p_isr->SprimeMin();
  m_isrspkey[1]=p_isr->SprimeMax();
  m_isrykey[0]=p_isr->YMin();
  m_isrykey[1]=p_isr->YMax();
  p_isr->SetSprimeMin(4.0*m_last[0]*m_last[0]);
  p_isr->SetSprimeMax(4.0*m_last[2]*m_last[3]);
  p_isr->SetYMin(log(m_last[0]/m_last[3]));
  p_isr->SetYMax(log(m_last[2]/m_last[0]));
}

void Simple_Chain::ResetISRRange()
{
  p_isr->SetYMin(m_isrykey[0]);
  p_isr->SetYMax(m_isrykey[1]);
  p_isr->SetSprimeMax(m_isrspkey[0]);
  p_isr->SetSprimeMin(m_isrspkey[1]);
}

bool Simple_Chain::CreateMomenta()
{
  PROFILE_HERE;
  m_filledblob=false;
  if (p_processes==NULL) {
    THROW(fatal_error,"Multiple interactions are not initialized");
  }
  m_inparticles.Clear();
  m_outparticles.Clear();
  if (m_processmap.find(m_selected)!=m_processmap.end()) {
    EXTRAXS::XS_Base *selected=m_processmap[m_selected];
    selected->Parent()->SetSelected(selected);
    p_processes->SetSelected(selected->Parent());
    double weight=1.;
    size_t pstrials=0, trials=0;
    Amisic_Histogram<double> *cur=m_differentials[m_selected];
    double max=cur->BinMax(m_last[0]);
    p_fsrinterface->SetTrigger(false);
    while (++pstrials<m_maxtrials) {
      Blob_Data_Base *data=selected->
	WeightedEvent(PHASIC::psm::no_lim_isr);
      if (data!=NULL) {
	weight=data->Get<PHASIC::Weight_Info>().weight;
	trials=data->Get<PHASIC::Weight_Info>().ntrial;
	delete data;
	if (weight>max) {
	  msg_Tracking()<<"Simple_Chain::CreateMomenta(): "
			<<"Weight exceeded maximum.\n"
			<<"   Setting new maximum "
			<<max<<" -> "<<weight<<std::endl;
	  m_differentials[m_selected]->SetBinMax(m_last[0],weight);
	}
	bool take(true);
	double mass(0.0);
	Vec4D sum;
	for (size_t j=0;j<selected->NIn();++j) {
	  sum+=selected->Momenta()[j];
	  mass+=selected->Flavours()[j].PSMass();
	  if (selected->Momenta()[j][0]<=selected->Flavours()[j].PSMass()) {
	    take=false;
	    break;
	  }
	}
	if (!take || sum.Mass()<mass) continue;
	mass=0.0;
	for (size_t j=selected->NIn();j<selected->NVector();++j) {
	  mass+=selected->Flavours()[j].PSMass();
	  if (selected->Momenta()[j][0]<=selected->Flavours()[j].PSMass()) {
	    take=false;
	    break;
	  }
	}
	if (!take || sum.Mass()<mass) continue;
	if (p_fsrinterface->Trigger()) {
	  double rn=ran.Get();
	  if (weight*m_maxreduction>=max*rn) {
	    if (weight*m_maxreduction<max) break;
	    double value=cur->BinExtra(m_last[0]);
	    if (value>0.0) {
	      if (value>=1.0 || (value<1.0 && value>ran.Get())) {
		cur->SetBinExtra(m_last[0],Max(0.0,value-1.0));
		m_spkey[3]=
		  Max(cur->BinExtra(m_last[0],1),
			      4.0*(1.0+s_epsilon)*m_last[0]*m_last[0]);
		m_ykey[2]=cur->BinExtra(m_last[0],2);
		double logtau=log(m_spkey[3]/m_spkey[2]);
		if (-logtau<m_ykey[2]) m_ykey[2]=-logtau;
		else if (m_ykey[2]<logtau) m_ykey[2]=logtau;
		msg_Debugging()<<"hit "<<m_selected<<" "<<m_last[0]<<" "
			       <<value<<" "<<cur->BinExtra(m_last[0])
			       <<" "<<m_spkey[3]<<" "<<m_ykey[2]<<"\n";
		SetISRRange();
		p_isr->SetLimits();
		selected->WeightedEvent(PHASIC::psm::no_lim_isr|
					PHASIC::psm::no_dice_isr);
		ResetISRRange();
		cur->AddBinExtra(m_last[0],1.0,3);
	      }
	      else {
		msg_Debugging()<<"no hit "<<m_selected<<" "
			       <<m_last[0]<<" "<<value<<" "
			       <<cur->BinExtra(m_last[0])<<" "
			       <<m_spkey[3]<<" "<<m_ykey[2]<<"\n";
		if (value<1.0) cur->SetBinExtra(m_last[0],0.0);
		return false;
	      }
	    }
	    else {
	      msg_Debugging()<<"set "<<m_selected<<" "<<m_last[0]<<" "
			     <<value<<" "<<m_spkey[3]<<" "
			     <<m_ykey[2]<<" "
			     <<weight*m_maxreduction/max<<"\n";
	      double overflow=weight*m_maxreduction/max;
	      if (overflow>m_maxreduction) {
		msg_Tracking()<<"Simple_Chain::CreateMomenta(): "
			      <<"overflow = "<<overflow<<" > "
			      <<"m_maxreduction = "<<m_maxreduction
			      <<std::endl;
		cur->SetBinMax(m_last[0],weight);
	      }
	      cur->SetBinExtra(m_last[0],overflow-1.0);
	      cur->SetBinExtra(m_last[0],m_spkey[3],1);
	      cur->SetBinExtra(m_last[0],m_ykey[2],2);
	      cur->SetBinExtra(m_last[0],1.0,3);
	    }
	    break;
	  }
	}
      }
    }
    for (size_t j=0;j<selected->NIn();++j) 
      m_last[2+j]-=2.0*selected->Momenta()[j][0]/m_ecms;
    selected->SetColours(selected->Momenta());
    Particle *particle;
    for (size_t j=0;j<selected->NIn();++j) {
      particle = new Particle(0,selected->Flavours()[j]);
      particle->SetMomentum(selected->Momenta()[j]);
      particle->SetFlow(1,selected->Colours()[j][0]);
      particle->SetFlow(2,selected->Colours()[j][1]);
      particle->SetStatus(part_status::active);
      m_inparticles.push_back(particle);
    }
    for (size_t j=selected->NIn();j<selected->NIn()+selected->NOut();++j) {
      particle = new Particle(0,selected->Flavours()[j]);
      particle->SetMomentum(selected->Momenta()[j]);
      particle->SetFlow(1,selected->Colours()[j][0]);
      particle->SetFlow(2,selected->Colours()[j][1]);
      particle->SetStatus(part_status::active);
      m_outparticles.push_back(particle);
    }
    m_filledblob=true;
    return true;
  }
  msg_Error()<<"Simple_Chain::CreateMomenta(..): "
	     <<"Cannot create momentum configuration."<<std::endl;
  return false;
}

bool Simple_Chain::DiceProcess()
{
  PROFILE_HERE;
  if (m_differentials.size()==0) return false;
  while (true) {
  if (!DiceOrderingParameter()) return false;
  if (m_dicedparameter) m_dicedparameter=false;
  else {
    m_dicedprocess=false;
    return true;
  }
  if (m_last[2]*m_last[3]<=m_last[0]*m_last[0]) {
    m_dicedprocess=false;
    return true;
  }
  p_fsrinterface->SetValue(m_last[0]);
  Sort_Map sorter;
  Sort_Map::key_type norm=0.0, cur=0.0;
  for (Amisic_Histogram_Map::iterator hit=m_differentials.begin();
       hit!=m_differentials.end();++hit) {
    cur=(*hit->second)(m_last[0]);
    sorter.insert(Sort_Map::value_type(cur,hit->first));
    norm+=cur;
  }
  if (norm==0.0) {
    msg_Error()<<METHOD<<"(): Warning. Zero cross section."<<std::endl;
    continue;
  }
  double rannr=ran.Get();
  cur=0.0;
  m_selected="";
  for (Sort_Map::iterator sit=sorter.begin();
       sit!=sorter.end();++sit) {
    for (;sit!=sorter.upper_bound(sit->first);++sit) {
      if ((cur+=sit->first/norm)>rannr) {
	m_selected=sit->second;
	break;
      }
    }
    if (m_selected!="") break;
  }
	SetISRRange();
	if (!CreateMomenta()) continue;
	ResetISRRange();
	m_dicedprocess=true;
	return m_filledblob;
  }
  THROW(critical_error,"Internal Error. Could not select any process.");
  return false;
}

bool Simple_Chain::DiceEnhanceFactor()
{
  if (p_profile==NULL) return true;
  double b=0.0;
  double last=(*p_total)(m_last[0]);
  do {
    b=p_profile->DiceImpactParameter();
    m_enhance=(*p_profile)(b)/p_profile->OMean();
  } while (exp(-m_enhance*last)<=ran.Get());
  msg_Tracking()<<"Simple_Chain::DiceEnhanceFactor(): { profile '"
		<<p_profile->Type()
		<<"'\n   m_last[0]  = "<<m_last[0]<<"\n   p(k_t^2)   = "<<last
		<<"\n   b          = "<<b
		<<"\n   e(b)       = "<<m_enhance<<"\n   e(b)_{min} = "
		<<p_profile->OMin()/p_profile->OMean()<<"\n   e(b)_{max} = "
		<<p_profile->OMax()/p_profile->OMean()<<"\n}"<<std::endl;
  return true;
}

bool Simple_Chain::DiceOrderingParameter()
{ 
  PROFILE_HERE;
  if (m_last[0]<=m_stop[0]) {
    msg_Error()<<"Simple_Chain::DiceOrderingParameter(): "
	       <<"Value exceeded minimum: last = "<<m_last[0]
	       <<" vs. stop = "<<m_stop[0]<<std::endl;
    s_stophard=true;
    return false;
  }
  if (s_cleaned) if (!DiceEnhanceFactor()) {
    s_stophard=true;
    return false;
  }
  m_last[0]=(*p_total)[(*p_total)
 		       (m_last[0])-log(ran.Get())/m_enhance]; 
  msg_Debugging()<<METHOD<<"(): new p_T = "<<m_last[0]<<"\n";
  s_cleaned=false;
  if (m_last[0]<=m_stop[0]) { 
    m_dicedparameter=false;
    s_stophard=true;
    return true;
  }
  m_dicedparameter=true;
  s_stophard=false;
  return true;
}

void Simple_Chain::Reset()
{
  for (unsigned int i=0;i<4;++i) m_last[i]=m_start[i];
  for (Amisic_Histogram_Map::const_iterator hit(m_differentials.begin());
       hit!=m_differentials.end();++hit) hit->second->StoreData();
}

void Simple_Chain::Update(const MI_Base *mibase)
{
  return;
}

bool Simple_Chain::ReadInStatus(const std::string &path) 
{
  msg_Info()<<METHOD<<"(): Reading status from '"
	    <<path<<m_pathextra<<"'."<<std::endl;
  p_gridcreator->SetOutputPath(path+m_pathextra);
  if (!p_gridcreator->ReadInGrid()) {
    msg_Error()<<METHOD<<"(): No status stored in '"
	       <<path<<m_pathextra<<"'"<<std::endl;
    return false;
  }
  return true;
}

void Simple_Chain::PrepareTerminate() 
{
  std::string path(rpa.gen.Variable("SHERPA_STATUS_PATH"));
  if (path=="") return;
  for (Amisic_Histogram_Map::const_iterator hit(m_differentials.begin());
       hit!=m_differentials.end();++hit) hit->second->RestoreData();
  CopyFile(InputPath()+m_selectorfile,path+"/"+m_selectorfile);
  path+="/"+m_pathextra;
  MakeDir(path,true);
  p_gridcreator->WriteOutGrid(String_Vector(),path);
}

bool Simple_Chain::VetoProcess(ATOOLS::Blob *blob)
{
  if (s_soft==NULL) return false;
  bool veto=ran.Get()>s_xsnd/s_xstot;
  if (veto) {
    s_soft->SetStart(m_stop[0],0); 
    s_soft->SetStart((*p_differential)(m_stop[0]),2); 
  }
  return s_stophard=veto;
}
