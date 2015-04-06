#include "AMISIC++/Model/Simple_Chain.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "EXTRA_XS/Main/Single_Process.H"
#include "AMISIC++/Tools/Semihard_QCD.H"
#include "AMISIC++/Tools/MPI_KFactor_Setter.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/My_Limits.H"
#include "PHASIC++/Channels/Vegas.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "PDF/Remnant/Remnant_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/Weight_Info.H"

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
  m_maxreduction(1.0), m_sigma_nd_fac(1.0),
  m_xsextension("_xs.dat"),
  m_resdir(""), m_ressuffix(""),
  p_model(NULL),
  p_beam(NULL), p_isr(NULL), p_profile(NULL), m_maxtrials(1000),
  m_ecms(rpa->gen.Ecms()), m_external(false), m_regulate(false)
{
  Init();
}

Simple_Chain::Simple_Chain(MODEL::Model_Base *const model,
			   BEAM::Beam_Spectra_Handler *const beam,
			   PDF::ISR_Handler *const isr):
  MI_Base("Simple Chain",MI_Base::HardEvent,5,4,1),
  p_differential(NULL), p_total(NULL), m_norm(1.0), m_enhance(1.0),
  m_maxreduction(1.0), m_sigma_nd_fac(1.0),
  m_xsextension("_xs.dat"),
  p_model(model),
  p_beam(beam), p_isr(isr), p_profile(NULL), m_maxtrials(1000),
  m_ecms(rpa->gen.Ecms()), m_external(true), m_regulate(false)
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
  m_spkey.Assign("s' isr",5,0,PHASIC::Phase_Space_Handler::GetInfo());
  m_ykey.Assign("y isr",3,0,PHASIC::Phase_Space_Handler::GetInfo());
  m_xkey.Assign("x isr",5,0,PHASIC::Phase_Space_Handler::GetInfo());
  m_isrspkey.Assign("s' isr mi",3,0,PHASIC::Phase_Space_Handler::GetInfo());
  m_isrykey.Assign("y isr mi",2,0,PHASIC::Phase_Space_Handler::GetInfo());
  p_remnants[1]=p_remnants[0]=NULL;
  p_gridcreator=NULL;
  m_xsec_output = std::string("MPI_Cross_Sections.dat");
}

Simple_Chain::~Simple_Chain()
{
  CleanUp();
  delete p_read;
}

void Simple_Chain::CleanUp() 
{
  if (p_gridcreator!=NULL) {
    exh->RemoveTerminatorObject(this);
    delete p_gridcreator;
    p_gridcreator=NULL;
  }
  for (size_t i=0; i<p_processes.size(); ++i) {
    p_processes[i]->SetFSRMode(3);
    p_processes[i]->CreateFSRChannels();
    delete p_fsrinterface[i];
  }
  p_fsrinterface.clear();
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
  std::string outputpath;
  MyStrStream converter;
  if (m_resdir=="") {
    std::string help[2];
    converter<<rpa->gen.Bunch(0);
    converter>>help[0];
    converter.clear();
    converter<<rpa->gen.Bunch(1);
    converter>>help[1];
    outputpath=std::string("MIG_")+help[0]+help[1]+
      std::string("_")+ToString(rpa->gen.Ecms());
    if (m_regulate) {
      outputpath+=std::string("_")+m_regulator[0];
      for (size_t i=0;i<m_regulation.size();++i) {
        outputpath+=std::string("_")+ToString(m_regulation[i]);
      }
    }
    if (p_isr->PDF(0)->Type()!=p_isr->PDF(1)->Type()) {
      outputpath+=std::string("_")+p_isr->PDF(0)->Type();
    }
    outputpath+=std::string("_")+p_isr->PDF(0)->Type()+std::string("_")
                +ToString(static_cast<MODEL::Running_AlphaS*>
                       (p_model->GetScalarFunction("alpha_S"))->Order())
                +m_ressuffix+"/";
  }
  else outputpath=m_resdir+m_ressuffix+"/";
  msg_Debugging()<<METHOD<<"(): path = "
                 <<OutputPath()+m_pathextra+outputpath<<std::endl;
  SetOutputPath(OutputPath()+m_pathextra+outputpath);
  m_pathextra=outputpath;
  return true;
}

bool Simple_Chain::ReadInData()
{
  Data_Reader *reader = new Data_Reader(" ",";","!","=");
  reader->AddComment("#");
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
      m_regulation=std::vector<double>(1,2.68);
    double exponent, scale;
    if (!reader->ReadFromFile(exponent,"RESCALE_EXPONENT")) exponent=0.244;
    if (!reader->ReadFromFile(scale,"REFERENCE_SCALE")) scale=1800.0;
    m_regulation[0]*=pow(m_ecms/scale,exponent);
  }
  m_heavy_flavour = reader->GetValue<int>("MI_HEAVY_FLAVOUR",1);
  if (!reader->ReadFromFile(m_error,"PS_ERROR")) m_error=1.e-2;
  if (!reader->ReadFromFile(m_pathextra,"PATH_EXTRA")) m_pathextra="";
  m_sigma_nd_fac = reader->GetValue<double>("SIGMA_ND_FACTOR",0.335);
  m_resdir = reader->GetValue<std::string>("MI_RESULT_DIRECTORY","");
  m_ressuffix = reader->GetValue<std::string>("MI_RESULT_DIRECTORY_SUFFIX","");
  GeneratePathName();
  delete reader;
  return true;
}

bool Simple_Chain::CreateGrid()
{
  bool vegas=PHASIC::Vegas::OnExternal();
  PHASIC::Vegas::SetOnExternal(m_vegas);
  double min=Min(m_stop[0],m_stop[4]);
  p_isr->SetFixedSprimeMin(4.0*min*min);
  p_isr->SetFixedSprimeMax(4.0*m_start[0]*m_start[0]);
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader(" ",";","!","=");
  reader->AddComment("#");
  reader->AddWordSeparator("\t");
  reader->SetInputPath(InputPath());
  reader->SetInputFile(InputFile(2));
  if (!reader->ReadFromFile(m_selectorfile,"MI_SELECTOR_FILE")) 
    m_selectorfile="MICuts.dat";
  delete reader;
  InitializeProcessList(Flavour(kf_jet), Flavour(kf_jet),
                        Flavour(kf_jet), Flavour(kf_jet));
  if (m_heavy_flavour && !Flavour(kf_jet).Includes(Flavour(kf_b))) {
    InitializeProcessList(Flavour(kf_jet), Flavour(kf_jet),
                          Flavour(kf_b), Flavour(kf_b, true));
  }
  if (m_heavy_flavour && !Flavour(kf_jet).Includes(Flavour(kf_c))) {
    InitializeProcessList(Flavour(kf_jet), Flavour(kf_jet),
                          Flavour(kf_c), Flavour(kf_c, true));
  }
  std::vector<EXTRAXS::Process_Group*> procs;
  for (size_t i=0; i<p_processes.size(); ++i) procs.push_back(p_processes[i]);
  p_gridcreator = new Grid_Creator(&m_differentials,procs);
  p_gridcreator->SetGridXMin(min);
  p_gridcreator->SetGridXMax(m_ecms/2.0);
  p_gridcreator->ReadInArguments(InputFile(),InputPath());
  p_gridcreator->SetXSExtension(m_xsextension);
  p_gridcreator->SetOutputPath(OutputPath());
  if (!p_gridcreator->InitializeCalculation()) {
    msg_Error()<<METHOD<<"(): Initialization failed! Abort."<<std::endl;
    return false;
  }
  if (!p_gridcreator->ReadInGrid()) {
    p_gridcreator->CreateGrid();
  }
  PHASIC::Vegas::SetOnExternal(vegas);
  exh->AddTerminatorObject(this);
  Reset();
  return true;
}

void Simple_Chain::InitializeProcessList(const Flavour& in1,
                                         const Flavour& in2,
                                         const Flavour& out1,
                                         const Flavour& out2)
{
  PHASIC::Process_Info pi;
  pi.m_ii.m_ps.push_back(PHASIC::Subprocess_Info(in1,"",""));
  pi.m_ii.m_ps.push_back(PHASIC::Subprocess_Info(in2,"",""));
  pi.m_fi.m_ps.push_back(PHASIC::Subprocess_Info(out1,"",""));
  pi.m_fi.m_ps.push_back(PHASIC::Subprocess_Info(out2,"",""));
  pi.m_maxcpl[1]=pi.m_mincpl[1]=0;
  pi.m_maxcpl[0]=pi.m_mincpl[0]=2;
  pi.m_scale=p_read->GetValue<std::string>("MPI_SCALE","MPI");
  pi.m_kfactor=p_read->GetValue<std::string>("MPI_KFACTOR","MPI");
  pi.m_coupling="Alpha_QCD 1";
  pi.m_mpiprocess=true;
  p_processes.push_back(new Semihard_QCD(p_read));
  p_processes.back()->Init(pi,p_beam,p_isr);
  msg_Info()<<METHOD<<"(): Init processes ";
  if (!p_processes.back()->Get<EXTRAXS::Process_Group>()->ConstructProcesses())
      THROW(fatal_error,"Cannot initialize MPI simulation.");
  msg_Info()<<" done."<<std::endl;
  p_processes.back()->SetScale(PHASIC::Scale_Setter_Arguments
                        (p_model,pi.m_scale,pi.m_coupling));
  p_processes.back()->SetKFactor(PHASIC::KFactor_Setter_Arguments(pi.m_kfactor));
  p_processes.back()->InitPSHandler(m_error,"","");
  p_processes.back()->SetGenerator(p_processes.back());
  for (size_t i(0);i<p_processes.back()->Size();++i)
    m_processmap[(*p_processes.back())[i]->Name()]=(*p_processes.back())[i];
}

bool Simple_Chain::SetUpInterface()
{
  for (size_t i=0; i<p_fsrinterface.size(); ++i) delete p_fsrinterface[i];
  p_fsrinterface.resize(p_processes.size());
  for (size_t i=0; i<p_fsrinterface.size(); ++i) {
    Flavour flavour[4]={p_processes[i]->Flavours()[0],
                        p_processes[i]->Flavours()[1],
                        p_processes[i]->Flavours()[2],
                        p_processes[i]->Flavours()[3]};
    if (p_fsrinterface[i]!=NULL) delete p_fsrinterface[i];
    p_fsrinterface[i] = new FSR_Channel(2,2,flavour,
                                        p_total->XAxis()->Variable()->Name());
    p_processes[i]->InitIntegrators();
    p_processes[i]->CreateISRChannels();
    p_processes[i]->SetFSRInterface(p_fsrinterface[i]);
    p_processes[i]->SetFSRMode(2);
    p_processes[i]->CreateFSRChannels();
  }
  return true;
}

void Simple_Chain::CalculateSigmaND()
{
  if(s_kftable.find(111)==s_kftable.end()) // if not initialized yet
    s_kftable[111]=new Particle_Info(111,0.134976,7.8486e-09,0,0,1,0,"pi","pi");
  double eps=0.0808, eta=-0.4525, X=21.70, Y=56.08, b=2.3;
  if (p_isr->Flav(0).IsAnti()^p_isr->Flav(1).IsAnti()) Y=98.39;
  double s=sqr(rpa->gen.Ecms());
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
  SetNorm(m_sigma_nd_fac*(xstot-xsel-2.0*xssd-xsdd)*1.0e9/rpa->Picobarn());
}

int Simple_Chain::CalculateTotal()
{
  if (m_differentials.size()==0) return -1;
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
	    <<"\\sigma_{hard} = "<<(m_sigmahard*rpa->Picobarn()/1.e9)
	    <<" mb\n   at "<<xaxis->Variable()->Name()<<"_{min} = "
	    <<m_stop[4]<<" GeV\n}"<<std::endl;
  CalculateSigmaND();
  if (m_sigmahard<m_norm) {
    msg_Error()<<"Simple_Chain::CalculateTotal(): "<<om::red
	       <<"\\sigma_{hard} = "
	       <<(m_sigmahard*rpa->Picobarn()/1.e9)
	       <<" mb < \\sigma_{nd} = "
	       <<(m_norm*rpa->Picobarn()/1.e9)
	       <<" mb !"<<om::reset<<std::endl;
    return 0;
  }
  p_total->Scale(1.0/m_norm);
  return 1;
}

bool Simple_Chain::Initialize()
{
  if (InputPath()=="" && InputFile()=="") return false;
  if (!rpa->gen.Beam1().IsHadron() ||
      !rpa->gen.Beam2().IsHadron()) return false;
  CleanUp();
  p_read = new Data_Reader(" ",";","!","=");
  p_read->AddComment("#");
  p_read->AddWordSeparator("\t");
  p_read->SetInterprete(true);
  p_read->SetInputPath(InputPath());
  p_read->SetInputFile(InputFile());
  if (!ReadInData()) return false;
  std::string xsfile=std::string("XS.dat");
  p_read->ReadFromFile(xsfile,"XS_FILE");
  SetInputFile(xsfile,1);
  double stop, exponent, scale, pt0, pt0exp;
  if (!p_read->ReadFromFile(pt0,"TURNOFF")) pt0=0.0;
  if (!p_read->ReadFromFile(stop,"SCALE_MIN")) stop=2.44;
  if (!p_read->ReadFromFile(pt0exp,"TURNOFF_EXPONENT")) pt0exp=0.244;
  if (!p_read->ReadFromFile(exponent,"RESCALE_EXPONENT")) exponent=0.244;
  if (!p_read->ReadFromFile(scale,"REFERENCE_SCALE")) scale=1800.0;
  MPI_KFactor_Setter::SetPT0(pt0*pow(m_ecms/scale,pt0exp));
  stop*=pow(m_ecms/scale,exponent);
  SetStop(stop,0);
  SetStop(stop,4); 
  if (m_regulate) {
    SetStop(rpa->gen.Accu()*stop,4);
    // // Uncomment for cross-check vs. PYHTIA
    // SetStop(0.08*stop,0);
  }
  if (!p_read->ReadFromFile(m_check,"CHECK_CONSISTENCY")) m_check=0;
  if (!p_read->ReadFromFile(m_vegas,"VEGAS_MI")) m_vegas=0;
  if (!p_read->ReadFromFile(m_maxreduction,"MI_MAX_REDUCTION")) 
    m_maxreduction=10.0;
  std::string function;
  std::vector<double> parameters;
  if (!p_read->ReadFromFile(function,"PROFILE_FUNCTION")) {
    function="Double_Gaussian";
    if (!p_read->VectorFromFile(parameters,"PROFILE_PARAMETERS")) {
      parameters.push_back(0.76);
      parameters.push_back(0.58);
    }
  }
  else {
    p_read->VectorFromFile(parameters,"PROFILE_PARAMETERS");
  }
  if (function!="None")
    p_profile = Profile_Function_Base::SelectProfile(function,parameters);
  int jetveto(1);
  if (!p_read->ReadFromFile(jetveto,"JET_VETO")) jetveto=1;
  SetJetVeto(jetveto);
  if (!CreateGrid()) {
    CleanUp();
    THROW(critical_error,"Grid creation failed.");
  }
  int res(CalculateTotal());
  if (res==0) {
    CleanUp();
    msg_Error()<<METHOD<<"(): Switching MPI simulation off."<<std::endl;
    return false;
  }
  if (res<0) {
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
  std::ofstream ofile;
  ofile.open(m_xsec_output.c_str());
  ofile<<"MPIs in Sherpa, Model = Amisic: \n"
       <<"   semihard xsec = "<<(m_sigmahard*rpa->Picobarn()/1.e9)<<" mb,\n"
       <<"   non-diffractive xsec = "<<(m_norm*rpa->Picobarn()/1.e9)<<" mb "
       <<"with nd factor = "<<m_sigma_nd_fac<<".\n";
  ofile.close();

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
  m_filledblob=false;
  if (p_processes.empty()) {
    THROW(fatal_error,"Multiple interactions are not initialized");
  }
  m_inparticles.Clear();
  m_outparticles.Clear();
  if (m_processmap.find(m_selected)!=m_processmap.end()) {
    p_xs=m_processmap[m_selected];
    double weight=1.;
    size_t pstrials=0, trials=0;
    Amisic_Histogram<double> *cur=m_differentials[m_selected];
    double max=cur->BinMax(m_last[0]);
    FSR_Channel* fsr = dynamic_cast<FSR_Channel*>
      (p_xs->Parent()->Integrator()->PSHandler()->FSRIntegrator()->Channel(0));
    if (fsr==NULL) THROW(fatal_error, "Internal error.");
    fsr->SetTrigger(false);
    while (++pstrials<m_maxtrials) {
      ATOOLS::Weight_Info *data=p_xs->
	OneEvent(0,PHASIC::psm::no_lim_isr);
      if (data!=NULL) {
	weight=data->m_weight;
	trials=data->m_ntrial;
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
	for (size_t j=0;j<p_xs->NIn();++j) {
	  sum+=p_xs->Integrator()->Momenta()[j];
	  mass+=p_xs->Flavours()[j].Mass();
	  if (p_xs->Integrator()->Momenta()[j][0]
	      <=p_xs->Flavours()[j].Mass()) {
	    take=false;
	    break;
	  }
	}
	if (!take || sum.Mass()<mass) continue;
	mass=0.0;
	for (size_t j=p_xs->NIn();
	     j<p_xs->NIn()+p_xs->NOut();++j) {
	  mass+=p_xs->Flavours()[j].Mass();
	  if (p_xs->Integrator()->Momenta()[j][0]
	      <=p_xs->Flavours()[j].Mass()) {
	    take=false;
	    break;
	  }
	}
	if (!take || sum.Mass()<mass) continue;
	if (fsr->Trigger()) {
	  double rn=ran->Get();
	  if (weight*m_maxreduction>=max*rn) {
	    if (weight*m_maxreduction<max) break;
	    double value=cur->BinExtra(m_last[0]);
	    if (value>0.0) {
	      if (value>=1.0 || (value<1.0 && value>ran->Get())) {
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
		p_isr->SetLimits(m_spkey.Doubles(),m_ykey.Doubles(),
				 m_xkey.Doubles());
		ATOOLS::Weight_Info *info=
		  p_xs->OneEvent(0,PHASIC::psm::no_lim_isr|
				 PHASIC::psm::no_gen_isr);
		if (info) delete info;
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
    if (pstrials==m_maxtrials) return false;
    for (size_t j=0;j<p_xs->NIn();++j) 
      m_last[2+j]-=2.0*p_xs->Integrator()->
	Momenta()[j][0]/m_ecms;
    Particle *particle;
    for (size_t j=0;j<p_xs->NIn();++j) {
      particle = new Particle(0,p_xs->Flavours()[j]);
      particle->SetMomentum(p_xs->Integrator()->Momenta()[j]);
      particle->SetStatus(part_status::active);
      m_inparticles.push_back(particle);
    }
    for (size_t j=p_xs->NIn();j<p_xs->NIn()+p_xs->NOut();++j) {
      particle = new Particle(0,p_xs->Flavours()[j]);
      particle->SetMomentum(p_xs->Integrator()->Momenta()[j]);
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

bool Simple_Chain::GenerateProcess()
{
  if (m_differentials.size()==0) return false;
  while (true) {
  if (!GenerateOrderingParameter()) return false;
  if (m_generatedparameter) m_generatedparameter=false;
  else {
    m_generatedprocess=false;
    return true;
  }
  if (m_last[2]*m_last[3]<=m_last[0]*m_last[0]) {
    m_generatedprocess=false;
    return true;
  }
  for (size_t i=0; i<p_fsrinterface.size(); ++i)
    p_fsrinterface[i]->SetValue(m_last[0]);
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
  double rannr=ran->Get();
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
	m_generatedprocess=true;
	return m_filledblob;
  }
  THROW(critical_error,"Internal Error. Could not select any process.");
  return false;
}

bool Simple_Chain::GenerateEnhanceFactor()
{
  if (p_profile==NULL) return true;
  double b=0.0;
  double last=(*p_total)(m_last[0]);
  do {
    b=p_profile->GenerateImpactParameter();
    m_enhance=(*p_profile)(b)/p_profile->OMean();
  } while (exp(-m_enhance*last)<=ran->Get());
  msg_Tracking()<<"Simple_Chain::GenerateEnhanceFactor(): { profile '"
		<<p_profile->Type()
		<<"'\n   m_last[0]  = "<<m_last[0]<<"\n   p(k_t^2)   = "<<last
		<<"\n   b          = "<<b
		<<"\n   e(b)       = "<<m_enhance<<"\n   e(b)_{min} = "
		<<p_profile->OMin()/p_profile->OMean()<<"\n   e(b)_{max} = "
		<<p_profile->OMax()/p_profile->OMean()<<"\n}"<<std::endl;
  return true;
}

bool Simple_Chain::GenerateOrderingParameter()
{ 
  if (m_last[0]<=m_stop[0]) {
    msg_Error()<<"Simple_Chain::GenerateOrderingParameter(): "
	       <<"Value exceeded minimum: last = "<<m_last[0]
	       <<" vs. stop = "<<m_stop[0]<<std::endl;
    s_stophard=true;
    return false;
  }
  if (s_cleaned) if (!GenerateEnhanceFactor()) {
    s_stophard=true;
    return false;
  }
  msg_Debugging()<<METHOD<<"(): old p_T = "<<m_last[0]<<", ";
  m_last[0]=(*p_total)[(*p_total)
 		       (m_last[0])-log(ran->Get())/m_enhance]; 
  msg_Debugging()<<"new p_T = "<<m_last[0]<<"\n";
  s_cleaned=false;
  if (m_last[0]<=m_stop[0]) { 
    m_generatedparameter=false;
    s_stophard=true;
    return true;
  }
  m_generatedparameter=true;
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
  std::string path(rpa->gen.Variable("SHERPA_STATUS_PATH"));
  if (path=="") return;
  for (Amisic_Histogram_Map::const_iterator hit(m_differentials.begin());
       hit!=m_differentials.end();++hit) hit->second->RestoreData();
  path+="/"+m_pathextra;
  MakeDir(path,true);
  p_gridcreator->WriteOutGrid(String_Vector(),path);
}

bool Simple_Chain::VetoProcess(ATOOLS::Blob *blob)
{
  if (s_soft==NULL) return false;
  bool veto=ran->Get()>s_xsnd/s_xstot;
  if (veto) {
    s_soft->SetStart(m_stop[0],0); 
    s_soft->SetStart((*p_differential)(m_stop[0]),2); 
  }
  return s_stophard=veto;
}
