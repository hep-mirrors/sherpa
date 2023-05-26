#include "PHASIC++/Main/Phase_Space_Handler.H"

#include "PHASIC++/Main/Phase_Space_Integrator.H"
#include "PHASIC++/Main/Channel_Creator.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Channels/FSR_Channels.H"
#include "PHASIC++/Channels/Rambo.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Single_Process.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Message.H"  
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Data_Writer.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Weight_Info.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace BEAM;
using namespace PDF;
using namespace std;

Integration_Info *PHASIC::Phase_Space_Handler::p_info=NULL;

Phase_Space_Handler::Phase_Space_Handler(Process_Integrator *proc,double error,
                                         const std::string eobs,
                                         const std::string efunc): m_name(proc->Process()->Name()), p_process(proc), p_active(proc),
      p_integrator(NULL), p_beamhandler(proc->Beam()), m_pspoint(Phase_Space_Point(this)),
      p_isrhandler(proc->ISR()), p_yfshandler(proc->YFS()), p_flavours(proc->Process()->Flavours()),
      m_nin(proc->NIn()), m_nout(proc->NOut()), m_nvec(m_nin + m_nout),
      m_initialized(false), m_sintegrator(0), m_killedpoints(0),
      m_printpspoint(false), m_enhanceObs(eobs), m_enhanceFunc(efunc) {
  RegisterDefaults();
  InitParameters(error);
  p_process->SetPSHandler(this);
  FSR_Channels * fsrchannels =
    new FSR_Channels(this,"fsr_"+p_process->Process()->Name());
  SetFSRIntegrator(fsrchannels);
  p_lab.resize(m_nvec);
}

Phase_Space_Handler::~Phase_Space_Handler() { delete p_integrator; }

bool Phase_Space_Handler::CreateIntegrators() {
  Channel_Creator channelcreator(this);
  if (channelcreator()) {
    m_pspoint.Init();
    m_psenhance.Init(this);
    m_psenhance.SetObservable(m_enhanceObs, p_process->Process());
    m_psenhance.SetFunction(m_enhanceFunc, p_process->Process());
    m_enhanceweight =
        m_psenhance.Factor(p_process->TotalXS());
    return true;
  } else
    THROW(fatal_error, "Creation of integrators failed.")
}

double Phase_Space_Handler::Integrate() 
{
  CheckSinglePoint();
  if (p_process->Points()>0)
    return p_process->TotalXS()*rpa->Picobarn();
  p_integrator = new Phase_Space_Integrator(this);
  if (m_nin==1) return p_integrator->CalculateDecay(m_error);
  if (m_nin==2) return p_integrator->Calculate(m_error,m_abserror,m_fin_opt);
  return 0.;
}

Weights_Map
Phase_Space_Handler::Differential(Process_Integrator *const process,
				  Variations_Mode varmode,
				  const psmode::code mode)
{
  m_cmode  = mode;
  p_active = process;
  m_wgtmap = 0.0;
  // check for failure to generate a meaningful phase space point
  if (!process->Process()->GeneratePoint() ||
      !m_pspoint(process,m_cmode))
    return 0.0;
  for (auto p : p_lab) {
    if (p.Nan()) return 0.0;
  }
  // phase space trigger, calculate and construct weights
  if (process->Process()->Trigger(p_lab)) {
    if (!p_active->Process()->Selector()->Pass()) return 0.0;
    m_psweight = CalculatePS();
    m_wgtmap   = CalculateME(varmode);
    m_wgtmap  *= m_psweight;
    m_wgtmap  *= (m_enhanceweight = m_psenhance.Factor(p_process->TotalXS()));
    m_wgtmap  *= (m_ISsymmetryfactor = m_pspoint.ISSymmetryFactor());
    p_lab      = process->Momenta();
    if (m_printpspoint || msg_LevelIsDebugging()) PrintIntermediate();
    ManageWeights(m_psweight*m_enhanceweight*m_ISsymmetryfactor);
  }
  // trigger failed, return 0.
  else ManageWeights(0.0);
  // stability checks may lead to event weight set to 0 in case of failure
  if (!CheckStability())
    m_wgtmap *= 0.;
  return m_wgtmap;
}


void Phase_Space_Handler::PrintIntermediate() {
  size_t precision(msg->Out().precision());
  msg->SetPrecision(15);
  msg_Out()<<"==========================================================\n"
	   <<p_active->Process()->Name()
	   <<"  ME = "<<m_wgtmap.Nominal()
           <<" ,  PS = "<<m_psweight
           <<" ,  enh = "<<m_enhanceweight
           <<"  ->  "
	   <<m_wgtmap.Nominal()*m_psweight*m_enhanceweight<<std::endl;
  if (p_active->Process()->GetSubevtList()) {
    NLO_subevtlist * subs(p_active->Process()->GetSubevtList());
    for (size_t i(0);i<subs->size();++i) msg_Out()<<(*(*subs)[i])<<"\n";
  }
  for (size_t i(0);i<p_lab.size();++i)
    msg_Out()<<"  p_lab["<<i<<"]=Vec4D"<<p_lab[i]<<";"<<std::endl;
  msg_Out()<<"==========================================================\n";
  msg->SetPrecision(precision);
}

void Phase_Space_Handler::ManageWeights(const double & factor) {
  if (factor!=0.0) {
    ME_Weight_Info* wgtinfo=p_active->Process()->GetMEwgtinfo();
    if (wgtinfo) { (*wgtinfo) *= factor; }
  }
  NLO_subevtlist* nlos=p_active->Process()->GetSubevtList();
  if (nlos) { (*nlos) *= factor; (*nlos).MultMEwgt(factor); }
}

bool Phase_Space_Handler::CheckStability() {
  // meaningful result - no problems down the line
  if (p_active->TotalXS() &&
      dabs(m_wgtmap.Nominal()/p_active->TotalXS())>dabs(m_thkill)) {
    if (m_thkill<0.0) {
      msg_Info()<<METHOD<<"(): Skip point in '"<<p_active->Process()->Name()<<"', "
		<<"weight = "<<m_wgtmap.Nominal()*rpa->Picobarn()<<", thkill = "<<m_thkill<<",\n"
		<<"   totalxs = "<<p_active->TotalXS()<<", result = "<<m_wgtmap.Nominal()<<".\n";
      return false;
    }
    // outputb tricky phase space point for further analysis, if necessary
    ATOOLS::MakeDir("stability");
    std::ofstream sf(("stability/"+p_active->Process()->Name()+
		      "_"+rpa->gen.Variable("RNG_SEED")).c_str(),
		     std::ios_base::app);
    sf.precision(16);
    sf<<"(P"<<m_killedpoints<<"){ # w = "
      <<m_wgtmap<<", ME = "<<m_wgtmap.Nominal()/m_psweight<<", PS = "<<m_psweight<<"\n";
    for (size_t i(0);i<p_lab.size();++i) sf<<"  p_lab["<<i<<"]=Vec4D"<<p_lab[i]<<";\n";
    sf<<"}(P"<<m_killedpoints<<");\n";
    ++m_killedpoints;
    ManageWeights(0.);
    return false;
  }
  return true;
}


Weight_Info *Phase_Space_Handler::OneEvent(Process_Base *const proc,
                                           Variations_Mode varmode,
                                           int mode)
{
  if (proc==NULL) THROW(fatal_error,"No process.");
  Process_Integrator *cur(proc->Integrator());
  p_isrhandler->SetRunMode(1);
  auto wgtmap = Differential(cur, varmode, (psmode::code)mode);
  if (wgtmap.IsZero() || IsBad(wgtmap.Nominal()))
    return NULL;
  cur->SetMomenta(p_lab);
  int fl1(0), fl2(0);
  double x1(0.0), x2(0.0), xf1(0.0), xf2(0.0), mu12(0.0), mu22(0.0), dxs(0.0);
  dxs=p_active->Process()->Get<PHASIC::Single_Process>()->LastXS();
  // if(p_yfshandler->GetMode()==1){
  //   p_yfshandler->CalculateBeta();
  //   wgtmap*=p_yfshandler->m_real;
  // }
  const int swap(p_isrhandler->Swap());
  fl1=(long int)p_active->Process()->Flavours()[swap];
  fl2=(long int)p_active->Process()->Flavours()[1-swap];
  x1=p_isrhandler->X1();
  x2=p_isrhandler->X2();
  xf1=p_isrhandler->XF1();
  xf2=p_isrhandler->XF2();
  mu12=p_isrhandler->MuF2(0);
  mu22=p_isrhandler->MuF2(1);
  auto res =
      new Weight_Info(wgtmap, dxs, 1.0, fl1, fl2, x1, x2, xf1, xf2, mu12, mu22);
  return res;
}

void Phase_Space_Handler::AddPoint(const double _value)
{
  p_process->AddPoint(_value);
  double value(_value);
  if (p_process->TotalXS()==0.0) value=(_value?1.0:0.0);
  if (value!=0.0) {
    m_pspoint.AddPoint(value);
    m_psenhance.AddPoint(value);
  }
}

void Phase_Space_Handler::WriteOut(const std::string &pID) 
{
  m_pspoint.WriteOut(pID);
  m_psenhance.WriteOut(pID);
  Data_Writer writer;
  writer.SetOutputPath(pID+"/");
  writer.SetOutputFile("Statistics.dat");
  writer.MatrixToFile(m_stats);
}

bool Phase_Space_Handler::ReadIn(const std::string &pID,const size_t exclude) 
{
  msg_Info()<<"Read in channels from directory: "<<pID<<std::endl;
  if (m_pspoint.ReadIn(pID,exclude)) {
    m_psenhance.ReadIn(pID);
    Data_Reader reader;
    reader.SetInputPath(pID+"/");
    reader.SetInputFile("Statistics.dat");
    std::vector<std::vector<double> > stats;
    if (reader.MatrixFromFile(stats,"")) m_stats=stats;
    return true;
  }
  return false;
}

void Phase_Space_Handler::RegisterDefaults() const
{
  Settings& settings = Settings::GetMainSettings();
  settings["IB_THRESHOLD_KILL"].SetDefault(-1.0e12);
  settings["ERROR"].SetDefault(0.01);
  settings["INTEGRATION_ERROR"].SetDefault(settings["ERROR"].Get<double>());
  settings["ABS_ERROR"].SetDefault(0.0);
  settings["FINISH_OPTIMIZATION"].SetDefault(true);
  settings["PRINT_PS_POINTS"].SetDefault(false);
  settings["PS_PT_FILE"].SetDefault("");
  settings["PS_POINT"].SetDefault("");
  settings["TCHANNEL_ALPHA"].SetDefault(0.9);
  settings["SCHANNEL_ALPHA"].SetDefault(0.5);
  settings["CHANNEL_EPSILON"].SetDefault(0.0);
  settings["THRESHOLD_EXPONENT"].SetDefault(0.5);
  settings["ENHANCE_XS"].SetDefault(0);
}

void Phase_Space_Handler::InitParameters(const double & error) { 
  Settings& s    = Settings::GetMainSettings();
  m_thkill       = s["IB_THRESHOLD_KILL"].Get<double>();
  m_error        = s["INTEGRATION_ERROR"].Get<double>();
  m_abserror     = s["ABS_ERROR"].Get<double>();
  m_fin_opt      = s["FINISH_OPTIMIZATION"].Get<bool>();
  m_printpspoint = s["PRINT_PS_POINTS"].Get<bool>();
  if (error>0.) { m_error = error; }
}
  
void Phase_Space_Handler::CheckSinglePoint()
{
  Settings& s = Settings::GetMainSettings();
  const std::string file{ s["PS_PT_FILE"].Get<std::string>() };
  if (file!="") {
    Data_Reader read_mom(" ",";","#","=");
    read_mom.SetInputFile(file);
    read_mom.AddIgnore("Vec4D");
    read_mom.RereadInFile();
    for (size_t i(0);i<p_lab.size();++i) {
      std::vector<std::string> vec;
      if (!read_mom.VectorFromFile(vec,"p_lab["+ToString(i)+"]"))
	THROW(fatal_error,"No ps points in file");
      if (vec.front()=="-") p_lab[i]=-ToType<Vec4D>(vec.back());
      else p_lab[i]=ToType<Vec4D>(vec.front());
      msg_Debugging()<<"p_lab["<<i<<"]=Vec4D"<<p_lab[i]<<";\n";
    }
  } else if (s["PS_POINT"].IsSetExplicitly()){
    std::istringstream point(s["PS_POINT"].Get<std::string>());
    std::string vec;
    for (size_t i = 0; i<p_lab.size(); ++i) {
      std::getline(point, vec);
      if (vec.empty())
        THROW(fatal_error, "Momentum missing for calculation. ")
      auto pos1 = vec.find_first_of('(');
      auto pos2 = vec.find_first_of(')');
      vec = vec.substr(pos1, pos2-pos1+1);
      if (vec[1] == '-') p_lab[i] = -ToType<Vec4D>(vec);
      else p_lab[i] = ToType<Vec4D>(vec);
    }
  } else
    return ;
  Process_Base *proc(p_active->Process());
  proc->Trigger(p_lab);
  CalculateME(Variations_Mode::nominal_only);
  msg->SetPrecision(16);
  msg_Out()<<"// "<<proc->Name()<<"\n";
  for (size_t i(0);i<p_lab.size();++i)
    msg_Out()<<"p_lab["<<i<<"]=Vec4D"<<p_lab[i]<<";"<<std::endl;
  for (int i=0; i<proc->Size(); ++i) {
    msg_Out()<<(*proc)[i]->Name()<<" ME = "<<(*proc)[i]->Get<Single_Process>()->LastXS()
             <<", ME with PDF = "<<(*proc)[i]->Get<Single_Process>()->Last()
             <<"; // in GeV^2, incl. symfacs"<<std::endl;
    if (proc->GetSubevtList()) {
      NLO_subevtlist * subs(proc->GetSubevtList());
      for (size_t i(0);i<subs->size();++i) msg_Out()<<(*(*subs)[i]);
    }
  }
  THROW(normal_exit,"Computed ME^2");
}

void Phase_Space_Handler::TestPoint(ATOOLS::Vec4D *const p,
				    ATOOLS::Vec4D_Vector cp,ATOOLS::Flavour_Vector fl,
				    const Subprocess_Info *info,size_t &n,
				    const ATOOLS::Mass_Selector* ms)
{
  size_t nin(fl.size());
  for (size_t i(0);i<nin;++i) msg_Debugging()<<fl[i]<<" ";
  msg_Debugging()<<"->";
  fl.resize(nin+info->m_ps.size());
  cp.resize(nin+info->m_ps.size());
  for (size_t i(0);i<info->m_ps.size();++i) {
    fl[nin+i]=info->m_ps[i].m_fl;
    msg_Debugging()<<" "<<fl[nin+i];
  }
  msg_Debugging()<<" {\n";
  if (info->m_ps.size()==1) {
    for (size_t i(0);i<nin;++i) cp.back()+=cp[i];
  }
  else {
    Single_Channel * TestCh = new Rambo(nin,info->m_ps.size(),&fl.front(),ms);
    TestCh->GeneratePoint(&cp.front(),(Cut_Data*)(NULL));
    delete TestCh;
    if (nin==1) {
      Poincare cms(cp.front());
      for (size_t i(1);i<cp.size();++i) cms.BoostBack(cp[i]);
    }
  }
  for (size_t i(0);i<info->m_ps.size();++i) {
    msg_Indent();
    if (info->m_ps[i].m_ps.empty()) {
      msg_Debugging()<<"p["<<n<<"] = "<<cp[nin+i]<<", m = "
		     <<sqrt(dabs(cp[nin+i].Abs2()))<<" ("<<fl[nin+i]<<")\n";
      p[n++]=cp[nin+i];
    }
    else {
      msg_Debugging()<<"P["<<nin+i<<"] = "<<cp[nin+i]<<", m = "
		     <<sqrt(dabs(cp[nin+i].Abs2()))<<" ("<<fl[nin+i]<<")\n";
      Vec4D_Vector ncp(1,cp[nin+i]);
      Flavour_Vector nfl(1,info->m_ps[i].m_fl);
      TestPoint(p,ncp,nfl,&info->m_ps[i],n,ms);
    }
  }
  msg_Debugging()<<"}\n";
}

void Phase_Space_Handler::TestPoint(ATOOLS::Vec4D *const p,
				    const Process_Info *info,
				    const ATOOLS::Mass_Selector* ms,
				    const int mode)
{
  DEBUG_FUNC(mode);
  Flavour_Vector fl_i(info->m_ii.GetExternal());
  Vec4D_Vector cp(fl_i.size());
  if (fl_i.size()==1) {
    double m(0.0);
    for (size_t j(0);j<fl_i[0].Size();++j) m+=ms->Mass(fl_i[0][j]);
    p[0]=cp[0]=Vec4D(m/fl_i[0].Size(),0.0,0.0,0.0);
    msg_Debugging()<<"p[0] = "<<p[0]<<"\n";
  }
  else {
    double m[2]={fl_i[0].Mass(),fl_i[1].Mass()};
    double E=rpa->gen.Ecms();
    if (info->m_fi.m_ps.size()==1 &&
	info->m_fi.m_ps[0].m_ps.empty()) {
      E=0.0;
      Flavour dfl(info->m_fi.m_ps.front().m_fl);
      for (size_t j(0);j<dfl.Size();++j) E+=ms->Mass(dfl[j]);
      E/=dfl.Size();
    }
    if (E<m[0]+m[1]) return;
    double x=1.0/2.0+(m[0]*m[0]-m[1]*m[1])/(2.0*E*E);
    p[0]=cp[0]=Vec4D(x*E,0.0,0.0,sqrt(sqr(x*E)-m[0]*m[0]));
    p[1]=cp[1]=Vec4D((1.0-x)*E,Vec3D(-p[0]));
    msg_Debugging()<<"p[0] = "<<p[0]<<"\np[1] = "<<p[1]<<"\n";
  }
  
  unsigned int osd_counter=0;
  for (size_t i=0;i<info->m_fi.GetDecayInfos().size();i++)
    if (info->m_fi.GetDecayInfos()[i]->m_osd) osd_counter++;
    
  if (osd_counter==info->m_fi.GetDecayInfos().size() || mode==1) {
    size_t n(fl_i.size());
    TestPoint(p,cp,fl_i,&info->m_fi,n,ms);
  }
  else {
    Flavour_Vector fl_f(info->m_fi.GetExternal());
    Flavour_Vector fl_tot(fl_i);
    fl_tot.insert(fl_tot.end(),fl_f.begin(),fl_f.end());
    //
    Single_Channel * TestCh = new Rambo(fl_i.size(),fl_f.size(),&fl_tot.front(),ms);
    TestCh->GeneratePoint(p,(Cut_Data*)(NULL));
    //
    delete TestCh;
  }
}

void Phase_Space_Handler::TestPoint(ATOOLS::Vec4D *const p,
				    const size_t &nin,const size_t &nout,
				    const Flavour_Vector &flavs,
				    const ATOOLS::Mass_Selector* ms)
{
  if (nin==1) {
    p[0]=Vec4D(flavs[0].Mass(),0.0,0.0,0.0);
    if (nout==1) { 
      p[1]=p[0]; 
      return;
    }
  }
  else {
    double m[2]={flavs[0].Mass(),flavs[1].Mass()};
    double E=0.5*rpa->gen.Ecms();
    if (E<m[0]+m[1]) return;
    double x=1.0/2.0+(m[0]*m[0]-m[1]*m[1])/(2.0*E*E);
    p[0]=Vec4D(x*E,0.0,0.0,sqrt(sqr(x*E)-m[0]*m[0]));
    p[1]=Vec4D((1.0-x)*E,Vec3D(-p[0]));
  }
  Single_Channel * TestCh = new Rambo(nin,nout,&flavs.front(),ms);
  TestCh->GeneratePoint(p,(Cut_Data*)(NULL));
  delete TestCh;
}
