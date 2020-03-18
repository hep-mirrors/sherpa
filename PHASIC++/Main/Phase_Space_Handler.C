#include "PHASIC++/Main/Phase_Space_Handler.H"

#include "PHASIC++/Main/Phase_Space_Integrator.H"
#include "PHASIC++/Main/Channel_Creator.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Channels/FSR_Channels.H"
#include "PHASIC++/Channels/ISR_Channels.H"
#include "PHASIC++/Channels/Beam_Channels.H"
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

Phase_Space_Handler::Phase_Space_Handler(Process_Integrator *proc,double error): 
  m_name(proc->Process()->Name()), p_process(proc), p_active(proc),
  p_integrator(NULL), 
  m_psvariationweights(false), 
  p_beamhandler(proc->Beam()), p_isrhandler(proc->ISR()), p_flavours(proc->Process()->Flavours()),
  m_nin(proc->NIn()), m_nout(proc->NOut()), m_nvec(m_nin+m_nout),
  m_dmode(1), m_initialized(0),
  m_sintegrator(0), m_killedpoints(0),
  m_printpspoint(false)
{
  RegisterDefaults();
  InitParameters(error);
  p_process->SetPSHandler(this);
  
  p_lab.resize(m_nvec);

  
  if (CreateIntegrators()) {
    m_pspoint.Init(this);
    m_psenhance.Init(this);
  }
}

Phase_Space_Handler::~Phase_Space_Handler()
{
  delete p_integrator;
}

bool Phase_Space_Handler::CreateIntegrators() {
  Channel_Creator channelcreator(this);
  if (!channelcreator()) THROW(fatal_error,"Creation of integration channels failed.");
  return true;
}
  
bool Phase_Space_Handler::InitIncoming() 
{
  return m_pspoint.MakeIncoming();
  //if (m_nin>1) {
  //  m_smin=ATOOLS::Max(sqr(p_process->ISRThreshold()),p_cuts->Smin());
  //}
  //m_initialized=1;
}

double Phase_Space_Handler::Integrate() 
{
  CheckSinglePoint();
  if (p_process->Points()>0 &&
      (p_process->TotalError()<dabs(m_error*p_process->TotalXS()) ||
       p_process->TotalError()<m_abserror)) 
    return p_process->TotalXS()*rpa->Picobarn();
  p_integrator = new Phase_Space_Integrator(this);
  if (!InitIncoming()) return 0;
  /*
    Phase_Space_Point() -> Print()
  */
  m_dmode=0;
  double res(0.0);
  if (m_nin==2) res=p_integrator->Calculate(m_error,m_abserror,m_fin_opt);
  if (m_nin==1) res=p_integrator->CalculateDecay(m_error);
  m_dmode=1;
  return res;
}

double Phase_Space_Handler::Differential(Process_Integrator *const process,
					 const psmode::code mode) 
{
  //msg_Out()<<"-------------------------------------------------------------\n"
  //	   <<"-------------------------------------------------------------\n"
  //	   <<"-------------------------------------------------------------\n"
  //	   <<METHOD<<" for ISR reset mode = "<<(mode&psmode::no_lim_isr)<<", "
  //	   <<(mode&psmode::no_gen_isr)<<", "
  //	   <<"beam handler on = "<<p_beamhandler->On()<<"\n";
  m_cmode  = mode;
  p_active = process;
  m_result = 0.0;
  // check for failure to generate a meaningful phase space point
  if (!process->Process()->GeneratePoint() ||
      !m_pspoint(process,m_cmode))                          return m_result;
  //msg_Out()<<" ### "<<p_lab[0]<<" + "<<p_lab[1]<<"\n"
  //	   <<"  -> "<<p_lab[2]<<" + "<<p_lab[3]<<".\n";
  for (size_t i(0);i<p_lab.size();++i) if (p_lab[i].Nan())  return m_result;
  if (process->Process()->Trigger(p_lab)) {
    if (!p_active->Process()->Selector()->Pass())           return m_result;
    m_result  = (m_meweight = CalculateME());
    m_result *= (m_psweight = CalculatePS());
    m_result *= (m_ISsymmetryfactor = m_pspoint.ISSymmetryFactor());
    p_lab=process->Momenta();
    if (m_printpspoint || msg_LevelIsDebugging()) PrintIntermediate();
    ManageWeights(m_psweight*m_ISsymmetryfactor);
  }
  if (!CheckStability()) { m_result = 0.; return 0.; }
  m_enhance = m_psenhance.Factor(p_process->Process(),p_process->TotalXS());
  return m_result*m_enhance;
}


void Phase_Space_Handler::PrintIntermediate() {
  size_t precision(msg->Out().precision());
  msg->SetPrecision(15);
  msg_Out()<<"==========================================================\n"
	   <<p_active->Process()->Name()
	   <<"  ME = "<<m_result<<" ,  PS = "<<m_psweight<<"  ->  "
	   <<m_result*m_psweight<<std::endl;
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
  if (m_psvariationweights && m_psvariation()) *m_psvariation() *= factor;
  ME_Weight_Info* wgtinfo=p_active->Process()->GetMEwgtinfo();
  if (wgtinfo) { (*wgtinfo) *= factor; }
  NLO_subevtlist* nlos=p_active->Process()->GetSubevtList();
  if (nlos) { (*nlos) *= factor; (*nlos).MultMEwgt(factor); }
}

bool Phase_Space_Handler::CheckStability() {
  // meaningful result - no problems down the line
  if (p_active->TotalXS() &&
      dabs(m_result/p_active->TotalXS())>dabs(m_thkill)) {
    if (m_thkill<0.0) {
      msg_Info()<<METHOD<<"(): Skip point in '"<<p_active->Process()->Name()<<"', "
		<<"weight = "<<m_result*rpa->Picobarn()<<", thkill = "<<m_thkill<<",\n"
		<<"   totalxs = "<<p_active->TotalXS()<<", result = "<<m_result<<".\n";
      return false;
    }
    // outputb tricky phase space point for further analysis, if necessary
    ATOOLS::MakeDir("stability");
    std::ofstream sf(("stability/"+p_active->Process()->Name()+
		      "_"+rpa->gen.Variable("RNG_SEED")).c_str(),
		     std::ios_base::app);
    sf.precision(16);
    sf<<"(P"<<m_killedpoints<<"){ # w = "
      <<m_result<<", ME = "<<m_result/m_psweight<<", PS = "<<m_psweight<<"\n";
    for (size_t i(0);i<p_lab.size();++i) sf<<"  p_lab["<<i<<"]=Vec4D"<<p_lab[i]<<";\n";
    sf<<"}(P"<<m_killedpoints<<");\n";
    ++m_killedpoints;
    ManageWeights(0.);
    return false;
  }
  return true;
}


Weight_Info *Phase_Space_Handler::OneEvent(Process_Base *const proc,int mode)
{
  if (!m_initialized) InitIncoming();
  if (proc==NULL) THROW(fatal_error,"No process.");
  Process_Integrator *cur(proc->Integrator());
  p_isrhandler->SetRunMode(1);
  double value=Differential(cur,(psmode::code)mode);
  bool zero(value==0.0);
  if (zero && m_psvariationweights)
    for (size_t i(0);i<m_psvariation()->GetNumberOfVariations();++i)
      if (m_psvariation()->GetVariationWeightAt(i)) zero=false;
  if (zero || IsBad(value)) return NULL;
  cur->SetMomenta(p_lab);
  int fl1(0), fl2(0);
  double x1(0.0), x2(0.0), xf1(0.0), xf2(0.0), mu12(0.0), mu22(0.0), dxs(0.0);
  dxs=p_active->Process()->Get<PHASIC::Single_Process>()->LastXS();
  fl1=(long int)p_active->Process()->Flavours()[0];
  fl2=(long int)p_active->Process()->Flavours()[1];
  x1=p_isrhandler->X1();
  x2=p_isrhandler->X2();
  xf1=p_isrhandler->XF1(0);
  xf2=p_isrhandler->XF2(0);
  mu12=p_isrhandler->MuF2(0);
  mu22=p_isrhandler->MuF2(1);
  return new Weight_Info(value,dxs,1.0,fl1,fl2,x1,x2,xf1,xf2,mu12,mu22);
}





void Phase_Space_Handler::AddPoint(const double _value)
{
  p_process->AddPoint(_value);
  double value(_value);
  if (p_process->TotalXS()==0.0) value=_value?1.0:0.0;
  if (value!=0.0) {
    double enhancexs = m_psenhance();
    m_pspoint.AddPoint(value*enhancexs);
    m_psenhance.AddPoint(value*enhancexs,p_process->Process());
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
  msg_Info()<<"Read in channels from directory : "<<pID<<std::endl;
  m_psenhance.ReadIn(pID);
  bool okay = true;
  Data_Reader reader;
  reader.SetInputPath(pID+"/");
  reader.SetInputFile("Statistics.dat");
  std::vector<std::vector<double> > stats;
  if (reader.MatrixFromFile(stats,"")) m_stats=stats;
  return okay;
}

bool Phase_Space_Handler::UpdateIntegrators()
{
  if (!m_sintegrator || m_nout==1) return false;
  double error=Process()->TotalVar()/Process()->TotalResult();
  msg_Info()<<om::blue
	    <<Process()->TotalResult()*rpa->Picobarn()
	    <<" pb"<<om::reset<<" +- ( "<<om::red
	    <<Process()->TotalVar()*rpa->Picobarn()
	    <<" pb = "<<error*100<<" %"<<om::reset<<" ) "
	    <<FSRIntegrator()->ValidN()<<" ( "
	    <<(FSRIntegrator()->ValidN()*1000/FSRIntegrator()->N())/10.0<<" % ) "<<std::endl;
  p_process->Process()->UpdateIntegrator(this);
  return true;
}

void Phase_Space_Handler::RegisterDefaults() const
{
  Settings& settings = Settings::GetMainSettings();
  settings["IB_THRESHOLD_KILL"].SetDefault(-1.0e12);
  settings["ERROR"].SetDefault(0.01);
  settings["INTEGRATION_ERROR"].SetDefault(settings["ERROR"].Get<double>());
  settings["ABS_ERROR"].SetDefault(0.0);
  settings["MAX_TRIALS"].SetDefault(1000000);
  settings["FINISH_OPTIMIZATION"].SetDefault(true);
  settings["PRINT_PS_POINTS"].SetDefault(false);
  settings["INT_MINALPHA"].SetDefault(0.0);
  settings["PS_PT_FILE"].SetDefault("");
  settings["TCHANNEL_ALPHA"].SetDefault(0.9);
  settings["SCHANNEL_ALPHA"].SetDefault(0.75);
  settings["CHANNEL_EPSILON"].SetDefault(0.0);
  settings["THRESHOLD_EPSILON"].SetDefault(1.5);
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
  msg_Out()<<METHOD<<"\n";
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
    Process_Base *proc(p_active->Process());
    proc->Trigger(p_lab);
    CalculateME();
    msg->SetPrecision(16);
    msg_Out()<<"// "<<proc->Name()<<"\n";
    for (size_t i(0);i<p_lab.size();++i)
      msg_Out()<<"p_lab["<<i<<"]=Vec4D"<<p_lab[i]<<";"<<std::endl;
    if (proc->Get<Single_Process>()) {
      msg_Out()<<"double ME = "<<proc->Get<Single_Process>()->LastXS()
	       <<"; // in GeV^2, incl. symfacs"<<std::endl;
      if (proc->GetSubevtList()) {
	NLO_subevtlist * subs(proc->GetSubevtList());
	for (size_t i(0);i<subs->size();++i) msg_Out()<<(*(*subs)[i]);
      }
    }
    THROW(normal_exit,"Computed ME^2");
  }
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


