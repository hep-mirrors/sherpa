#include "PHASIC++/Main/Process_Integrator.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "PHASIC++/Channels/Extra_Emission_Generator.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Smart_Pointer.C"

using namespace PHASIC;
using namespace ATOOLS;

namespace ATOOLS { template class SP(Process_Integrator); }

static int s_whbins(1000);
static int s_omitnlosuffix(0);
 
Process_Integrator::Process_Integrator(Process_Base *const proc):
  p_proc(proc), p_pshandler(NULL),
  p_beamhandler(NULL), p_isrhandler(NULL), p_eeg(NULL),
  m_nin(0), m_nout(0), m_corenout(0), m_smode(1),
  m_threshold(0.), m_overflow(0.), m_enhancefac(1.0), m_maxeps(0.0),
  m_xinfo(std::vector<double>(4)), m_n(0), m_rbn(0), m_itmin(0), 
  m_max(0.), m_rbmax(0.0), m_rbsum(0.0), m_rbsum2(0.0), m_totalxs(0.), 
  m_totalsum (0.), m_totalsumsqr(0.), m_totalerr(0.), m_ssum(0.), 
  m_ssumsqr(0.), m_smax(0.), m_ssigma2(0.), m_wmin(0.), m_vmean(0.), m_sn(0), 
  m_svn(0), m_son(1), m_swaped(false), m_writeout(false),
  m_efunc("1"), p_whisto(NULL), p_colint(NULL), p_helint(NULL)
{
  m_colorscheme=cls::sum;
  m_helicityscheme=hls::sum;
}

bool Process_Integrator::Initialize
(BEAM::Beam_Spectra_Handler *const beamhandler,
 PDF::ISR_Handler *const isrhandler)
{
  m_nin=p_proc->NIn();
  m_corenout=m_nout=p_proc->NOut();
  p_momenta.resize(m_nin+m_nout);
  p_beamhandler=beamhandler;
  p_isrhandler=isrhandler;
  static bool minit(false);
  if (!minit) {
    Data_Reader read(" ",";","!","=");
    int smode;
    if (read.ReadFromFile(smode,"IB_SMODE")) {
      m_smode=smode;
      msg_Info()<<METHOD<<"(): Set sum mode = "<<m_smode<<".\n";
    }
    if (!read.ReadFromFile(s_whbins,"IB_WHBINS")) s_whbins=1000;
    else msg_Info()<<METHOD<<"(): Set weight histo bin number = "<<s_whbins<<".\n";
    if (read.ReadFromFile(s_omitnlosuffix,"RESULT_OMIT_NLO_SUFFIX")) {
      if (s_omitnlosuffix) 
	msg_Info()<<METHOD<<"(): store/read results without Process_Info suffix.\n";
    }
    else s_omitnlosuffix=0;
    minit=true;
  }
  return true;
}

Process_Integrator::~Process_Integrator()
{
  if (p_eeg!=NULL) delete p_eeg;
  if (p_whisto!=NULL) delete p_whisto;
}

double Process_Integrator::Sigma2()
{ 
  Process_Integrator *p(p_proc->Parent()->Integrator());
  if (m_sn!=p->m_sn) {
    msg_Error()<<METHOD<<"(): Inconsistent summation for '"
	       <<p_proc->Name()<<"' \\in '"<<p->Process()->Name()
	       <<"', m_sn = "<<m_sn<<" vs. p->m_sn = "
	       <<p->m_sn<<"."<<std::endl;
    if (msg_LevelIsTracking()) exh->GenerateStackTrace(std::cout);
  }
  if (m_sn<2) return 0.0;
  return sqr(p->m_ssum/m_sn)/
    ((p->m_ssumsqr/m_sn-sqr(p->m_ssum/m_sn))/(m_sn-1));
}

double Process_Integrator::TotalSigma2()
{ 
  return m_ssigma2+Sigma2();
}

double Process_Integrator::TotalResult()
{ 
  if (m_ssigma2==0.0) return m_ssum/m_sn; 
  switch (m_smode) {
  case 2: {
    if (m_sn<2) return 1.0/(m_totalsum/m_ssigma2); 
    double s2(Sigma2());
    return 1.0/((m_totalsum+s2/m_ssum*m_sn)/(m_ssigma2+s2));
  }
  default: {
    if (m_sn<2) return m_totalsum/m_ssigma2; 
    double s2(Sigma2());
    return (m_totalsum+s2*m_ssum/m_sn)/(m_ssigma2+s2);
  }
  }
  THROW(fatal_error,"Invalid summation mode");
  return 0.0;
}

double Process_Integrator::TotalVar() 
{
  if (m_nin==1 && m_nout==2) return 0.;
  switch (m_smode) {
  case 2: {
    double s2(m_totalsumsqr);
    if (m_sn>1) {
      double vij2(sqr(m_ssum/m_sn)*(m_sn-1)/
		  (m_ssumsqr/m_sn-sqr(m_ssum/m_sn)));
      s2+=sqr(Sigma2())/vij2*sqr(m_sn/m_ssum);
    }
    return sqr(TotalResult())*sqrt(s2)/TotalSigma2();
  }
  default: {
    double s2(m_totalsumsqr);
    if (m_sn>1) {
      double vij2(sqr(m_ssum/m_sn)*(m_sn-1)/
		  (m_ssumsqr/m_sn-sqr(m_ssum/m_sn)));
      s2+=sqr(Sigma2())/vij2*sqr(m_ssum/m_sn);
    }
    return sqrt(s2)/TotalSigma2();
  }
  }
  THROW(fatal_error,"Invalid summation mode");
  return 0.0;
}

void Process_Integrator::OptimizeSubResult(const double &s2)
{
  if (m_smode>0 && s2>m_wmin) {
    double vij2(sqr(m_ssum/m_sn)*(m_sn-1)/
		(m_ssumsqr/m_sn-sqr(m_ssum/m_sn)));
    m_ssigma2+=s2; 
    switch (m_smode) {
    case 1:
      m_totalsum+=s2*m_ssum/m_sn;
      m_totalsumsqr+=sqr(s2)/vij2*sqr(m_ssum/m_sn);
      break;
    case 2:
      m_totalsum+=s2*m_sn/m_ssum;
      m_totalsumsqr+=sqr(s2)/vij2*sqr(m_sn/m_ssum);
      break;
    }
    if (s2/m_son>m_wmin) m_wmin=s2/m_son;
    m_ssum=m_ssumsqr=0.0;
    m_son=m_svn=m_sn=0;
  }
  m_son++;
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->OptimizeSubResult(s2);
}

double Process_Integrator::RemainTimeFactor(double maxerr) 
{
  if (m_sn<1000) return 0.;
  return (sqr(1./maxerr)-m_ssigma2)/Sigma2();
}

void Process_Integrator::SetMomenta(const Vec4D_Vector &p) 
{ 
  p_momenta=p;
  if (p_proc->Selected()!=p_proc) 
    p_proc->Selected()->Integrator()->p_momenta=p;
}

void Process_Integrator::InitWeightHistogram() 
{
  if (p_whisto) {
    delete p_whisto; 
    p_whisto=0;
  }
  if (!m_writeout) return;

  double av(TotalResult());
  if (!av>0.) {
    msg_Error()<<"Process_Integrator::InitWeightHistogram(): "
		       <<"No valid result: "<<av<<std::endl;
    return;
  }
  if (av<.3) av/=10.;
  av = exp(log(10.)*int(log(av)/log(10.)+0.5));
  p_whisto = new Histogram(10,av*1.e-4,av*1.e6,s_whbins);
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->InitWeightHistogram();
}

bool Process_Integrator::ReadInXSecs(const std::string &path)
{
  std::string fname(p_proc->Name());
  if (s_omitnlosuffix) {
    size_t pos=fname.find("EW");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    pos=fname.find("QCD");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
  }
  std::string name;
  std::ifstream from((path+"/"+fname).c_str());
  if (!from.good()) return false;
  from.precision(12);
  from>>name>>m_totalxs>>m_max>>m_totalerr>>m_totalsum>>m_totalsumsqr>>m_n
      >>m_ssum>>m_ssumsqr>>m_ssigma2>>m_sn>>m_wmin>>m_son>>m_vmean>>m_svn;
  if (name!=fname) THROW(fatal_error,"Corrupted results file");
  msg_Tracking()<<"Found result: xs for "<<name<<" : "
		<<m_totalxs*rpa.Picobarn()<<" pb"
		<<" +- ( "<<m_totalerr*rpa.Picobarn()<<" pb = "
		<<m_totalerr/m_totalxs*100.<<" % ) max: "
		<<m_max*rpa.Picobarn()<<std::endl;
  bool res(true);
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      if (!(*p_proc)[i]->Integrator()->ReadInXSecs(path)) res=false;
  SetMax(m_max);
  return res;
}

void Process_Integrator::ReadInHistogram(std::string dir)
{
  std::string filename = dir+"/"+p_proc->Name();
  std::ifstream from(filename.c_str());
  if (!from.good()) return;
  if (p_whisto) delete p_whisto; 
  p_whisto = new Histogram(filename);	
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->ReadInHistogram(dir);
}

void Process_Integrator::WriteOutXSecs(const std::string &path)
{
  std::string fname(p_proc->Name());
  if (s_omitnlosuffix) {
    size_t pos=fname.find("EW");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    pos=fname.find("QCD");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
  }
  std::ofstream outfile((path+"/"+fname).c_str());
  if (outfile) m_writeout=1;
  outfile.precision(12);
  outfile<<fname<<"  "<<m_totalxs<<"  "<<m_max<<"  "
	 <<m_totalerr<<" "<<m_totalsum<<" "<<m_totalsumsqr<<" "
	 <<m_n<<" "<<m_ssum<<" "<<m_ssumsqr<<" "<<m_ssigma2<<" "
	 <<m_sn<<" "<<m_wmin<<" "<<m_son<<" "
	 <<m_vmean<<" "<<m_svn<<std::endl; 
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->WriteOutXSecs(path);
}

void Process_Integrator::WriteOutHistogram(std::string dir)
{
  if (p_whisto) p_whisto->Output(dir+"/"+p_proc->Name());	
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->WriteOutHistogram(dir);
}

void Process_Integrator::SetTotal(const int mode)  
{ 
  if (p_proc->IsGroup()) {
    m_selectionweight=0.0;
    m_max=0.0;
    msg_Indent();
    for (size_t i(0);i<p_proc->Size();++i) {
      (*p_proc)[i]->Integrator()->SetTotal(msg_LevelIsTracking());
      m_max+=(*p_proc)[i]->Integrator()->Max();
      m_selectionweight+=(*p_proc)[i]->Integrator()->SelectionWeight();
    }
  }
  m_totalxs=TotalResult();
  if (!p_proc->IsGroup()) m_selectionweight=dabs(m_totalxs);
  m_totalerr=TotalVar();
  if (mode) {
    if (p_proc->NIn()==1) {
      msg_Info()<<om::bold<<p_proc->Name()<<om::reset<<" : "<<om::blue<<om::bold
                <<m_totalxs<<" GeV"<<om::reset<<" +- ( "<<om::red
                <<m_totalerr<<" GeV = "<<m_totalerr/m_totalxs*100.
                <<" %"<<om::reset<<" ) "<<om::bold<<" exp. eff: "
                <<om::red<<(100.*m_totalxs/m_max)<<" %"<<om::reset<<std::endl;
    }
    else {
      msg_Info()<<om::bold<<p_proc->Name()<<om::reset<<" : "<<om::blue<<om::bold
                <<m_totalxs*rpa.Picobarn()<<" pb"<<om::reset<<" +- ( "<<om::red
                <<m_totalerr*rpa.Picobarn()<<" pb = "<<m_totalerr/m_totalxs*100.
                <<" %"<<om::reset<<" ) "<<om::bold<<" exp. eff: "
                <<om::red<<(100.*m_totalxs/m_max)<<" %"<<om::reset<<std::endl;
    }
  }
}

double Process_Integrator::GetMaxEps(double epsilon)
{
  if (!p_whisto) return m_max;
  double pxs = TotalResult()*epsilon*p_whisto->Fills();
  double cutxs = 0.;
  double cnt = 0.;

  for (int i=p_whisto->Nbin()+1;i>0;i--) {
    cutxs+= p_whisto->Value(i) * 
      exp(log(10.)*(p_whisto->Xmin()+(i-0.5)*p_whisto->BinSize()));
    cnt+= p_whisto->Value(i);
    double twmax = exp(log(10.)*(p_whisto->Xmin()+(i-1)*p_whisto->BinSize()));
    if (cutxs-cnt*twmax>pxs) {
      return ATOOLS::Min(exp(log(10.)*(p_whisto->Xmin()+i*p_whisto->BinSize())),m_max);
    }
  }
  return m_max;
}

void Process_Integrator::SetUpEnhance() 
{
  m_selectionweight=SelectionWeight()*m_enhancefac;
  if (m_enhancefac<0.0) m_selectionweight=p_proc->Size()*dabs(m_enhancefac);
  if (m_maxeps>0.0) {
    double max(GetMaxEps(m_maxeps));
    msg_Info()<<"Max reduction for '"<<p_proc->Name()<<"': "<<Max()/max
	      <<" ( \\epsilon = "<<m_maxeps<<")"<<std::endl;
    SetMax(max);
  }
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->SetUpEnhance();
}

void Process_Integrator::SetEnhanceFunction(const std::string &efunc)
{ 
  m_efunc=efunc; 
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->SetEnhanceFunction(efunc);
}

void Process_Integrator::SetEnhanceFactor(const double &efac)
{ 
  m_enhancefac=efac; 
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->SetEnhanceFactor(efac);
}

void Process_Integrator::SetMaxEpsilon(const double &maxeps)
{ 
  m_maxeps=maxeps; 
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->SetMaxEpsilon(maxeps);
}

void Process_Integrator::AddPoint(const double value) 
{
  m_n++;
  m_sn++;
  if (value!=0.0) ++m_svn;
  m_ssum    += value;
  m_ssumsqr += value*value;
  if (value>m_max)  m_max  = value;
  if (value>m_smax) m_smax = value;
  if (p_whisto) {
    if(value!=0.) p_whisto->Insert(value);
    else p_whisto->SetFills(p_whisto->Fills()+1);
  }
  if (!p_proc->IsGroup()) {
    p_proc->AddPoint(value);
  }
  else {
    if (p_proc->Last()==0.0 || value==0.0)
      for (size_t i(0);i<p_proc->Size();++i)
	(*p_proc)[i]->Integrator()->AddPoint(0.0);
    else
      for (size_t i(0);i<p_proc->Size();++i)
	(*p_proc)[i]->Integrator()->
	  AddPoint(value*(*p_proc)[i]->Last()/p_proc->Last());
  }
}

void Process_Integrator::AddRBPoint(const double rb) 
{
  if (p_proc->IsGroup()) THROW(fatal_error,"Invalid call");
  ++m_rbn;
  m_rbsum+=rb;
  m_rbsum2+=rb*rb;
  m_rbmax=ATOOLS::Max(m_rbmax,rb);
}

void Process_Integrator::SetMax(const double max,const int flag) 
{
  if (flag==1) m_max=max;
  if (!p_proc->IsGroup()) return;
  if (flag==1) {
    for (size_t i(0);i<p_proc->Size();++i) 
      (*p_proc)[i]->Integrator()->SetMax(max/(double)p_proc->Size(),flag);
    return;
  }
  double sum(0.0);
  m_max=0.0;
  for (size_t i(0);i<p_proc->Size();++i) {
    sum+=(*p_proc)[i]->Integrator()->TotalXS();
    m_max+=(*p_proc)[i]->Integrator()->Max();
  }
  if (m_totalxs!=0.) {
    if (!ATOOLS::IsEqual(sum,m_totalxs,1e-11)) {
      msg_Error().precision(12);
      msg_Error()<<METHOD<<"(): Summation does not agree for '"
		 <<p_proc->Name()<<".\n  sum = "<<sum
		 <<" vs. total = "<<m_totalxs<<" ("
		 <<((sum-m_totalxs)/m_totalxs)<<")"<<std::endl;
      msg_Error().precision(6);
    }
    m_totalxs=sum;
  }
} 

void Process_Integrator::InitEEG()
{
  if (p_eeg) {
    delete p_eeg;
    p_eeg=NULL;
  }
  if (p_proc->Info().m_fi.m_nloqcdtype!=nlo_type::lo){
    p_eeg = new Extra_Emission_Generator(this,true);
  }
  if (p_proc->Info().m_fi.m_nloewtype!=nlo_type::lo){
    p_eeg = new Extra_Emission_Generator(this,false);
  }
}

void Process_Integrator::Reset()
{
  m_n=0;
  m_vmean=m_totalxs=m_totalsum=m_totalsumsqr=m_totalerr=0.0;
  m_smax=m_max=m_wmin=m_ssigma2=m_ssumsqr=m_ssum=0.0;
  m_svn=m_sn=0;
  m_son=1;
  m_vsmax.clear(); 
  m_vsn.clear();   
  m_vsum.clear(); 
  m_vsvn.clear();   
  if (p_eeg) ResetRB();
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i) 
      (*p_proc)[i]->Integrator()->Reset();
}

void Process_Integrator::ResetRB()
{
  m_rbn=0;
  m_rbmax=m_rbsum2=m_rbsum=0.0;
}

void Process_Integrator::ResetMax(int flag) 
{
  if (p_proc->IsGroup()) {
    m_max=0.0;
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->ResetMax(flag);
    return;
  }
  if (flag==3) {
    m_vsmax.clear();
    m_vsn.clear();
    m_vsum.clear();
    m_vsvn.clear();
    m_vmean=m_max=0.0;
    return;
  }
  if (flag==0) {
    if (m_vsmax.size()>1) {
      m_vsmax.erase(m_vsmax.begin());
      m_vsn.erase(m_vsn.begin());
      m_vsum.erase(m_vsum.begin());
      m_vsvn.erase(m_vsvn.begin());
    }
    if (m_vsmax.empty()) {
      m_vsmax.push_back(m_max);
      m_vsn.push_back(m_n);
      m_vsum.push_back(m_ssum);
      m_vsvn.push_back(m_svn);
    }
    m_vsmax.back() = ATOOLS::Max(m_smax,m_vsmax.back());
    m_vsn.back()   = m_n;
    m_vsum.back()  = m_ssum;
    m_vsvn.back()  = m_svn;
  }
  else {
    if (flag==2 && m_vsmax.size()==4) {
      m_vsmax.erase(m_vsmax.begin());
      m_vsn.erase(m_vsn.begin());
      m_vsum.erase(m_vsum.begin());
      m_vsvn.erase(m_vsvn.begin());
    }
    m_vsmax.push_back(m_smax);
    m_vsn.push_back(m_n);
    m_vsum.push_back(m_ssum);
    m_vsvn.push_back(m_svn);
    if (flag==2) m_smax = 0.;
  }
  m_max=0.0;
  double sum(0.0), svn(0.0);
  for (size_t i=0;i<m_vsmax.size();i++) {
    m_max=ATOOLS::Max(m_max,m_vsmax[i]);
    sum+=m_vsum[i];
    svn+=m_vsvn[i];
  }
  m_vmean=sum/svn;
} 

void Process_Integrator::SetPSHandler(const SP(Phase_Space_Handler) &pshandler)
{
  p_pshandler=pshandler;
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->SetPSHandler(pshandler);
} 

void Process_Integrator::Optimize()
{
  if (p_eeg) p_eeg->Optimize();
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->Optimize();
} 

void Process_Integrator::EndOptimize()
{
  if (p_eeg) p_eeg->EndOptimize();
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->EndOptimize();
} 

void Process_Integrator::WriteOutEEG(const std::string &pid)
{
  if (p_eeg) p_eeg->WriteOut(pid);
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->WriteOutEEG(pid);
}

bool Process_Integrator::ReadInEEG(const std::string &pid)
{
  if (p_eeg) p_eeg->ReadIn(pid);
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->ReadInEEG(pid);
  return true;
}

void Process_Integrator::OptimizeResult()
{
  OptimizeSubResult(Sigma2());
} 

void Process_Integrator::SetMomenta()   
{ 
  p_momenta=p_pshandler->LabPoint();
}

double Process_Integrator::TriggerEfficiency()
{
  if (p_pshandler!=NULL) {
    Multi_Channel *fsr=p_pshandler->FSRIntegrator();
    return fsr->ValidN()/(double)fsr->N();
  }
  return 1.0;
}

void Process_Integrator::SetISRThreshold(const double threshold) 
{
  m_threshold=threshold;
}

void Process_Integrator::SwapInOrder() 
{
  if (!m_swaped) {
    p_proc->SwapInOrder();
    std::swap<Vec4D>(p_momenta[0],p_momenta[1]);
    for (size_t i(0);i<p_momenta.size();++i) 
      p_momenta[i]=Vec4D(p_momenta[i][0],-p_momenta[i]);
    m_swaped=true;
  }
}

void Process_Integrator::RestoreInOrder() 
{
  if (m_swaped) {
    p_proc->SwapInOrder();
    std::swap<Vec4D>(p_momenta[0],p_momenta[1]);
    for (size_t i(0);i<p_momenta.size();++i) 
      p_momenta[i]=Vec4D(p_momenta[i][0],-p_momenta[i]);
    m_swaped=false;
  }
}

void Process_Integrator::StoreResults(const int mode)
{
  if (m_resultpath.length()==0) return;
  if (m_totalxs!=0.0 && mode==0) return;
  SetTotal(0); 
  std::string fname(p_proc->Name());
  if (s_omitnlosuffix) {
    size_t pos=fname.find("EW");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    pos=fname.find("QCD");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
  }
  MakeDir(m_resultpath+"/XS_"+fname,0); 
  MakeDir(m_resultpath+"/WD_"+p_proc->Name(),0); 
  MakeDir(m_resultpath+"/MC_"+fname,0); 
  MakeDir(m_resultpath+"/RB_"+p_proc->Name(),0); 
  WriteOutXSecs(m_resultpath+"/XS_"+fname);
  WriteOutHistogram(m_resultpath+"/WD_"+p_proc->Name());
  p_pshandler->WriteOut(m_resultpath+"/MC_"+fname);
  WriteOutRB(m_resultpath+"/RB_"+p_proc->Name());
}

void Process_Integrator::ReadResults()
{
  if (m_resultpath.length()==0) return;
  std::string fname(p_proc->Name());
  if (s_omitnlosuffix) {
    size_t pos=fname.find("EW");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    pos=fname.find("QCD");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
  }
  if (!ReadInXSecs(m_resultpath+"/XS_"+fname)) return;
  ReadInHistogram(m_resultpath+"/WD_"+p_proc->Name());
  p_pshandler->ReadIn(m_resultpath+"/MC_"+fname);
  ReadInRB(m_resultpath+"/RB_"+p_proc->Name());
  SetTotal(0); 
}

bool Process_Integrator::WriteOutRB(const std::string &path)
{
  if (p_eeg) {
    std::ofstream rbmax((path+"/"+p_proc->Name()+".rb").c_str());
    if (!rbmax.good()) return false;
    rbmax.precision(12);
    rbmax<<m_rbsum<<" "<<m_rbsum2<<" "<<m_rbmax<<" "<<m_rbn<<std::endl;
    if (m_rbmax<1.0) {
      msg_Error()<<METHOD<<"(): ( (R/B)_{ME} / (R/B)_{PS} )_{max} < "
		 <<m_rbmax<<". Set to 1."<<std::endl;
      m_rbmax=1.0;
    }
  }
  bool res(true);
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      if (!(*p_proc)[i]->Integrator()->WriteOutRB(path)) res=false;
  return res;
}

bool Process_Integrator::ReadInRB(const std::string &path)
{
  if (p_eeg) {
    std::ifstream rbmax((path+"/"+p_proc->Name()+".rb").c_str());
    if (!rbmax.good()) return false;
    rbmax.precision(12);
    rbmax>>m_rbsum>>m_rbsum2>>m_rbmax>>m_rbn;
    if (m_rbmax<1.0) {
      msg_Error()<<METHOD<<"(): ( (R/B)_{ME} / (R/B)_{PS} )_{max} < "
		 <<m_rbmax<<". Set to 1."<<std::endl;
      m_rbmax=1.0;
    }
  }
  bool res(true);
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      if (!(*p_proc)[i]->Integrator()->ReadInRB(path)) res=false;
  return res;
}

void Process_Integrator::PrepareTerminate()  
{
  if (rpa.gen.BatchMode()) return;
  SetTotal();
  StoreResults();
}

