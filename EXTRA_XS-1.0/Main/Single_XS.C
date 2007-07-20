#include "Single_XS.H"
#include "ISR_Handler.H"
#include "Phase_Space_Handler.H"
#include "Run_Parameter.H"
#include "Regulator_Base.H"
#include "Message.H"

using namespace EXTRAXS;
using namespace ATOOLS;

Single_XS::Single_XS(const size_t nin,const size_t nout,
		     const ATOOLS::Flavour *flavours,
		     const PHASIC::scl::scheme scalescheme,
		     const int kfactorscheme,
		     BEAM::Beam_Spectra_Handler *const beamhandler,
		     PDF::ISR_Handler *const isrhandler,
		     ATOOLS::Selector_Data *const selectordata,
		     XS_Model_Base *const model):
  XS_Base(nin,nout,flavours,scalescheme,kfactorscheme,
	  beamhandler,isrhandler,selectordata,model)
{
  p_selected=this;
}

Single_XS::Single_XS(const size_t nin,const size_t nout,
		     const ATOOLS::Flavour *flavours,
		     XS_Model_Base *const model):
  XS_Base(nin,nout,flavours,model)
{
  p_selected=this;
}

void Single_XS::WriteOutXSecs(std::ofstream &outfile)
{
  outfile.precision(12);
  outfile<<m_name<<"  "<<m_totalxs<<"  "<<m_max<<"  "<<m_totalerr<<" "
	 <<m_totalsum<<" "<<m_totalsumsqr<<" "<<m_n<<" "
	 <<m_ssum<<" "<<m_ssumsqr<<" "<<m_ssigma2<<" "<<m_sn<<" "<<m_wmin<<" "<<m_son<<std::endl; 
}

bool Single_XS::CalculateTotalXSec(const std::string &resultpath,
				   const bool create) 
{ 
  Reset();
  if (p_isrhandler) {
    if (m_nin==2) {
      if (p_flavours[0].Mass()!=p_isrhandler->Flav(0).Mass() ||
	  p_flavours[1].Mass()!=p_isrhandler->Flav(1).Mass()) {
	p_isrhandler->SetPartonMasses(p_flavours);
      }
    }
  }
  p_pshandler->InitCuts();
  if (p_isrhandler)
    p_isrhandler->SetSprimeMin(p_pshandler->Cuts()->Smin());
  CreateFSRChannels();
  if (!m_channels) {
    p_pshandler->CreateIntegrators();
    CreateISRChannels();
    m_channels = true;
  }
  std::string filename(resultpath+"/"+m_name+".xstotal"), name;
  if (resultpath!="") {
    double totalxs,totalerr,max,sum,sqrsum,ssum,ssqrsum,ss2,wmin;
    long int n,sn,son;
    std::ifstream from(filename.c_str());
    from>>name>>totalxs>>max>>totalerr>>sum>>sqrsum>>n
	>>ssum>>ssqrsum>>ss2>>sn>>wmin>>son;
    if (name==m_name) {
      m_totalxs  = totalxs;
      m_totalerr = totalerr; 
      m_max      = max;
      m_n        = n;
      m_totalsum = sum;
      m_totalsumsqr = sqrsum;
      m_vsmax.clear(); 
      m_vsmax.push_back(max);
      m_vsn.clear();   
      m_vsn.push_back(n);
      m_sn = sn; m_smax=0.;
      m_ssum     = ssum;
      m_ssumsqr  = ssqrsum;
      m_ssigma2  = ss2;
      m_wmin     = wmin;
      m_son      = son;
      msg_Info()<<"Found result : xs for "<<m_name<<" : "
		<<m_totalxs*rpa.Picobarn()<<" pb"
		<<" +/- "<<m_totalerr/m_totalxs*100.
		<<"%, max : "<<m_max<<std::endl;
      p_pshandler->ReadIn(resultpath+"/MC_"+m_name);
    }
    from.close();
    if (p_pshandler->BeamIntegrator()) 
      p_pshandler->BeamIntegrator()->Print();
    if (p_pshandler->ISRIntegrator()) 
      p_pshandler->ISRIntegrator()->Print();
    if (p_pshandler->FSRIntegrator()) 
      p_pshandler->FSRIntegrator()->Print();
  }
  m_resultpath=resultpath;
  m_resultfile=filename;
  exh->AddTerminatorObject(this);
  p_pshandler->InitIncoming();
  double var=TotalVar();
  m_totalxs=p_pshandler->Integrate()/rpa.Picobarn(); 
  if (!(IsZero((m_totalxs-TotalResult())/(m_totalxs+TotalResult())))) {
    msg_Error()<<"Result of PS-Integrator and internal summation do not coincide!"<<std::endl
	       <<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<TotalResult()<<std::endl;
  }
  if (m_totalxs>0.) {
    SetTotal();
    if (var==TotalVar()) {
      ATOOLS::exh->RemoveTerminatorObject(this);
      return 1;
    }
    if (resultpath!=std::string("")) {
      std::ofstream to;
      to.open(filename.c_str(),std::ios::out);
      to.precision(12);
      msg_Info()<<"Store result : xs for "<<m_name<<" : ";
      WriteOutXSecs(to);
      if (m_nin==2) msg_Info()<<m_totalxs*ATOOLS::rpa.Picobarn()<<" pb";
      if (m_nin==1) msg_Info()<<m_totalxs<<" GeV";
      msg_Info()<<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<std::endl
		<<"  max : "<<m_max<<std::endl;
      p_pshandler->WriteOut(resultpath+std::string("/MC_")
			    +m_name,create);
      to.close();
    }
    ATOOLS::exh->RemoveTerminatorObject(this);
    return 1;
  }
  ATOOLS::exh->RemoveTerminatorObject(this);
  return 0;
}

ATOOLS::Blob_Data_Base *Single_XS::OneEvent() 
{ 
  return p_activepshandler->OneEvent(); 
}

ATOOLS::Blob_Data_Base *Single_XS::WeightedEvent(const int mode) 
{ 
  return p_activepshandler->WeightedEvent(mode); 
}

void Single_XS::SetTotal(const int mode) 
{
  m_totalxs=TotalResult();//m_totalsum/m_n; 
  m_totalerr=TotalVar();//sqrt((m_n*m_totalsumsqr-ATOOLS::sqr(m_totalsum))/(m_n-1))/m_n;
  if (mode&1)
    msg_Info()<<om::bold<<m_name<<om::reset<<" : "<<om::blue<<om::bold
	      <<m_totalxs*rpa.Picobarn()<<" pb"<<om::reset
	      <<" +/- "<<om::reset<<om::blue<<m_totalerr/m_totalxs*100.
	      <<" %,"<<om::reset<<om::bold<<" exp. eff: "
	      <<om::red<<(100.*m_totalxs/m_max)<<" %"<<om::reset<<std::endl;
}

void Single_XS::OptimizeResult()
{
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
//   cout<<"Weights (actual/min) "<<ssigma2<<" "<<m_wmin<<endl;
  m_son++;
}

void Single_XS::ResetMax(int flag)
{
  if (flag==3) {
    m_vsmax.clear();
    m_vsn.clear();
    m_max=0.0;
    return;
  }
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
//   cout<<"Single_Process::ResetMax "<<flag<<": "<<endl;
//   for (size_t i=0;i<m_vsmax.size();i++) cout<<m_vsmax[i]<<" "<<m_vsn[i]<<endl;
  for (size_t i=0;i<m_vsmax.size();i++) m_max=ATOOLS::Max(m_max,m_vsmax[i]);
}

double Single_XS::Differential(const double s,const double t,const double u)
{
  m_lastdxs=(*p_regulator)((*this)(s,t,u));
  if (m_lastdxs<=0.) return m_lastdxs=m_last=0.;
  if (p_isrhandler && m_nin==2) { 
    if (p_isrhandler->On()) m_lastlumi=p_isrhandler->Weight(p_flavours); 
    else m_lastlumi=1.;
  }
  else m_lastlumi=1.;
  return m_last=m_lastdxs*m_lastlumi*KFactor();
}

double Single_XS::Differential2() 
{
  if (p_isrhandler && m_nin==2) {
    if (p_flavours[0]==p_flavours[1] || p_isrhandler->On()==0) return 0.;
    double tmp=m_lastdxs*p_isrhandler->Weight2(p_flavours); 
    m_last+=tmp;
    return tmp*KFactor();
  }
  return 0;
}

double Single_XS::operator()(const double s,const double t,const double u) 
{
  msg_Error()<<"Single_XS::operator()("<<s<<","<<t<<","<<u<<"): "
		     <<"Virtual method called!"<<std::endl;
  return 0.;
}

void Single_XS::SetISR(PDF::ISR_Handler *const isrhandler) 
{ 
  p_isrhandler=isrhandler; 
}

XS_Base *const Single_XS::operator[](const size_t i) const 
{
  return (XS_Base*)this;
}

size_t Single_XS::Size() const
{ 
  return 1; 
}

bool Single_XS::Tests()
{
  p_activepshandler->TestPoint(p_momenta);
  return true;
}
