#include "PHASIC++/Main/Process_Integrator.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Smart_Pointer.C"
#include <stdio.h>
#include <string.h>
#include <errno.h>

using namespace PHASIC;
using namespace ATOOLS;

namespace ATOOLS { template class SP(Process_Integrator); }

static int s_whbins(1000);
static int s_genresdir(0);
 
Process_Integrator::Process_Integrator(Process_Base *const proc):
  p_proc(proc), p_pshandler(NULL),
  p_beamhandler(NULL), p_isrhandler(NULL),
  m_nin(0), m_nout(0), m_smode(0), m_swmode(0),
  m_threshold(0.), m_enhancefac(1.0), m_maxeps(0.0), m_rsfac(1.0),
  m_n(0), m_itmin(0), m_max(0.), m_totalxs(0.), 
  m_totalsum (0.), m_totalsumsqr(0.), m_totalerr(0.), m_ssum(0.), 
  m_ssumsqr(0.), m_smax(0.), m_ssigma2(0.), m_wmin(0.),
  m_mssum(0.), m_mssumsqr(0.), m_msn(0.), m_sn(0), m_son(1),
  m_writeout(false),
  p_whisto(NULL), p_colint(NULL), p_helint(NULL)
{
  m_colorscheme=cls::sum;
  m_helicityscheme=hls::sum;
}

bool Process_Integrator::Initialize
(BEAM::Beam_Spectra_Handler *const beamhandler,
 PDF::ISR_Handler *const isrhandler)
{
  m_nin=p_proc->NIn();
  m_nout=p_proc->NOut();
  p_momenta.resize(m_nin+m_nout);
  p_beamhandler=beamhandler;
  p_isrhandler=isrhandler;
  Data_Reader read(" ",";","!","=");
  m_swmode=read.GetValue<int>("SELECTION_WEIGHT_MODE", 0);
  static bool minit(false);
  if (!minit) {
    int smode;
    if (read.ReadFromFile(smode,"IB_SMODE")) {
      m_smode=smode;
      msg_Info()<<METHOD<<"(): Set sum mode = "<<m_smode<<".\n";
    }
    if (!read.ReadFromFile(s_whbins,"IB_WHBINS")) s_whbins=1000;
    else msg_Info()<<METHOD<<"(): Set weight histo bin number = "<<s_whbins<<".\n";
    if (read.ReadFromFile(s_genresdir,"GENERATE_RESULT_DIRECTORY"));
    else s_genresdir=1;
    minit=true;
  }
  return true;
}

Process_Integrator::~Process_Integrator()
{
  if (p_whisto!=NULL) delete p_whisto;
}

double Process_Integrator::SelectionWeight(const int mode) const
{
  if (!p_proc->IsGroup()) {
    if (mode!=0) return m_max*m_enhancefac;
    if (m_n+m_sn==0.0) return -1.0;
    if (m_totalxs==0.0) return 0.0;
    double selweight = m_swmode==0 ?
      sqrt((m_n+m_sn-1) * sqr(TotalVar()) + sqr(TotalResult())) :
      dabs(m_totalxs);
    return selweight*m_enhancefac;
  }
  double sw(0.0);
  for (size_t i(0);i<p_proc->Size();++i) {
    sw+=dabs((*p_proc)[i]->Integrator()->SelectionWeight(mode));
  }
  return sw;
}

double Process_Integrator::Sigma2() const
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
  return 1.0/
    ((p->m_ssumsqr/m_sn-sqr(p->m_ssum/m_sn))/(m_sn-1));
}

double Process_Integrator::TotalSigma2() const
{ 
  return m_ssigma2+Sigma2();
}

double Process_Integrator::TotalResult() const
{ 
  if (m_smode==0) return m_n+m_sn?(m_totalsum+m_ssum)/(m_n+m_sn):0.0;
  if (m_ssigma2==0.0) return m_sn?m_ssum/m_sn:0.0; 
  if (m_sn<2) return m_ssigma2?m_totalsum/m_ssigma2:0.0; 
  double s2(Sigma2());
  return m_ssigma2+s2?(m_totalsum+s2*m_ssum/m_sn)/(m_ssigma2+s2):0.0;
}

double Process_Integrator::TotalVar() const
{
  if (m_nin==1 && m_nout==2) return 0.;
  if (m_smode==0) {
    if (m_n+m_sn<2) return TotalResult();
    return sqrt(((m_totalsumsqr+m_ssumsqr)/(m_n+m_sn)
		 -sqr((m_totalsum+m_ssum)/(m_n+m_sn)))/(m_n+m_sn-1.0));
  }
  double s2(m_totalsumsqr);
  if (m_sn>1) {
    double vij2((m_sn-1)/
		(m_ssumsqr/m_sn-sqr(m_ssum/m_sn)));
    s2+=sqr(Sigma2())/vij2;
  }
  if (s2<0.0) return 0.0;
  return sqrt(s2)/TotalSigma2();
}

void Process_Integrator::OptimizeSubResult(const double &s2)
{
  m_n+=m_sn;
  if (m_smode==0) {
    m_totalsum+=m_ssum;
    m_totalsumsqr+=m_ssumsqr;
  }
  else {
    double vij2((m_sn-1)/
		(m_ssumsqr/m_sn-sqr(m_ssum/m_sn)));
    m_ssigma2+=s2; 
    m_totalsum+=s2*m_ssum/m_sn;
    m_totalsumsqr+=sqr(s2)/vij2;
  }
  m_ssum=m_ssumsqr=0.0;
  m_sn=0;
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->OptimizeSubResult(s2);
}

double Process_Integrator::RemainTimeFactor(double maxerr) 
{
  if (m_sn<1000) return 0.;
  return (sqr(1./(maxerr*TotalResult()))-m_ssigma2)/Sigma2();
}

void Process_Integrator::SetMomenta(const Vec4D_Vector &p)
{
  p_momenta=p;
  if (p_proc->Selected() && p_proc->Selected()!=p_proc)
    p_proc->Selected()->Integrator()->SetMomenta(p);
}

void Process_Integrator::SetMomenta(const Cluster_Amplitude &ampl)
{
  if (p_momenta.size()!=ampl.Legs().size()) {
    msg_Error()<<METHOD<<"("<<this<<"){\n  "
        <<"Cannot Set Momenta of Cluster_Amplitude "<<&ampl
        <<" because dimensions do not match.\n}\n";
    return;
  }
  for (size_t i(0);i<ampl.NIn();++i)
    p_momenta[i]=-ampl.Leg(i)->Mom();
  for (size_t i(ampl.NIn());i<p_momenta.size();++i)
    p_momenta[i]=ampl.Leg(i)->Mom();
  if (p_proc->Selected() && p_proc->Selected()!=p_proc)
    THROW(fatal_error,"Invalid function call");
}

void Process_Integrator::InitWeightHistogram()
{
  if (p_whisto) {
    delete p_whisto; 
    p_whisto=0;
  }

  double av(dabs(TotalResult()));
  if (!(av>0.)) {
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
  size_t vn;
  std::string name, dummy;
  My_In_File from(path+"/"+fname);
  if (!from.Open()) return false;
  from->precision(16);
  *from>>name>>m_totalxs>>m_max>>m_totalerr>>m_totalsum>>m_totalsumsqr
      >>m_n>>m_ssum>>m_ssumsqr>>m_smax>>m_ssigma2>>m_sn>>m_wmin
      >>m_son>>dummy>>dummy>>vn;
  if (name!=fname) THROW(fatal_error,"Corrupted results file");
  if (vn>100) {
    msg_Error()<<METHOD<<"(): Invalid vn in '"<<fname<<"'."<<std::endl;
  }
  else {
  m_vsmax.resize(vn);
  m_vsum.resize(vn);
  m_vsn.resize(vn);
  for (size_t i(0);i<m_vsn.size();++i)
    *from>>m_vsmax[i]>>m_vsum[i]>>m_vsn[i]>>dummy;
  }
  msg_Tracking()<<"Found result: xs for "<<name<<" : "
		<<m_totalxs*rpa->Picobarn()<<" pb"
		<<" +- ( "<<m_totalerr*rpa->Picobarn()<<" pb = "
		<<m_totalerr/m_totalxs*100.<<" % ) max: "
		<<m_max*rpa->Picobarn()<<std::endl;
  if (!p_proc->ReadIn(path)) return false;
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
  if (!FileExists(filename)) return;
  if (p_whisto) delete p_whisto; 
  p_whisto = new Histogram(filename);	
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->ReadInHistogram(dir);
}

void Process_Integrator::WriteOutXSecs(const std::string &path)
{
  std::string fname(p_proc->Name());
  My_Out_File outfile(path+"/"+fname);
  if (outfile.Open()) m_writeout=1;
  outfile->precision(16);
  *outfile<<fname<<"  "<<m_totalxs<<"  "<<m_max<<"  "
	 <<m_totalerr<<" "<<m_totalsum<<" "<<m_totalsumsqr<<" "
	 <<m_n<<" "<<m_ssum<<" "<<m_ssumsqr<<" "<<m_smax<<" "
	 <<m_ssigma2<<" "<<m_sn<<" "<<m_wmin<<" "<<m_son<<" "
	 <<-1<<" "<<-1<<"\n"<<m_vsn.size()<<"\n";
  for (size_t i(0);i<m_vsn.size();++i)
    *outfile<<m_vsmax[i]<<" "<<m_vsum[i]<<" "
	   <<m_vsn[i]<<" "<<-1<<"\n";
  p_proc->WriteOut(path);
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
    m_max=0.0;
    msg_Indent();
    for (size_t i(0);i<p_proc->Size();++i) {
      (*p_proc)[i]->Integrator()->SetTotal(msg_LevelIsTracking());
      m_max+=(*p_proc)[i]->Integrator()->Max();
    }
  }
  m_totalxs=TotalResult();
  m_totalerr=TotalVar();
  if (mode && m_totalxs!=0.0) {
    if (p_proc->NIn()==1) {
      msg_Info()<<om::bold<<p_proc->Name()<<om::reset<<" : "<<om::blue<<om::bold
                <<m_totalxs<<" GeV"<<om::reset<<" +- ( "<<om::red
                <<m_totalerr<<" GeV = "<<dabs(m_totalerr/m_totalxs)*100.
                <<" %"<<om::reset<<" ) "<<om::bold<<" exp. eff: "
                <<om::red<<(100.*dabs(m_totalxs/m_max))<<" %"<<om::reset<<std::endl;
    }
    else {
      msg_Info()<<om::bold<<p_proc->Name()<<om::reset<<" : "<<om::blue<<om::bold
                <<m_totalxs*rpa->Picobarn()<<" pb"<<om::reset<<" +- ( "<<om::red
                <<m_totalerr*rpa->Picobarn()<<" pb = "<<dabs(m_totalerr/m_totalxs)*100.
                <<" %"<<om::reset<<" ) "<<om::bold<<" exp. eff: "
                <<om::red<<(100.*dabs(m_totalxs/m_max))<<" %"<<om::reset<<std::endl;
    }
  }
}

double Process_Integrator::GetMaxEps(double epsilon)
{
  if (!p_whisto) return m_max;
  double res = dabs(TotalResult());
  double pxs = res*epsilon*p_whisto->Fills();
  double cutxs = 0.;
  double cnt = 0.;

  for (int i=p_whisto->Nbin()+1;i>0;i--) {
    cutxs+= p_whisto->Value(i) * 
      exp(log(10.)*(p_whisto->Xmin()+(i-0.5)*p_whisto->BinSize()));
    cnt+= p_whisto->Value(i);
    double twmax = exp(log(10.)*(p_whisto->Xmin()+(i-1)*p_whisto->BinSize()));
    if (cutxs-cnt*twmax>pxs) {
      return Min(exp(log(10.)*(p_whisto->Xmin()+i*p_whisto->BinSize())),dabs(m_max));
    }
  }
  return m_max;
}

void Process_Integrator::SetUpEnhance(const int omode) 
{
  if (m_maxeps>0.0 && !p_proc->IsGroup()) {
    double max(GetMaxEps(m_maxeps));
    if (omode)
      msg_Info()<<"  reduce max for "<<p_proc->Name()<<" to "
		<<max/Max()<<" ( eps = "<<m_maxeps<<" ) "<<std::endl;
    SetMax(max);
  }
  if (p_proc->IsGroup()) {
    double oldmax(m_max);
    m_max=0.0;
    for (size_t i(0);i<p_proc->Size();++i) {
      (*p_proc)[i]->Integrator()->SetUpEnhance(msg_LevelIsTracking());
      m_max+=(*p_proc)[i]->Integrator()->Max();
    }
    if (omode || p_proc->Parent()==p_proc)
      if (p_whisto)
    msg_Info()<<"  reduce max for "<<p_proc->Name()<<" to "
	      <<m_max/oldmax<<" ( eps = "<<m_maxeps<<" ) "<<std::endl;
  }

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

void Process_Integrator::SetRSEnhanceFactor(const double &rsfac)
{ 
  m_rsfac=rsfac; 
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->SetRSEnhanceFactor(rsfac);
}

void Process_Integrator::AddPoint(const double value) 
{
#ifdef USING__MPI
  m_msn++;
  m_mssum    += value;
  m_mssumsqr += sqr(value);
#else
  m_sn++;
  m_ssum    += value;
  m_ssumsqr += sqr(value);
#endif
  double cur=value*p_pshandler->Enhance();
  double max=dabs(cur)/dabs(p_proc->Last())*
    ATOOLS::Max(p_proc->LastPlus(),-p_proc->LastMinus());
  if (max>m_max)  m_max  = max;
  if (max>m_smax) m_smax = max;
  if (p_whisto) {
    if(cur!=0.) p_whisto->Insert(max,1.0/p_pshandler->Enhance()); /*TODO*/
    else p_whisto->Insert(1.0,0.0);
  }
  p_proc->AddPoint(value);
  if (p_proc->IsGroup()) {
    if (p_proc->Last()==0.0 || value==0.0)
      for (size_t i(0);i<p_proc->Size();++i)
	(*p_proc)[i]->Integrator()->AddPoint(0.0);
    else
      for (size_t i(0);i<p_proc->Size();++i)
	(*p_proc)[i]->Integrator()->
	  AddPoint(value*(*p_proc)[i]->Last()/p_proc->Last());
  }
}

void Process_Integrator::SetMax(const double max) 
{
  m_max=max;
  if (!p_proc->IsGroup()) return;
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

void Process_Integrator::Reset()
{
  m_n=0;
  m_totalxs=m_totalsum=m_totalsumsqr=m_totalerr=0.0;
  m_smax=m_max=m_wmin=m_ssigma2=m_ssumsqr=m_ssum=0.0;
  m_sn=0;
  m_son=1;
  m_vsmax.clear(); 
  m_vsn.clear();   
  m_vsum.clear(); 
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i) 
      (*p_proc)[i]->Integrator()->Reset();
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
    m_max=0.0;
    return;
  }
  if (flag==0) {
    if (m_vsmax.size()>1) {
      m_vsmax.erase(m_vsmax.begin());
      m_vsn.erase(m_vsn.begin());
      m_vsum.erase(m_vsum.begin());
    }
    if (m_vsmax.empty()) {
      m_vsmax.push_back(m_max);
      m_vsn.push_back(m_n);
      m_vsum.push_back(m_ssum);
    }
    m_vsmax.back() = ATOOLS::Max(m_smax,m_vsmax.back());
    m_vsn.back()   = m_n;
    m_vsum.back()  = m_ssum;
  }
  else {
    if (flag==2 && m_vsmax.size()==4) {
      m_vsmax.erase(m_vsmax.begin());
      m_vsn.erase(m_vsn.begin());
      m_vsum.erase(m_vsum.begin());
    }
    m_vsmax.push_back(m_smax);
    m_vsn.push_back(m_n);
    m_vsum.push_back(m_ssum);
    if (flag==2) m_smax = 0.;
  }
  m_max=0.0;
  for (size_t i=0;i<m_vsmax.size();i++) {
    m_max=ATOOLS::Max(m_max,m_vsmax[i]);
  }
} 

void Process_Integrator::SetPSHandler(const SP(Phase_Space_Handler) &pshandler)
{
  p_pshandler=pshandler;
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->SetPSHandler(pshandler);
} 

void Process_Integrator::MPICollect
(std::vector<double> &sv,std::vector<double> &mv,size_t &i)
{
  sv.resize(3*(i+1));
  mv.resize(2*(i+1));
  sv[3*i+0]=m_msn;
  sv[3*i+1]=m_mssum;
  sv[3*i+2]=m_mssumsqr;
  mv[2*i+0]=m_max;
  mv[2*i+1]=m_smax;
  ++i;
  if (p_proc->IsGroup())
    for (size_t j(0);j<p_proc->Size();++j)
      (*p_proc)[j]->Integrator()->MPICollect(sv,mv,i);
}

void Process_Integrator::MPIReturn
(std::vector<double> &sv,std::vector<double> &mv,size_t &i)
{
  m_msn=sv[3*i+0];
  m_mssum=sv[3*i+1];
  m_mssumsqr=sv[3*i+2];
  m_max=mv[2*i+0];
  m_smax=mv[2*i+1];
  ++i;
  if (p_proc->IsGroup())
    for (size_t j(0);j<p_proc->Size();++j)
      (*p_proc)[j]->Integrator()->MPIReturn(sv,mv,i);
}

void Process_Integrator::MPISync(const int mode)
{
  if (p_whisto) p_whisto->MPISync();
#ifdef USING__MPI
  if (mode==0) {
    size_t i(0), j(0);
    std::vector<double> sv, mv;
    MPICollect(sv,mv,i);
    if (MPI::COMM_WORLD.Get_size()) {
      mpi->MPIComm()->Allreduce
	(MPI_IN_PLACE,&sv[0],sv.size(),MPI::DOUBLE,MPI::SUM);
      mpi->MPIComm()->Allreduce
	(MPI_IN_PLACE,&mv[0],mv.size(),MPI::DOUBLE,MPI::MAX);
    }
    MPIReturn(sv,mv,j);
  }
  m_sn+=m_msn;
  m_ssum+=m_mssum;
  m_ssumsqr+=m_mssumsqr;
  m_msn=m_mssum=m_mssumsqr=0.0;
#endif
  p_proc->MPISync(mode);
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->MPISync(1);
}

void Process_Integrator::OptimizeResult()
{
  OptimizeSubResult(Sigma2());
} 

void Process_Integrator::EndOptimize()
{
  p_proc->EndOptimize();
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->EndOptimize();
} 

void Process_Integrator::SetISRThreshold(const double threshold) 
{
  m_threshold=threshold;
}

void Process_Integrator::StoreBackupResults()
{
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()) return;
#endif
  if (!FileExists(m_resultpath+".db")) return;
  if (!Copy(m_resultpath+".db",m_resultpath+".db.bak",true))
    msg_Error()<<METHOD<<"(): Copy error. "
	       <<strerror(errno)<<"."<<std::endl;
}

void Process_Integrator::StoreResults(const int mode)
{
  if (m_msn) MPISync();
  if (m_resultpath.length()==0) return;
  if (m_totalxs!=0.0 && mode==0) return;
  SetTotal(0); 
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()) return;
#endif
  My_In_File::ExecDB(m_resultpath+"/","begin");
  std::string fname(p_proc->Name());
  WriteOutXSecs(m_resultpath+"/"+p_proc->Generator()->Name()+"/XS_"+fname);
  WriteOutHistogram(m_resultpath+"/"+p_proc->Generator()->Name()+"/WD_"+fname);
  p_pshandler->WriteOut(m_resultpath+"/"+p_proc->Generator()->Name()+"/MC_"+fname);
  My_In_File::ExecDB(m_resultpath+"/","commit");
  StoreBackupResults();
}

void Process_Integrator::ReadResults()
{
  if (m_resultpath.length()==0) return;
  std::string fname(p_proc->Name());
  if (!ReadInXSecs(m_resultpath+"/"+p_proc->Generator()->Name()+"/XS_"+fname)) return;
  ReadInHistogram(m_resultpath+"/"+p_proc->Generator()->Name()+"/WD_"+fname);
  p_pshandler->ReadIn(m_resultpath+"/"+p_proc->Generator()->Name()+"/MC_"+fname);
  SetTotal(0); 
}

void Process_Integrator::PrepareTerminate()  
{
  if (rpa->gen.BatchMode()&1) return;
  SetTotal();
  StoreResults();
}

