#include "PHASIC++/Main/Process_Integrator.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <errno.h>

using namespace PHASIC;
using namespace ATOOLS;

static int s_whbins(600);
static int s_genresdir(0);

Process_Integrator::Process_Integrator(Process_Base *const proc):
  p_proc(proc),
  p_beamhandler(NULL), p_isrhandler(NULL),
  m_nin(0), m_nout(0), m_swmode(0),
  m_threshold(0.), m_enhancefac(1.0), m_maxeps(0.0), m_rsfac(1.0),
  m_itmin(0), m_start(0), m_totalxs(0.), m_totalerr(0.), m_max(0.),
  m_sn(0.), m_ssum(0.), m_ssumsqr(0.), m_smax(0.),
  m_msn(0.), m_msv(0.), m_mssum(0.), m_mssumsqr(0.),
  m_writeout(false),
  p_whisto(NULL)
{
  m_colorscheme=cls::sum;
  m_helicityscheme=hls::sum;
}

bool Process_Integrator::Initialize
(BEAM::Beam_Spectra_Handler *const beamhandler,
 PDF::ISR_Handler *const isrhandler)
{
  Settings& s = Settings::GetMainSettings();
  m_nin=p_proc->NIn();
  m_nout=p_proc->NOut();
  p_momenta.resize(m_nin+m_nout);
  p_beamhandler=beamhandler;
  p_isrhandler=isrhandler;
  m_swmode=s["SELECTION_WEIGHT_MODE"].SetDefault(0).Get<int>();
  if (s["HDF5_UNWEIGHT"].SetDefault(0).Get<int>()!=0) m_swmode=2;
  if (m_swmode!=0) msg_Info()<<METHOD<<"(): Set selection weight mode "<<m_swmode<<".\n";
  static bool minit(false);
  if (!minit) {
    // weight histo bin number
    s_whbins = s["IB_WHBINS"].SetDefault(100).Get<int>();
    minit=true;
  }
  return true;
}

Process_Integrator::~Process_Integrator()
{
  if (p_whisto!=NULL) delete p_whisto;
}

double Process_Integrator::N() const
{
  double n(0);
  for (size_t i(m_start);i<m_vn.size();++i) n+=m_vn[i];
  return n;
}

double Process_Integrator::Sum() const
{
  double sum(0);
  for (size_t i(m_start);i<m_vsum.size();++i) sum+=m_vsum[i];
  return sum;
}

double Process_Integrator::SumSqr() const
{
  double sumsqr(0);
  for (size_t i(m_start);i<m_vsumsqr.size();++i) sumsqr+=m_vsumsqr[i];
  return sumsqr;
}

double Process_Integrator::GetMax() const
{
  double max(m_smax);
  for (size_t i(m_start);i<m_vmax.size();++i)
    max=ATOOLS::Max(max,m_vmax[i]);
  return max;
}

double Process_Integrator::SelectionWeight(const int mode) const
{
  if (!p_proc->IsGroup()) {
    if (mode!=0 || m_swmode==2) return m_max*m_enhancefac;
    if (N()==0) return -1.0;
    if (m_totalxs==0.0) return 0.0;
    double selweight=m_swmode?dabs(m_totalxs):sqrt(SumSqr()/N());
    return selweight*m_enhancefac;
  }
  double sw(0.0);
  for (size_t i(0);i<p_proc->Size();++i) {
    sw+=dabs((*p_proc)[i]->Integrator()->SelectionWeight(mode));
  }
  return sw;
}

double Process_Integrator::Result() const
{ 
  double n(N());
  if (n==0) return 0.0;
  return Sum()/n;
}

double Process_Integrator::Variance() const
{
  if (m_nin==1 && m_nout==2) return 0.;
  double n(N());
  if (n<2) return Result();
  return sqrt((SumSqr()/n-sqr(Sum()/n))/(n-1.0));
}

double Process_Integrator::RemainTimeFactor(double maxerr) 
{
  if (m_sn<1000) return 0.;
  return sqr(Variance()/(maxerr*Result()));
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
  double av(dabs(Result()));
  if (IsBad(av)) {
    msg_Error()<<METHOD<<"(): Average = "<<av
	       <<" in "<<p_proc->ResultsName()<<std::endl;
    return;
  }
  /* If av=0, then subprocess at hand does not contribute.
     In this case, set av to arbitrary value to avoid nans in
     following histogram */
  if (IsZero(av)) av=1.;
  av = exp(log(10.)*int(log(av)/log(10.)+0.5));
  p_whisto = new Histogram(10,av*1.e-9,av*1.e6,s_whbins);
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->InitWeightHistogram();
}

bool Process_Integrator::ReadInXSecs(const std::string &path)
{
  std::string fname(p_proc->ResultsName());
  size_t vn;
  std::string name, dummy;
  My_In_File from(path+"/"+fname);
  if (!from.Open()) return false;
  from->precision(16);
  *from>>name>>m_totalxs>>m_max>>m_totalerr
       >>m_sn>>m_sv>>m_ssum>>m_ssumsqr>>m_smax
       >>vn>>m_start;
  if (name!=fname) THROW(fatal_error,"Corrupted results file");
  m_vn.resize(vn);
  m_vsum.resize(vn);
  m_vsumsqr.resize(vn);
  m_vmax.resize(vn);
  m_vw.resize(vn);
  for (size_t i(0);i<m_vn.size();++i)
    *from>>m_vn[i]>>m_vsum[i]>>m_vsumsqr[i]>>m_vmax[i]>>m_vw[i];
  msg_Tracking()<<"Found result: xs for "<<name<<" : "
		<<m_totalxs*rpa->Picobarn()<<" pb"
		<<" +- ( "<<m_totalerr*rpa->Picobarn()<<" pb = "
		<<m_totalerr/m_totalxs*100.<<" % ) max: "
		<<m_max*rpa->Picobarn()<<std::endl;
  if (!p_proc->ReadIn(path)) return false;
  if (p_colint!=NULL) p_colint->ReadIn(path+"/"+fname+"_Color");
  bool res(true);
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      if (!(*p_proc)[i]->Integrator()->ReadInXSecs(path)) res=false;
  SetMax(m_max);
  return res;
}

void Process_Integrator::ReadInHistogram(std::string dir)
{
  std::string filename = dir+"/"+p_proc->ResultsName();
  if (!FileExists(filename)) return;
  if (p_whisto) delete p_whisto; 
  p_whisto = new Histogram(filename);	
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->ReadInHistogram(dir);
}

void Process_Integrator::WriteOutXSecs(const std::string &path)
{
  std::string fname(p_proc->ResultsName());
  My_Out_File outfile(path+"/"+fname);
  if (outfile.Open()) m_writeout=1;
  outfile->precision(16);
  *outfile<<fname<<"  "<<m_totalxs<<"  "<<m_max<<"  "
	  <<m_totalerr<<" "<<m_sn<<" "<<m_sv<<" "
	  <<m_ssum<<" "<<m_ssumsqr<<" "<<m_smax<<"\n"
	  <<m_vn.size()<<" "<<m_start<<"\n";
  for (size_t i(0);i<m_vn.size();++i)
    *outfile<<m_vn[i]<<" "<<m_vsum[i]<<" "<<m_vsumsqr[i]
	    <<" "<<m_vmax[i]<<" "<<m_vw[i]<<"\n";
  p_proc->WriteOut(path);
  if (p_colint!=NULL) p_colint->WriteOut(path+"/"+fname+"_Color");
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->WriteOutXSecs(path);
}

void Process_Integrator::WriteOutHistogram(std::string dir)
{
  if (p_whisto) p_whisto->Output(dir+"/"+p_proc->ResultsName());
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->WriteOutHistogram(dir);
}

void Process_Integrator::SetTotal(const int mode)  
{ 
  if (m_msn) MPISync();
  if (p_proc->IsGroup()) {
    m_max=0.0;
    msg_Indent();
    for (size_t i(0);i<p_proc->Size();++i) {
      (*p_proc)[i]->Integrator()->SetTotal(msg_LevelIsTracking());
      m_max+=(*p_proc)[i]->Integrator()->Max();
    }
  }
  m_totalxs=Result();
  m_totalerr=Variance();
  SetMax(GetMax());
  if (mode && m_totalxs!=0.0) {
    if (p_proc->NIn()==1) {
      msg_Info()<<om::bold<<p_proc->ResultsName()<<om::reset<<" : "<<om::blue<<om::bold
                <<m_totalxs<<" GeV"<<om::reset<<" +- ( "<<om::red
                <<m_totalerr<<" GeV = "<<dabs(m_totalerr/m_totalxs)*100.
                <<" %"<<om::reset<<" ) "<<om::bold<<" exp. eff: "
                <<om::red<<dabs(m_totalxs/m_max)<<om::reset<<std::endl;
    }
    else {
      msg_Info()<<om::bold<<p_proc->ResultsName()<<om::reset<<" : "<<om::blue<<om::bold
                <<m_totalxs*rpa->Picobarn()<<" pb"<<om::reset<<" +- ( "<<om::red
                <<m_totalerr*rpa->Picobarn()<<" pb = "<<dabs(m_totalerr/m_totalxs)*100.
                <<" %"<<om::reset<<" ) "<<om::bold<<" exp. eff: "
                <<om::red<<dabs(m_totalxs/m_max)<<om::reset<<std::endl;
    }
  }
}

double Process_Integrator::GetMaxEps(double epsilon)
{
  if (!p_whisto) return m_max;
  if (epsilon<=-1.) {
    int nsamples(-epsilon), npoints(p_whisto->Fills());
    npoints*=-(epsilon-int(epsilon));
#ifdef USING__MPI
    nsamples=std::max(1,nsamples/mpi->Size());
#endif
    double nonzero(0.);
    for (size_t i(0);i<p_whisto->Nbin();++i) nonzero+=p_whisto->Value(i);
    nonzero*=npoints/p_whisto->Fills();
    std::vector<double> maxs(nsamples,0.0);
    for (size_t j(0);j<nsamples;++j)
      for (size_t i(0);i<nonzero;++i) {
	double x=p_whisto->GeneratePoint(ran->Get());
	if (x>maxs[j]) maxs[j]=x;
      }
    std::sort(maxs.begin(),maxs.end(),std::less<double>());
#ifdef USING__MPI
    mpi->Allreduce(&maxs[maxs.size()/2],1,MPI_DOUBLE,MPI_MAX);
#endif
    return maxs[maxs.size()/2];
  }

  double res = dabs(Result());
  double pxs = res*dabs(epsilon)*p_whisto->Fills();
  double cutxs = 0.;
  double cnt = 0.;

  for (int i=p_whisto->Nbin()+1;i>0;i--) {
    cutxs+= p_whisto->Value(i) * 
      exp(log(10.)*(p_whisto->Xmin()+(i-0.5)*p_whisto->BinSize()));
    cnt+= p_whisto->Value(i);
    double twmax = exp(log(10.)*(p_whisto->Xmin()+(i-1)*p_whisto->BinSize()));
    if (epsilon<0.) {
      if (cnt*twmax>pxs) {
	return Min(exp(log(10.)*(p_whisto->Xmin()+i*p_whisto->BinSize())),dabs(m_max));
      }
    }
    else {
      if (cutxs-cnt*twmax>pxs) {
	return Min(exp(log(10.)*(p_whisto->Xmin()+i*p_whisto->BinSize())),dabs(m_max));
      }
    }
  }
  return m_max;
}

void Process_Integrator::SetUpEnhance(const int omode) 
{
  if (m_maxeps!=0.0 && !p_proc->IsGroup()) {
    double max(GetMaxEps(m_maxeps));
    if (omode || msg->LevelIsTracking())
      msg_Info()<<"  reduce max for "<<p_proc->ResultsName()<<" to "
		<<max/Max()<<" ( eps = "<<m_maxeps<<" -> exp. eff "
                <<dabs(m_totalxs/max)<<" ) "<<std::endl;
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
    msg_Info()<<"  reduce max for "<<p_proc->ResultsName()<<" to "
	      <<m_max/oldmax<<" ( eps = "<<m_maxeps<<" -> exp. eff "
	      <<dabs(m_totalxs/m_max)<<" ) "<<std::endl;
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
  double enhance = p_pshandler->EnhanceWeight();
#ifdef USING__MPI
  m_msn++;
  if (value) m_msv++;
  m_mssum    += value/enhance;
  m_mssumsqr += sqr(value/enhance);
#else
  m_sn++;
  if (value) m_sv++;
  m_ssum    += value/enhance;
  m_ssumsqr += sqr(value/enhance);
#endif
  double enh=p_pshandler!=NULL?p_pshandler->Enhance():1;
  double cur=value*enh;
  double max=dabs(cur)/dabs(p_proc->Last())*
    ATOOLS::Max(p_proc->LastPlus(),-p_proc->LastMinus());
  if (max>m_max)  m_max  = max;
  if (max>m_smax) m_smax = max;
  if (p_whisto) {
    if(cur!=0.) p_whisto->Insert(max,1.0/enh); /*TODO*/
    else p_whisto->Insert(1.0,0.0);
  }
  if (p_colint!=NULL) p_colint->AddPoint(value);
  p_proc->AddPoint(value);
  if (p_helint!=NULL) p_helint->AddPoint(value);
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
		 <<p_proc->ResultsName()<<".\n  sum = "<<sum
		 <<" vs. total = "<<m_totalxs<<" ("
		 <<((sum-m_totalxs)/m_totalxs)<<")"<<std::endl;
      msg_Error().precision(6);
    }
    m_totalxs=sum;
  }
} 

void Process_Integrator::Reset(const int mode)
{
  m_totalxs=m_totalerr=m_max=0.0;
  m_sn=m_sv=m_ssum=m_ssumsqr=m_smax=0.0;
  m_msn=m_msv=m_mssum=m_mssumsqr=0.0;
  m_vn.clear();
  m_vsum.clear(); 
  m_vsumsqr.clear();
  m_vmax.clear();
  m_vw.clear();
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i) 
      (*p_proc)[i]->Integrator()->Reset(mode);
}

void Process_Integrator::SetPSHandler(Phase_Space_Handler * pshandler) {
  //  const std::shared_ptr<Phase_Space_Handler> &pshandler) {
  p_pshandler=pshandler;
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->SetPSHandler(pshandler);
} 

void Process_Integrator::MPICollect
(std::vector<double> &sv,std::vector<double> &mv,size_t &i)
{
  sv.resize(4*(i+1));
  mv.resize(2*(i+1));
  sv[4*i+0]=m_msn;
  sv[4*i+1]=m_msv;
  sv[4*i+2]=m_mssum;
  sv[4*i+3]=m_mssumsqr;
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
  m_msn=sv[4*i+0];
  m_msv=sv[4*i+1];
  m_mssum=sv[4*i+2];
  m_mssumsqr=sv[4*i+3];
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
    if (mpi->Size()) {
      mpi->Allreduce(&sv[0],sv.size(),MPI_DOUBLE,MPI_SUM);
      mpi->Allreduce(&mv[0],mv.size(),MPI_DOUBLE,MPI_MAX);
    }
    MPIReturn(sv,mv,j);
  }
  m_sn+=m_msn;
  m_sv+=m_msv;
  m_ssum+=m_mssum;
  m_ssumsqr+=m_mssumsqr;
  m_msn=m_msv=m_mssum=m_mssumsqr=0.0;
#endif
  if (p_colint!=NULL) p_colint->MPISync();
  p_proc->MPISync(mode);
  if (p_helint!=NULL) p_helint->MPISync();
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->MPISync(1);
}

void Process_Integrator::OptimizeResult(const double &s2)
{
  m_vn.push_back(m_sn);
  m_vsum.push_back(m_ssum);
  m_vsumsqr.push_back(m_ssumsqr);
  m_vmax.push_back(m_smax);
  m_vw.push_back(s2);
  m_sn=m_sv=m_ssum=m_ssumsqr=m_smax=0.0;
  if (p_helint!=NULL) p_helint->Optimize();
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->OptimizeResult(s2);
}

void Process_Integrator::OptimizeResult()
{
  if (m_msn) MPISync();
  double s2((m_ssumsqr/m_sn-sqr(m_ssum/m_sn))/(m_sn-1));
  if (m_sn) OptimizeResult(1/s2);
} 

void Process_Integrator::EndOptimize()
{
  m_start=m_vn.size();
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
  if (!FileExists(m_resultpath+".zip")) return;
  if (!Copy(m_resultpath+".zip",m_resultpath+".zip~",true))
    msg_Error()<<METHOD<<"(): Copy error. "
	       <<strerror(errno)<<"."<<std::endl;
}

void Process_Integrator::StoreResults(const int mode)
{
  if (m_resultpath.length()==0) return;
  if (m_totalxs!=0.0 && mode==0) return;
  std::string fname(p_proc->ResultsName());
  WriteOutXSecs(m_resultpath+"/"+p_proc->Generator()->Name()+"/XS_"+fname);
  WriteOutHistogram(m_resultpath+"/"+p_proc->Generator()->Name()+"/WD_"+fname);
  if (p_pshandler!=NULL) p_pshandler->WriteOut(m_resultpath+"/"+p_proc->Generator()->Name()+"/MC_"+fname);
  if (p_helint!=NULL) p_helint->WriteOut(m_resultpath+"/"+p_proc->Generator()->Name()+"/HI_"+fname);
  My_In_File::CloseDB(m_resultpath+"/",0);
  StoreBackupResults();
}

void Process_Integrator::ReadResults()
{
  if (m_resultpath.length()==0) return;
  std::string fname(p_proc->ResultsName());
  if (!ReadInXSecs(m_resultpath+"/"+p_proc->Generator()->Name()+"/XS_"+fname)) return;
  ReadInHistogram(m_resultpath+"/"+p_proc->Generator()->Name()+"/WD_"+fname);
  p_pshandler->ReadIn(m_resultpath+"/"+p_proc->Generator()->Name()+"/MC_"+fname);
  if (p_helint!=NULL) p_helint->ReadIn(m_resultpath+"/"+p_proc->Generator()->Name()+"/HI_"+fname);
  SetTotal(0);
}

void Process_Integrator::PrepareTerminate()  
{
  if (rpa->gen.BatchMode()&1) return;
  SetTotal();
  StoreResults();
}

