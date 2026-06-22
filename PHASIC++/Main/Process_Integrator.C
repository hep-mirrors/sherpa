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

static int s_whbins(1300);
static int s_genresdir(0);

Process_Integrator::Process_Integrator(Process_Base *const proc):
  p_proc(proc), p_pshandler(NULL),
  p_beamhandler(NULL), p_isrhandler(NULL),
  m_nin(0), m_nout(0), m_smode(0), m_swmode(0),
  m_threshold(0.), m_enhancefac(1.0), m_maxeps(0.0), m_rsfac(1.0),
  m_n(0), m_itmin(0), m_itmax(1000000), m_max(0.), m_weightmax(0.), m_totalxs(0.),
  m_totalsum (0.), m_totalsumsqr(0.), m_totalerr(0.), m_ssum(0.), m_ssumenh(0.),
  m_ssumenhabs(0.), m_ssumsqr(0.), m_ssumenhsqr(0.), m_smax(0.), m_ssigma2(0.), m_wmin(0.),
  m_mssum(0.), m_mssumabs(0.), m_mssumsqr(0.), m_mssumenh(0.), m_mssumenhabs(0.), m_mssumenhsqr(0.),
  m_msn(0.), m_msncut(0.), m_sn(0), m_sncut(0), m_son(1),
  m_external_selectionweight(-1), // -1 is the sentinel meaning "not set"
  m_writeout(false),
  p_whisto_pos(NULL),
  p_whisto_neg(NULL)
{
  m_colorscheme=cls::sum;
  m_helicityscheme=hls::sum;
}

bool Process_Integrator::Initialize
(BEAM::Beam_Spectra_Handler *const beamhandler,
 PDF::ISR_Handler *const isrhandler,
 YFS::YFS_Handler *const yfshandler)
{
  Settings& s = Settings::GetMainSettings();
  m_nin=p_proc->NIn();
  m_nout=p_proc->NOut();
  p_momenta.resize(m_nin+m_nout);
  p_beamhandler=beamhandler;
  p_isrhandler=isrhandler;
  p_yfshandler=yfshandler;
  m_swmode=s["SELECTION_WEIGHT_MODE"].SetDefault(0).Get<int>();
  static bool minit(false);
  if (!minit) {
    // weight histo bin number
    s_whbins = s["IB_WHBINS"].SetDefault(1300).Get<int>();
    minit=true;
  }
  return true;
}

Process_Integrator::~Process_Integrator()
{
  if (p_whisto_pos!=NULL) delete p_whisto_pos;
  if (p_whisto_neg!=NULL) delete p_whisto_neg;
}

double Process_Integrator::SelectionWeight(const int mode) const
{
  if (!p_proc->IsGroup()) {
    if (mode!=0) {
      if (m_external_selectionweight != -1) {
	return m_external_selectionweight;
      }
      // m_selweight_xs == 0 when m_ssumenhabs == 0 (no events passed cuts or all weights
      // are zero). Using the pre-computed m_selweight_xs instead of evaluating
      // TotalResult()*m_meanenhfunc/sqrt(m_effevperev) at call time avoids 0/0 for
      // processes where positive and negative weights cancel (zero total cross section
      // but non-zero absolute cross section, e.g. interference terms).
      if (m_effi < 0.) {
        THROW(critical_error, METHOD+"(): m_effi="+to_string(m_effi)+" (uninitialized sentinel) for '"
	      +p_proc->Name()+"'; SetUpEnhance may not have been called.");
      }
      if (m_effi == 0. || m_selweight_xs == 0.) return 0.;
      return m_selweight_xs * m_enhancefac / m_effi;
    }
    // Return -1 for completely unsampled processes (no integration data yet).
    // Via dabs(-1)=1 in Process_Group::OneEvent this gives them a uniform
    // placeholder weight so they still get sampled during optimisation.
    // Returning 0 instead would permanently exclude them from selection.
    if (m_n+m_sn==0.0) return -1.0;
    // Process has been sampled but has zero cross section: exclude from selection.
    if (m_totalxs==0.0) return 0.0;
    double selweight = m_swmode==0 ?
      sqrt((m_n+m_sn-1) * sqr(TotalVar()) + sqr(TotalResult())) :
      dabs(m_totalxs);
    return selweight*m_enhancefac;
  }
  double sw(0.0);
  for (size_t i(0);i<p_proc->Size();++i) {
    // dabs() is intentional: for mode==0 an unsampled child returns -1 as a placeholder
    // (see comment above), and dabs(-1)=1 gives it a uniform weight in the group sum.
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
    if (msg_LevelIsTracking()) DO_STACK_TRACE;
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
  // Simple averaging mode: mean over all events accumulated so far across all cycles and the current step.
  if (m_smode==0) return m_n+m_sn?(m_totalsum+m_ssum)/(m_n+m_sn):0.0; // 0 if no events sampled yet
  // Variance-weighted mode, but no historical data yet (first cycle): use current step only.
  if (m_ssigma2==0.0) return m_sn?m_ssum/m_sn:0.0; // 0 if current step is also empty
  // Current step has fewer than 2 events, so Sigma2 is undefined: fall back to history-only estimate.
  if (m_sn<2) return m_ssigma2?m_totalsum/m_ssigma2:0.0; // 0 if there is no history either
  // General case: precision-weighted combination of historical result (weight m_ssigma2)
  // and current step's result (weight s2=Sigma2). Both weights zero means no data at all.
  double s2(Sigma2());
  return m_ssigma2+s2?(m_totalsum+s2*m_ssum/m_sn)/(m_ssigma2+s2):0.0;
}

std::vector<double> Process_Integrator::TotalEffiAndEffEvPerEv(bool unweighted) const
{
  std::vector<double> totaleffiandeffevperev(2,-1.0);
  if (p_proc->IsGroup()) {
    //weight effevperev with selection weight
    double sum_selw=0.0;
    double sum_effi=0.0;
    double sum_xsec=0.0;
    double sum_xsec_abs=0.0;
    for (size_t i(0);i<p_proc->Size();++i) {
      //need to sum up weighted with wsel
      std::vector<double> proci_totaleffiandeffevperev = (*p_proc)[i]->Integrator()->TotalEffiAndEffEvPerEv(unweighted);
      double proci_effevperev = proci_totaleffiandeffevperev[1];
      if (proci_effevperev==-1.0) {
        msg_Error() << "WARNING:"<<METHOD<<"(): Subprocess "<<(*p_proc)[i]->Name()
                     <<" has undefined TotalEffiAndEffEvPerEv. SetUpEnhance might not have been called.\n";
        continue;
      }
      double proci_effi = proci_totaleffiandeffevperev[0];
      // A subprocess with zero efficiency or zero effective events per event contributes
      // nothing to the selection-weight-weighted average and its proci_selw is undefined.
      if (proci_effevperev == 0. || proci_effi == 0.) continue;
      double proci_xsec = (*p_proc)[i]->Integrator()->GetSSumEnh()/m_sn;
      double proci_selw = dabs(proci_xsec)/sqrt(proci_effevperev)/proci_effi;

      sum_xsec += proci_xsec;
      sum_xsec_abs += dabs(proci_xsec);
      sum_selw += proci_selw;
      sum_effi+=proci_effi*proci_selw;
      //this is from whisto, but want to be more correct by not relying on bin width approximation
      //sum_selw += (*p_proc)[i]->Integrator()->SelectionWeight(wmode);
      //sum_effi+=(*p_proc)[i]->Integrator()->Efficiency()*(*p_proc)[i]->Integrator()->SelectionWeight(wmode);
    }
    if (sum_selw!=0) {
      totaleffiandeffevperev[0] = sum_effi/sum_selw;
      totaleffiandeffevperev[1] = pow(sum_xsec_abs/sum_effi,2)*pow(sum_xsec/sum_xsec_abs,2);
    }
    return totaleffiandeffevperev;
  }
  if (m_sncut == 0) return totaleffiandeffevperev;
  double meanweight = m_ssumenh/m_sncut;
  double meanabsweight = m_ssumenhabs/m_sncut;
  double meansqrweight = m_ssumenhsqr/m_sncut;
  if (unweighted) {
    totaleffiandeffevperev[0] = dabs(meanabsweight)/m_smax;
    totaleffiandeffevperev[1] = pow(meanweight,2)/pow(meanabsweight,2);
  } else {
    totaleffiandeffevperev[0] = 1;
    totaleffiandeffevperev[1] = pow(meanweight,2)/meansqrweight;
  }
  return totaleffiandeffevperev;
}

double Process_Integrator::TotalVar() const
{
  if (m_nin==1 && m_nout==2) return 0.;
  if (m_smode==0) {
    if (m_n+m_sn<2) return TotalResult();
    return sqrt(dabs((m_totalsumsqr+m_ssumsqr)/(m_n+m_sn)
		 -sqr((m_totalsum+m_ssum)/(m_n+m_sn)))/(m_n+m_sn-1.0));
  }
  double s2(m_totalsumsqr);
  if (m_sn>1) {
    double vij2((m_sn-1)/
		dabs(m_ssumsqr/m_sn-sqr(m_ssum/m_sn)));
    s2+=sqr(Sigma2())/vij2;
  }
  if (s2<0.0) return 0.0;
  return sqrt(s2)/TotalSigma2();
}

void Process_Integrator::OptimizeSubResult(const double &s2)
{
  //called after each optimisation cycle
  m_n+=m_sn;
  //only m_totalsum and m_totalsumsqr are needed accumulated over all optimsation cycles to get the final cross section result and error estimate.
  if (m_smode==0) {
    m_totalsum+=m_ssum;
    m_totalsumsqr+=m_ssumsqr;
  }
  else {
    double vij2((m_sn-1)/
		dabs(m_ssumsqr/m_sn-sqr(m_ssum/m_sn)));
    m_ssigma2+=s2; 
    m_totalsum+=s2*m_ssum/m_sn;
    m_totalsumsqr+=sqr(s2)/vij2;
  }
  m_ssum=m_ssumenh=m_ssumenhabs=m_ssumsqr=m_ssumenhsqr=0.0;
  m_sn=0;
  m_sncut=0;
  if (p_colint!=NULL) p_colint->Optimize();
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->OptimizeSubResult(s2);
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
  if (p_whisto_pos) {
    delete p_whisto_pos;
    p_whisto_pos=0;
  }
  if (p_whisto_neg) {
    delete p_whisto_neg;
    p_whisto_neg=0;
  }

  double av(dabs(TotalResult()));
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
  p_whisto_pos = new Histogram(10,av*1.e-100,av*1.e30,s_whbins);
  p_whisto_neg = new Histogram(10,av*1.e-100,av*1.e30,s_whbins);
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->InitWeightHistogram();
}

bool Process_Integrator::ReadInXSecs(const std::string &path)
{
  std::string fname(p_proc->ResultsName());
  std::string integration_sherpa_version(p_pshandler->IntegrationSherpaVersion());
  size_t vn;
  std::string name, dummy;
  My_In_File from(path+"/"+fname);
  if (!from.Open()) return false;
  from->precision(16);
  if (integration_sherpa_version=="before 3.1") {
    *from>>name>>m_totalxs>>m_weightmax>>m_totalerr>>m_totalsum>>m_totalsumsqr
	 >>m_n>>m_ssum>>m_ssumsqr>>m_smax>>m_ssigma2>>m_sn>>m_wmin
	 >>m_son>>dummy>>dummy>>vn;
    int wmode = ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE"));
    if (wmode!=0) msg_Error() << "WARNING: Exact information about enhancement is missing for: "<< p_proc->ResultsName() << ". Please use Sherpa 3.1 or newer to consider them. Assuming that there is no enhancement during sampling. There is no ehancement backward-compatibility implemented yet." << std::endl;
    m_ssumenh = m_ssum;
    m_ssumenhsqr = m_ssumsqr;
    m_ssumenhabs = dabs(m_ssumenh);
    msg_Tracking() << "Note: Information about cut efficiency missing for: "<< p_proc->ResultsName() << ". Please use Sherpa 3.1 or newer to separate it. Exp. eff. includes cut efficiency." << std::endl;
    m_sncut = m_sn;
  } else {
    *from>>name>>m_totalxs>>m_weightmax>>m_totalerr>>m_totalsum>>m_totalsumsqr
	 >>m_n>>m_ssum>>m_ssumenh>>m_ssumenhabs>>m_ssumsqr>>m_ssumenhsqr>>m_smax>>m_ssigma2>>m_sn>>m_sncut>>m_wmin
	 >>m_son>>dummy>>dummy>>vn;
  }
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
		<<m_weightmax*rpa->Picobarn()<<std::endl;
  if (!p_proc->ReadIn(path)) return false;
  if (p_colint!=NULL) p_colint->ReadIn(path+"/"+fname+"_Color");
  bool res(true);
  if (p_proc->IsGroup()) {
    for (size_t i(0);i<p_proc->Size();++i)
      if (!(*p_proc)[i]->Integrator()->ReadInXSecs(path)) res=false;
    
    //consistency check of m_totalxs
    double sum(0.0);
    for (size_t i(0);i<p_proc->Size();++i) {
      sum+=(*p_proc)[i]->Integrator()->TotalXS();
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
  return res;
}

void Process_Integrator::ReadInHistogram(std::string dir)
{
  std::string integration_sherpa_version(p_pshandler->IntegrationSherpaVersion());
  std::string filename = dir+"/"+p_proc->ResultsName();
  if (integration_sherpa_version=="before 3.1") {
    int wmode = ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE"));
    if (wmode!=0) {//Todo: also check for Timing_statistics here
      msg_Error() << "WARNING: Information of the weight sign distribution from integration are missing for: "<< p_proc->ResultsName() << ". Please use Sherpa 3.1 or newer to consider them. Assuming that all weights have the same sign - unweighted sampling might have an inefficient process selection. Also all estimates regarding the number of effective events are likely off." << std::endl;
      msg_Error() << "WARNING: Information for partial unweighting might be incomplete for: "<< p_proc->ResultsName() << ". Please use Sherpa 3.1 or newer to be able to choose any max_epsilon. Partial unweighting has no effect above a certain point." << std::endl;
    }
    if (!FileExists(filename)) return;
    if (p_whisto_pos) delete p_whisto_pos;
    if (p_whisto_neg) delete p_whisto_neg;
    p_whisto_pos = new Histogram(filename);
    p_whisto_neg = new Histogram(p_whisto_pos);
    if (TotalResult()>0)
      p_whisto_neg->Reset();
    else
      p_whisto_pos->Reset();
  } else {
    if (!FileExists(filename+"_pos")) return;
    if (p_whisto_pos) delete p_whisto_pos;
    if (p_whisto_neg) delete p_whisto_neg;
    p_whisto_pos = new Histogram(filename+"_pos");
    p_whisto_neg = new Histogram(filename+"_neg");
  }
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
  *outfile<<fname<<"  "<<m_totalxs<<"  "<<m_weightmax<<"  "
	 <<m_totalerr<<" "<<m_totalsum<<" "<<m_totalsumsqr<<" "
	 <<m_n<<" "<<m_ssum<<" "<<m_ssumenh<<" "<<m_ssumenhabs<<" "<<m_ssumsqr<<" "<<m_ssumenhsqr<<" "<<m_smax<<" "
	 <<m_ssigma2<<" "<<m_sn<<" "<<m_sncut<<" "<<m_wmin<<" "<<m_son<<" "
	 <<-1<<" "<<-1<<"\n"<<m_vsn.size()<<"\n";
  for (size_t i(0);i<m_vsn.size();++i)
    *outfile<<m_vsmax[i]<<" "<<m_vsum[i]<<" "
	   <<m_vsn[i]<<" "<<-1<<"\n";
  p_proc->WriteOut(path);
  if (p_colint!=NULL) p_colint->WriteOut(path+"/"+fname+"_Color");
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->WriteOutXSecs(path);
}

void Process_Integrator::WriteOutHistogram(std::string dir)
{
  if (p_whisto_pos) p_whisto_pos->Output(dir+"/"+p_proc->ResultsName()+"_pos");
  if (p_whisto_neg) p_whisto_neg->Output(dir+"/"+p_proc->ResultsName()+"_neg");
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->WriteOutHistogram(dir);
}

void Process_Integrator::SetTotal(const int mode)  
{
  if (p_proc->IsGroup()) {
    m_weightmax=0.0;
    msg_Indent();
    for (size_t i(0);i<p_proc->Size();++i) {
      (*p_proc)[i]->Integrator()->SetTotal(msg_LevelIsTracking());
      m_weightmax+=(*p_proc)[i]->Integrator()->WeightMax();
    }
  }
  double totalxs=TotalResult();
  double totalerr=TotalVar();
  if (m_totalxs==0.0) Reset(0);
  m_totalxs=totalxs;
  m_totalerr=totalerr;
  if (mode && m_totalxs!=0.0) {
    if (p_proc->NIn()==1) {
      msg_Info()<<om::bold<<p_proc->ResultsName()<<om::reset<<" : "<<om::blue<<om::bold
                <<m_totalxs<<" GeV"<<om::reset<<" +- ( "<<om::red
                <<m_totalerr<<" GeV = "<<dabs(m_totalerr/m_totalxs)*100.;
    }
    else {
      msg_Info()<<om::bold<<p_proc->ResultsName()<<om::reset<<" : "<<om::blue<<om::bold
                <<m_totalxs*rpa->Picobarn()<<" pb"<<om::reset<<" +- ( "<<om::red
                <<m_totalerr*rpa->Picobarn()<<" pb = "<<dabs(m_totalerr/m_totalxs)*100.;
    }
    std::vector<double> meaneffiandeffevperev = TotalEffiAndEffEvPerEv(1);
    msg_Info()<<" %"<<om::reset<<" ) "<<om::bold<<" exp. eff: "
	      <<om::red<<(100.*meaneffiandeffevperev[0])<<" % "<<om::reset
	      <<om::bold<<" Neff/N: "<<om::brown<<meaneffiandeffevperev[1]<<om::reset<<std::endl;
  }
}

void Process_Integrator::SetFullUnweightingStats(const double w_max)
{
  // m_effi is the unweighting efficiency for events which already passed cuts.
  // w_max==0 with m_sncut>0 is physically impossible (non-zero-weight events imply
  // m_weightmax>0), so it signals a broken upstream max-weight calculation.
  if (w_max == 0. && m_sncut > 0) {
    THROW(critical_error, METHOD+"(): w_max=0 but m_sncut="+to_string(m_sncut)+" for '"
	  +p_proc->Name()+"'; upstream max-weight may be wrong.");
  }
  m_effi         = (m_sncut > 0 && w_max != 0.) ? m_ssumenhabs/m_sncut/w_max : 0.;
  m_effevperev   = m_ssumenhabs != 0. ? pow(m_ssumenh/m_ssumenhabs,2) : 0.;
  m_selweight_xs = m_sn > 0 ? m_ssumenhabs/m_sn : 0.;
}

double Process_Integrator::GetMaxEps(double epsilon)
{
  if (!p_whisto_pos) {
    // No histogram: fall back to full unweighting against m_weightmax.
    SetFullUnweightingStats(m_weightmax);
    return m_weightmax;
  }
  //construct whisto having all pos and neg weights
  std::unique_ptr<Histogram> p_whisto(new Histogram(p_whisto_pos));
  *p_whisto += *p_whisto_neg;
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
    // Stat-based max estimate; fall back to full unweighting against it.
    double fin_w_max = maxs[maxs.size()/2];
    SetFullUnweightingStats(fin_w_max);
    return fin_w_max;
  }

  //save last, first filled bin and more for next loop
  double whisto_abs_sum = 0;
  double whisto_sum = 0;
  double whisto_sum2 = 0;
  //number of fills differs from ->Fills(), because the first does not include over and underflow events (ratio is cut effi, if all w are in whisto range)
  int whisto_fills = 0;
  int last_filled_bin = p_whisto->Nbin()+2;
  int first_filled_bin = 0;
  //add it forward, because have to add small numbers first for numerical stability
  for (int i=1;i<p_whisto->Nbin()+2;i++) {
    if (p_whisto->Value(i)!=0) {
      //bin middle times count is added
      double w = exp(log(10.)*(p_whisto->Xmin()+(i-0.5)*p_whisto->BinSize()));
      if (exp(log(10.)*(p_whisto->Xmin()+i*p_whisto->BinSize())) > m_weightmax) w = m_weightmax;
      whisto_abs_sum += p_whisto->Value(i) * w;
      whisto_sum2 += p_whisto->Value(i) * w * w;
      whisto_sum += p_whisto_pos->Value(i) * w;
      whisto_sum -= p_whisto_neg->Value(i) * w;
      whisto_fills += p_whisto->Value(i);
      last_filled_bin = i;
      if (first_filled_bin==0) {
	first_filled_bin = i;
      }
    }
  }

  if (first_filled_bin<10) msg_Error() << "WARNING: The lower bin edge of whisto might be to high for subprocess " << p_proc->ResultsName() << " to cover all weights!" << std::endl;
  if (last_filled_bin>p_whisto->Nbin()-10) msg_Error() << "WARNING: The upper bin edge of whisto might be to low for subprocess " << p_proc->ResultsName() << " to cover all weights!" << std::endl;

  if (whisto_abs_sum==0) {
    //fallback for empty WD histogram: fall back to full unweighting against m_weightmax.
    int wmode = ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE"));
    if (wmode!=0) msg_Error()<<"Warning: The weight histogram for "<< p_proc->ResultsName() << " is empty. Partial unweighting is not possible."<<std::endl;
    SetFullUnweightingStats(m_weightmax);
    return m_weightmax;
  }

  //cutxs: sum of |w| which are cut atm
  double cutxs_abs = 0.;
  //cutxs: sum of w which are cut atm
  double cutxs = 0.;
  //cutxs: sum of w^2 which are cut atm
  double cutxs2 = 0.;
  //cnt: number of events cut atm
  double cnt = 0.;

  double pxs = whisto_abs_sum*(1-epsilon);
  for (int i=first_filled_bin-1;i<last_filled_bin+1;i++) {
    //bin middle times count is added
    double w = exp(log(10.)*(p_whisto->Xmin()+(i-0.5)*p_whisto->BinSize()));
    if (exp(log(10.)*(p_whisto->Xmin()+i*p_whisto->BinSize())) > m_weightmax) w = m_weightmax;
    cutxs_abs+= p_whisto->Value(i) * w;
    cutxs+= p_whisto_pos->Value(i) * w;
    cutxs-= p_whisto_neg->Value(i) * w;
    cutxs2+= p_whisto->Value(i) * pow(w,2);
    cnt+= p_whisto->Value(i);
    if (cutxs_abs>=pxs) {
      double fin_w_max = Min(exp(log(10.)*(p_whisto->Xmin()+i*p_whisto->BinSize())),dabs(m_weightmax));
      double mean_efficiency = (whisto_fills-cnt+cutxs_abs/fin_w_max)/whisto_fills;
      // selweight_num2 = (m_selweight_xs * total_fills)^2 = squared numerator of m_selweight_xs;
      // equivalently the denominator of m_effevperev (see equations below)
      double selweight_num2 = (cutxs_abs * fin_w_max + (whisto_sum2 - cutxs2))
                              * mean_efficiency * (double)whisto_fills;
      //without enhancement function: whisto_sum/p_whisto->Fills() = dabs(TotalResult())
      //with enhancement function:    whisto_sum/p_whisto->Fills() = dabs(TotalResult())*m_meanenhfunc
      //defined such that: dabs(TotalResult())*m_meanenhfunc = mean_w*cut_effi = whisto_sum/whisto_fills * (whisto_fills/p_whisto->Fills())
      //m_meanenhfunc = whisto_sum/p_whisto->Fills()/dabs(TotalResult());
      m_effi = mean_efficiency;
      //according to equation (39) in arXiv:2506.06203v2
      m_effevperev = pow(whisto_sum, 2) / selweight_num2;
      // m_selweight_xs is algebraically equivalent to |TotalResult()*m_meanenhfunc/sqrt(m_effevperev)|
      // when there is no sign cancellation.
      // Sign-cancel-safe: selweight_num2 uses only absolute/squared weight sums,
      // not the signed whisto_sum, so m_selweight_xs stays finite when positive
      // and negative weights cancel.
      double total_fills = p_whisto->Fills();
      m_selweight_xs = (total_fills > 0. && selweight_num2 > 0.)
          ? sqrt(selweight_num2) / total_fills
          : 0.;
      return fin_w_max;
    }
  }
  // Loop ran to completion: all weights fit under m_weightmax, full unweighting.
  SetFullUnweightingStats(m_weightmax);
  return m_weightmax;
}

void Process_Integrator::SetUpEnhance(const int omode) 
{
  if (p_proc->IsGroup()) {
    //this is the only place where m_max is calculated for groups
    m_max=0.0;
    m_effi=0.0;
    m_effevperev=0.0;
    double effevperev_signed=0.0;
    int wmode = ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE"));
    for (size_t i(0);i<p_proc->Size();++i) {
      (*p_proc)[i]->Integrator()->SetUpEnhance(msg_LevelIsTracking());
      m_max+=(*p_proc)[i]->Integrator()->Max();
      //need to sum up weighted with wsel: m_effi
      m_effi+=(*p_proc)[i]->Integrator()->Efficiency()*(*p_proc)[i]->Integrator()->SelectionWeight(wmode);
      //interested in average effevperev after unweighting: multiply with efficiency
      double child_effevperev = (*p_proc)[i]->Integrator()->EffEvPerEv();
      // After SetUpEnhance, EffEvPerEv is always >= 0: zero-XS processes give 0,
      // processes with data give pow(ssumenh/ssumenhabs,2) > 0. A value of -1 is the
      // uninitialized sentinel and means SetUpEnhance was not called on this child.
      if (child_effevperev >= 0.) {
        double effevperev = sqrt(child_effevperev)*(*p_proc)[i]->Integrator()->Efficiency()*(*p_proc)[i]->Integrator()->SelectionWeight(wmode);
        m_effevperev+=effevperev;
        effevperev_signed+=effevperev*((*p_proc)[i]->Integrator()->TotalResult() >= 0. ? 1. : -1.);
      } else {
        THROW(critical_error, METHOD+"(): EffEvPerEv="+to_string(child_effevperev)
                   +" (uninitialized) for '"+(*p_proc)[i]->Name()
	      +"'; SetUpEnhance may not have been called.");
      }
    }
    // If every child subprocess has zero efficiency or zero weight, the group as a
    // whole has no effective contribution and m_effevperev is meaningless and set to 0 for 0 selection weight.
    if (m_effi != 0. && m_effevperev != 0.)
      m_effevperev = pow(m_effevperev/m_effi,2)*pow(effevperev_signed/m_effevperev,2);
    else
      m_effevperev = 0.;
    double sel_wgt = SelectionWeight(wmode);
    // sel_wgt is the sum of all children's selection weights. It is zero only
    // when every child has zero cross section, in which case m_effi is also 0.
    m_effi = sel_wgt != 0. ? m_effi/sel_wgt : 0.;
    if (omode || p_proc->Parent()==p_proc)
      if (p_whisto_pos)
    msg_Info()<<"  reduce max for "<<p_proc->ResultsName()<<" to "
	      <<m_max/m_weightmax<<" ( eps = "<<m_maxeps<<" -> exp. eff "
	      << m_effi<<", Neff/N "
	      << m_effevperev << " ) "<<std::endl;
    return;
  }
  // m_meanenhfunc = mean_enhanced_weight / mean_unenhanced_weight.
  // mean_enhanced_weight = m_ssumenh/m_sn, mean_unenhanced_weight = TotalResult()
  // m_sn is used instead of m_sncut, to have an enhanced cross section compatible with TotalResult()
  m_meanenhfunc = TotalResult() ? m_ssumenh/m_sn/dabs(TotalResult()) : 1;
  if (m_maxeps!=0.0) {
    double max(GetMaxEps(m_maxeps));
    if (omode || msg->LevelIsTracking()) {
      msg_Info()<<"  reduce max for "<<p_proc->ResultsName()<<" to "
                <<max/m_weightmax<<" ( eps = "<<m_maxeps<<" -> exp. eff "
                <<  m_effi <<", Neff/N " << m_effevperev << " ) "<<std::endl;
    }
    SetMax(max);
  } else {
    //should rely on m_weightmax here -> to be compatible with partially unweighting the m_effi etc needs to be set
    SetMax(m_weightmax);
    SetFullUnweightingStats(m_max);
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
  if (dabs(value)!=0) m_msncut++;
  m_mssumenh    += value;
  m_mssumenhabs += dabs(value);
  m_mssumenhsqr += sqr(value);
  m_mssum    += value/enhance;
  m_mssumabs += dabs(value/enhance);
  m_mssumsqr += sqr(value/enhance);
#else
  m_sn++;
  if (dabs(value)!=0) m_sncut++;
  m_ssumenh    += value;
  m_ssumenhabs += dabs(value);
  m_ssumenhsqr += sqr(value);
  m_ssum    += value/enhance;
  m_ssumsqr += sqr(value/enhance);
#endif
  double max=dabs(value)/dabs(p_proc->Last())*
    ATOOLS::Max(p_proc->LastPlus(),-p_proc->LastMinus());
  if (max>m_weightmax)  m_weightmax  = max;
  if (max>m_smax) m_smax = max;
  if (p_whisto_pos) {
    if(value>0.)  p_whisto_pos->Insert(max);
    else if(value<0.)  p_whisto_neg->Insert(max);
    else p_whisto_pos->Insert(1.0,0.0);
  }
  if (p_colint!=NULL) p_colint->AddPoint(value);
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
  //sets m_max for unweighting
  m_max=max;
} 

void Process_Integrator::Reset(const int mode)
{
  m_n=0;
  m_totalxs=m_totalsum=m_totalsumsqr=m_totalerr=0.0;
  m_smax=m_weightmax=m_max=m_wmin=m_ssigma2=m_ssumsqr=m_ssum=m_ssumenh=m_ssumenhabs=m_ssumenhsqr=0.0;
  m_sn=0;
  m_sncut=0;
  m_son=1;
  m_vsmax.clear(); 
  m_vsn.clear();   
  m_vsum.clear(); 
  if (p_proc->IsGroup() && mode==1)
    for (size_t i(0);i<p_proc->Size();++i) 
      (*p_proc)[i]->Integrator()->Reset(mode);
}

void Process_Integrator::ResetMax(const int flag)
{
  if (p_proc->IsGroup()) {
    m_weightmax=0.0;
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->ResetMax(flag);
    return;
  }
  if (flag==3) {
    m_vsmax.clear();
    m_vsn.clear();
    m_vsum.clear();
    m_weightmax=0.0;
    return;
  }
  if (flag==0) {
    if (m_vsmax.size()>1) {
      m_vsmax.erase(m_vsmax.begin());
      m_vsn.erase(m_vsn.begin());
      m_vsum.erase(m_vsum.begin());
    }
    if (m_vsmax.empty()) {
      m_vsmax.push_back(m_weightmax);
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
  m_weightmax=0.0;
  if (flag==1) return; //restart m_weightmax determination after "full optimisation"
  for (size_t i=0;i<m_vsmax.size();i++) {
    m_weightmax=ATOOLS::Max(m_weightmax,m_vsmax[i]);
  }
} 

void Process_Integrator::SetPSHandler(Phase_Space_Handler * pshandler) {
  p_pshandler=pshandler;
  if (p_proc->IsGroup())
    for (size_t i(0);i<p_proc->Size();++i)
      (*p_proc)[i]->Integrator()->SetPSHandler(pshandler);
}

void Process_Integrator::SetPSHandler(const double &maxerr,const std::string eobs,const std::string efunc) {
  p_ownpshandler.reset(new Phase_Space_Handler(this, maxerr, eobs, efunc));
  SetPSHandler(p_ownpshandler.get());
}

void Process_Integrator::MPICollect
(std::vector<double> &sv,std::vector<double> &mv,size_t &i)
{
  sv.resize(8*(i+1));
  mv.resize(2*(i+1));
  sv[8*i+0]=m_msn;
  sv[8*i+1]=m_msncut;
  sv[8*i+2]=m_mssum;
  sv[8*i+3]=m_mssumabs;
  sv[8*i+4]=m_mssumenh;
  sv[8*i+5]=m_mssumenhabs;
  sv[8*i+6]=m_mssumenhsqr;
  sv[8*i+7]=m_mssumsqr;
  mv[2*i+0]=m_weightmax;
  mv[2*i+1]=m_smax;
  ++i;
  if (p_proc->IsGroup())
    for (size_t j(0);j<p_proc->Size();++j)
      (*p_proc)[j]->Integrator()->MPICollect(sv,mv,i);
}

void Process_Integrator::MPIReturn
(std::vector<double> &sv,std::vector<double> &mv,size_t &i)
{
  m_msn=sv[8*i+0];
  m_msncut=sv[8*i+1];
  m_mssum=sv[8*i+2];
  m_mssumabs=sv[8*i+3];
  m_mssumenh=sv[8*i+4];
  m_mssumenhabs=sv[8*i+5];
  m_mssumenhsqr=sv[8*i+6];
  m_mssumsqr=sv[8*i+7];
  m_weightmax=mv[2*i+0];
  m_smax=mv[2*i+1];
  ++i;
  if (p_proc->IsGroup())
    for (size_t j(0);j<p_proc->Size();++j)
      (*p_proc)[j]->Integrator()->MPIReturn(sv,mv,i);
}

void Process_Integrator::MPISync(const int mode)
{
  if (p_whisto_pos) p_whisto_pos->MPISync();
  if (p_whisto_neg) p_whisto_neg->MPISync();
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
  m_sncut+=m_msncut;
  m_ssumenh+=m_mssumenh;
  m_ssumenhabs+=m_mssumenhabs;
  m_ssumenhsqr+=m_mssumenhsqr;
  m_ssum+=m_mssum;
  m_ssumsqr+=m_mssumsqr;
  m_msn=m_msncut=m_mssumenh=m_mssumenhabs=m_mssumenhsqr=m_mssum=m_mssumabs=m_mssumsqr=0.0;
#endif
  if (p_colint!=NULL) p_colint->MPISync();
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
  if (!FileExists(m_resultpath+".zip")) return;
  if (!Copy(m_resultpath+".zip",m_resultpath+".zip~",true))
    msg_Error()<<METHOD<<"(): Copy error. "
	       <<strerror(errno)<<"."<<std::endl;
}

void Process_Integrator::StoreResults(const int mode)
{
  if (m_msn) MPISync();
  if (m_resultpath.length()==0) return;
  if (m_totalxs!=0.0 && mode==0) return;
  SetTotal(0);
  std::string fname(p_proc->ResultsName());
  WriteOutXSecs(m_resultpath+"/"+p_proc->Generator()->Name()+"/XS_"+fname);
  WriteOutHistogram(m_resultpath+"/"+p_proc->Generator()->Name()+"/WD_"+fname);
  p_pshandler->WriteOut(m_resultpath+"/"+p_proc->Generator()->Name()+"/MC_"+fname);
  My_In_File::CloseDB(m_resultpath+"/",0);
  StoreBackupResults();
}

void Process_Integrator::ReadResults()
{
  if (m_resultpath.length()==0) return;
  std::string fname(p_proc->ResultsName());
  //needs to be called before the next three, to set p_pshandler->IntegrationSherpaVersion())
  p_pshandler->ReadInProperties(m_resultpath+"/"+p_proc->Generator()->Name()+"/MC_"+fname);
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

