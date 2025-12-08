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
  m_n(0), m_itmin(0), m_itmax(1000000), m_max(0.), m_totalxs(0.), 
  m_totalsum (0.), m_totalsumsqr(0.), m_totalerr(0.), m_ssum(0.), 
  m_ssumsqr(0.), m_smax(0.), m_ssigma2(0.), m_wmin(0.),
  m_mssum(0.), m_mssumsqr(0.), m_msn(0.), m_sn(0), m_son(1),
  m_external_selectionweight(-1),
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
  m_timing_statistics_large_weight_fraction=s["TIMING_STATISTICS_LARGE_WEIGHT_FRACTION"].SetDefault(0.001).Get<double>();
  m_timing_statistics=s["TIMING_STATISTICS"].SetDefault(0).Get<int>();
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
      //Sherpa manual - only works with completely written out whisto - otherweise m_meanenhfunc ("mean enhancement function") wrong
      return dabs(TotalResult()*m_meanenhfunc*m_enhancefac/m_effi/m_effevperev);
    }
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
  m_n+=m_sn;
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
  m_ssum=m_ssumsqr=0.0;
  m_sn=0;
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
  if (!FileExists(filename+"_pos")) return;
  if (p_whisto_pos) delete p_whisto_pos;
  if (p_whisto_neg) delete p_whisto_neg;
  p_whisto_pos = new Histogram(filename+"_pos");
  p_whisto_neg = new Histogram(filename+"_neg");
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
	 <<m_totalerr<<" "<<m_totalsum<<" "<<m_totalsumsqr<<" "
	 <<m_n<<" "<<m_ssum<<" "<<m_ssumsqr<<" "<<m_smax<<" "
	 <<m_ssigma2<<" "<<m_sn<<" "<<m_wmin<<" "<<m_son<<" "
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
    m_max=0.0;
    msg_Indent();
    for (size_t i(0);i<p_proc->Size();++i) {
      (*p_proc)[i]->Integrator()->SetTotal(msg_LevelIsTracking());
      m_max+=(*p_proc)[i]->Integrator()->Max();
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
                <<m_totalerr<<" GeV = "<<dabs(m_totalerr/m_totalxs)*100.
                <<" %"<<om::reset<<" ) "<<om::bold<<" exp. eff: "
                <<om::red<<(100.*dabs(m_totalxs/m_max))<<" %"<<om::reset<<std::endl;
    }
    else {
      msg_Info()<<om::bold<<p_proc->ResultsName()<<om::reset<<" : "<<om::blue<<om::bold
                <<m_totalxs*rpa->Picobarn()<<" pb"<<om::reset<<" +- ( "<<om::red
                <<m_totalerr*rpa->Picobarn()<<" pb = "<<dabs(m_totalerr/m_totalxs)*100.
                <<" %"<<om::reset<<" ) "<<om::bold<<" exp. eff: "
                <<om::red<<(100.*dabs(m_totalxs/m_max))<<" %"<<om::reset<<std::endl;
    }
  }
}

double Process_Integrator::GetMaxEps(double epsilon)
{
  if (!p_whisto_pos) return m_max;
  //construct whisto having all pos and neg weights
  ATOOLS::Histogram *p_whisto = new Histogram(p_whisto_pos);
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
    return maxs[maxs.size()/2];
  }

  //save last, first filled bin and more for next loop
  double whisto_abs_sum = 0;
  double whisto_sum = 0;
  //number of fills differs from ->Fills(), because the first does not include over and underflow events (ratio is cut effi, if all w are in whisto range)
  int whisto_fills = 0;
  int last_filled_bin = p_whisto->Nbin()+2;
  int first_filled_bin = 0;
  for (int i=p_whisto->Nbin()+1;i>0;i--) {
    if (p_whisto->Value(i)!=0) {
      //bin middle times count is added
      double w = exp(log(10.)*(p_whisto->Xmin()+(i-0.5)*p_whisto->BinSize()));
      whisto_abs_sum += p_whisto->Value(i) * w;
      whisto_sum += p_whisto_pos->Value(i) * w;
      whisto_sum -= p_whisto_neg->Value(i) * w;
      whisto_fills += p_whisto->Value(i);
      last_filled_bin = i;
      if (first_filled_bin==0) {
	first_filled_bin = i;
      }
    }
  }

  if (last_filled_bin<10) msg_Error() << "WARNING: The lower bin edge of whisto might be to high for subprocess " << p_proc->ResultsName() << " to cover all weights!" << std::endl;
  if (first_filled_bin>p_whisto->Nbin()-10) msg_Error() << "WARNING: The upper bin edge of whisto might be to low for subprocess " << p_proc->ResultsName() << " to cover all weights!" << std::endl;

  //cutxs: sum of |w| which are cut atm
  double cutxs_abs = 0.;
  //cutxs: sum of w which are cut atm
  double cutxs = 0.;
  //cutxs: sum of w^2 which are cut atm
  double cutxs2 = 0.;
  //cnt: number of events cut atm
  double cnt = 0.;

  if (m_timing_statistics) {
    //determine alpha and beta for target_fraction -> done in next loop. For its calculation whisto_sum, whisto_sum2 and whisto_fills are needed.
    double target_fraction = m_timing_statistics_large_weight_fraction;
    double whisto_sum2_fraction = 0;
    double whisto_sum_fraction = 0;
    double whisto_sum_abs_fraction = 0;
    int whisto_fills_fraction = 0;
    double fraction_pxs = target_fraction*whisto_sum;
    for (int i=first_filled_bin+1;i>last_filled_bin-1;i--) {
      if (p_whisto->Value(i)!=0) {
        double w = exp(log(10.)*(p_whisto->Xmin()+(i-0.5)*p_whisto->BinSize()));
        whisto_sum_abs_fraction += p_whisto->Value(i)*w;
        whisto_sum_fraction += p_whisto_pos->Value(i)*w;
        whisto_sum_fraction -= p_whisto_neg->Value(i)*w;
        whisto_sum2_fraction += p_whisto->Value(i)*w*w;
        whisto_fills_fraction += p_whisto->Value(i);
        //it is rounded up the target_fraction -> for fraction of 0 largest weight is picked
        if (whisto_sum_fraction>=fraction_pxs) break;
      }
    }


    std::vector<double> epsilon_max_values = rpa->gen.EpsilonValues();
    vector<double> wmax_manual_list(epsilon_max_values.size(),-1-1*(first_filled_bin==0));
    vector<double> alpha_manual_list(epsilon_max_values.size(),-1);
    vector<double> alpha_manual_fraction_list(epsilon_max_values.size(),-1);
    vector<double> efficiency_manual_list(epsilon_max_values.size(),-1);
    int curr_epsilon_max_manual_index = 0;
    double next_pxs_manual = whisto_sum*exp(log(10.)*epsilon_max_values[curr_epsilon_max_manual_index]);
    for (int i=first_filled_bin+1;i>last_filled_bin-1;i--) {
      double prev_cutxs2 = cutxs2;
      double prev_cutxs = cutxs;
      double prev_cutxs_abs = cutxs_abs;
      double prev_cnt = cnt;
      //bin middle times count is added
      double w = exp(log(10.)*(p_whisto->Xmin()+(i-0.5)*p_whisto->BinSize()));
      cutxs_abs+= p_whisto->Value(i) * w;
      cutxs+= p_whisto_pos->Value(i) * w;
      cutxs-= p_whisto_neg->Value(i) * w;
      cutxs2+= p_whisto->Value(i) * pow(w,2);
      cnt+= p_whisto->Value(i);
      //manual epsilon_max definition
      while (cutxs>=next_pxs_manual and (curr_epsilon_max_manual_index!=epsilon_max_values.size()-1 or last_filled_bin==i)) {
        double fin_w_max = Min(exp(log(10.)*(p_whisto->Xmin()+i*p_whisto->BinSize())),dabs(m_max));
        double this_cutxs2 = cutxs2;
        double this_cutxs = cutxs;
        double this_cutxs_abs = cutxs_abs;
        double this_cnt = cnt;
        if (cutxs==next_pxs_manual) {//needed for epsilon=1
          fin_w_max = Min(exp(log(10.)*(p_whisto->Xmin()+(i-1)*p_whisto->BinSize())),dabs(m_max));
        } else {
          this_cutxs2 = prev_cutxs2;
          this_cutxs = prev_cutxs;
          this_cutxs_abs = prev_cutxs_abs;
          this_cnt = prev_cnt;
        }
        //double mean_overweight = (1*(whisto_fills-this_cnt)+this_cutxs/fin_w_max)/whisto_fills;
        double mean_efficiency = (1*this_cnt+(whisto_abs_sum-this_cutxs_abs)/fin_w_max)/whisto_fills;
        msg_Debugging() << "Manual: " << epsilon_max_values[curr_epsilon_max_manual_index] << ": "<< mean_efficiency << std::endl;
        msg_Debugging() << "Manual: " << "this_cnt: "<< this_cnt << "whisto_fills: "<< whisto_fills<< "whisto_sum-this_cutxs: "<< whisto_sum-this_cutxs << "(whisto_sum-this_cutxs)/fin_w_max: "<< (whisto_sum-this_cutxs)/fin_w_max << std::endl;
        //according to equation (39) in arXiv:2506.06203v2
        double alpha = pow(whisto_sum/fin_w_max,2)/(((whisto_abs_sum-this_cutxs_abs+this_cutxs2/fin_w_max)/fin_w_max)*mean_efficiency*whisto_fills);
        double mean_efficiency_fraction = (1*this_cnt+(whisto_sum_abs_fraction-this_cutxs_abs)/fin_w_max)/whisto_fills_fraction;
        if (this_cutxs_abs>whisto_sum_abs_fraction) {
          mean_efficiency_fraction = 1;
        }
        double alpha_fraction = (whisto_abs_sum/whisto_fills/mean_efficiency)/(whisto_sum_abs_fraction/whisto_fills_fraction/mean_efficiency_fraction); //is beta
        msg_Debugging() << "Manual fraction beta: " << alpha_fraction  << std::endl;
        double temp_alpha_fraction = alpha_fraction;
        if (this_cutxs<whisto_sum_fraction) {//alpha is const else
          alpha_fraction = alpha_fraction*pow(whisto_sum_fraction/fin_w_max,2)/(((whisto_sum_abs_fraction-this_cutxs_abs+this_cutxs2/fin_w_max)/fin_w_max)*mean_efficiency_fraction*whisto_fills_fraction);
        } else {
          alpha_fraction = alpha_fraction*pow(whisto_sum_fraction,2)/(whisto_sum2_fraction*whisto_fills_fraction);
        }
        msg_Debugging() << "Manual fraction alpha: " << alpha_fraction/temp_alpha_fraction  << std::endl;
        msg_Debugging() << "Manual fraction beta*alpha: " << alpha_fraction  << std::endl;
        wmax_manual_list[curr_epsilon_max_manual_index] = fin_w_max;
        alpha_manual_list[curr_epsilon_max_manual_index] = alpha;
        alpha_manual_fraction_list[curr_epsilon_max_manual_index] = alpha_fraction;
        efficiency_manual_list[curr_epsilon_max_manual_index] = mean_efficiency*whisto_fills/p_whisto->Fills();//directly multiply with cut efficiency for manual definition
        if (curr_epsilon_max_manual_index==epsilon_max_values.size()-1) break;
        curr_epsilon_max_manual_index += 1;
        next_pxs_manual = whisto_sum*exp(log(10.)*epsilon_max_values[curr_epsilon_max_manual_index]);
      }
    }
    rpa->gen.SetXsecMap(p_proc->ResultsName(), whisto_sum/p_whisto->Fills()*m_enhancefac);//to be consistent with above psel: whisto_sum/p_whisto->Fills()*m_enhancefac
    rpa->gen.SetEfficiencyManualMap(p_proc->ResultsName(), efficiency_manual_list);
    rpa->gen.SetAlphaManualMap(p_proc->ResultsName(), alpha_manual_list);
    rpa->gen.SetAlphaManualFractionMap(p_proc->ResultsName(), alpha_manual_fraction_list);
    rpa->gen.SetWmaxManualMap(p_proc->ResultsName(), wmax_manual_list);


    cutxs_abs = 0.;
    cutxs = 0.;
    cutxs2 = 0.;
    cnt = 0.;
  }

  double pxs = whisto_sum*epsilon; //use this definition, because very slight differences, but most impotrant: mean enhancement missing!
  //double pxs = dabs(TotalResult())*p_whisto->Fills()*epsilon; //use this definition, because very slight differences, but most impotrant: mean enhancement missing!
  for (int i=first_filled_bin+1;i>last_filled_bin-1;i--) {
    double prev_cutxs2 = cutxs2;
    double prev_cutxs = cutxs;
    double prev_cutxs_abs = cutxs_abs;
    double prev_cnt = cnt;
    //bin middle times count is added
    double w = exp(log(10.)*(p_whisto->Xmin()+(i-0.5)*p_whisto->BinSize()));
    cutxs_abs+= p_whisto->Value(i) * w;
    cutxs+= p_whisto_pos->Value(i) * w;
    cutxs-= p_whisto_neg->Value(i) * w;
    cutxs2+= p_whisto->Value(i) * pow(w,2);
    cnt+= p_whisto->Value(i);
    if (cutxs>=pxs) {
      double fin_w_max = Min(exp(log(10.)*(p_whisto->Xmin()+i*p_whisto->BinSize())),dabs(m_max));
      if (cutxs==pxs) {//needed for epsilon=1
	fin_w_max = Min(exp(log(10.)*(p_whisto->Xmin()+(i-1)*p_whisto->BinSize())),dabs(m_max));
      } else {
	cutxs2 = prev_cutxs2;
	cutxs = prev_cutxs;
	cutxs_abs = prev_cutxs_abs;
	cnt = prev_cnt;
      }
      double mean_efficiency = (1*cnt+(whisto_abs_sum-cutxs_abs)/fin_w_max)/whisto_fills;
      //according to equation (39) in arXiv:2506.06203v2
      double effevperev = pow(whisto_sum/fin_w_max,2)/(((whisto_abs_sum-cutxs_abs+cutxs2/fin_w_max)/fin_w_max)*mean_efficiency*whisto_fills);
      rpa->gen.SetChosenEfficiencyMap(p_proc->ResultsName(), mean_efficiency);
      rpa->gen.SetChosenAlphaMap(p_proc->ResultsName(), effevperev);
      //without enhancement function: whisto_sum/p_whisto->Fills() = dabs(TotalResult())
      //with enhancement function:    whisto_sum/p_whisto->Fills() = dabs(TotalResult())*m_meanenhfunc
      //defined such that: dabs(TotalResult())*m_meanenhfunc = mean_w*cut_effi = whisto_sum/whisto_fills * (whisto_fills/p_whisto->Fills())
      m_meanenhfunc = whisto_sum/p_whisto->Fills()/dabs(TotalResult());
      m_effi = mean_efficiency;
      m_effevperev = effevperev;
      return fin_w_max;
    }
  }
  return m_max;
}

void Process_Integrator::SetUpEnhance(const int omode) 
{
  if (m_maxeps!=0.0 && !p_proc->IsGroup()) {
    double max(GetMaxEps(m_maxeps));
    if (omode || msg->LevelIsTracking()) {
      msg_Info()<<"  reduce max for "<<p_proc->ResultsName()<<" to "
		<<max/Max()<<" ( eps = "<<m_maxeps<<" -> exp. eff "
                <<  m_effi <<" ) "<<std::endl;
    }
    SetMax(max);
  }
  if (p_proc->IsGroup()) {
    double oldmax(m_max);
    m_max=0.0;
    m_effi=0.0;
    int wmode = ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE"));
    for (size_t i(0);i<p_proc->Size();++i) {
      (*p_proc)[i]->Integrator()->SetUpEnhance(msg_LevelIsTracking());
      m_max+=(*p_proc)[i]->Integrator()->Max();
      //need to sum up weighted with wsel: m_effi
      m_effi+=(*p_proc)[i]->Integrator()->Efficiency()*(*p_proc)[i]->Integrator()->SelectionWeight(wmode);
    }
    if (omode || p_proc->Parent()==p_proc)
      if (p_whisto_pos)
    msg_Info()<<"  reduce max for "<<p_proc->ResultsName()<<" to "
	      <<m_max/oldmax<<" ( eps = "<<m_maxeps<<" -> exp. eff "
	      << m_effi/SelectionWeight(wmode)<<" ) "<<std::endl;
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
  m_mssum    += value/enhance;
  m_mssumsqr += sqr(value/enhance);
#else
  m_sn++;
  m_ssum    += value/enhance;
  m_ssumsqr += sqr(value/enhance);
#endif
  double max=dabs(value)/dabs(p_proc->Last())*
    ATOOLS::Max(p_proc->LastPlus(),-p_proc->LastMinus());
  if (max>m_max)  m_max  = max;
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
  m_n=0;
  m_totalxs=m_totalsum=m_totalsumsqr=m_totalerr=0.0;
  m_smax=m_max=m_wmin=m_ssigma2=m_ssumsqr=m_ssum=0.0;
  m_sn=0;
  m_son=1;
  m_vsmax.clear(); 
  m_vsn.clear();   
  m_vsum.clear(); 
  if (p_proc->IsGroup() && mode==1)
    for (size_t i(0);i<p_proc->Size();++i) 
      (*p_proc)[i]->Integrator()->Reset(mode);
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
  m_ssum+=m_mssum;
  m_ssumsqr+=m_mssumsqr;
  m_msn=m_mssum=m_mssumsqr=0.0;
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

