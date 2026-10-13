#include "PHASIC++/Main/Timing_Statistics.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace {

  // names of Timing_Statistics::Type in the keys of the report
  const char *s_typenames[] = {"trial", "total", "PS", "ME", "gen", "overw", "maxoverw",
                               "kept", "keptcopy", "copy", "overhead_after",
                               "overhead_after_kept", "rerun", "pass"};

  double GetEffiUnc(long int n_num, long int n_denom)
  {
    //from gaussion error propagation with n_denom=n_num+rest
    return sqrt(double(n_num)*(n_denom-n_num)/pow(n_denom,3));
  }

  double GetUnc(double w, double w2, long int n)
  {
    return sqrt((w2-pow(w,2)/((double)n))/((double)n*((double)n-1.0)));
  }

  void printEffi(std::string sub_name, std::string name, std::string num_type, std::string denom_type, std::map<std::string, long int> &number_map)
  {
    long int n_num = number_map["n_"+num_type+"_"+sub_name];
    long int n_denom = number_map["n_"+denom_type+"_"+sub_name];
    double unc = GetEffiUnc(n_num, n_denom);
    msg_Info()<< sub_name << " : "<< name << " efficiency: "<<n_num/(double)n_denom<<" +-"<< unc <<std::endl;
  }

  void printTime(std::string sub_name, std::string time_type, std::map<std::string, double> &time_map, std::map<std::string, long int> &number_map)
  {
    double w = time_map["sum_"+time_type+"_"+sub_name];
    double w2 = time_map["sum2_"+time_type+"_"+sub_name];
    long int n = number_map["n_"+time_type+"_"+sub_name];
    double unc = GetUnc(w, w2, n);
    msg_Info()<< sub_name << " : "<<"Avg time per "<< time_type << " weight: "<<w/(double)n*1000<<" ms +-"<< unc*1000. <<" ms"<<std::endl;
  }

}

Timing_Statistics &Timing_Statistics::Instance()
{
  static Timing_Statistics timing;
  return timing;
}

Timing_Statistics::Timing_Statistics()
{
  Settings& s = Settings::GetMainSettings();
  m_level=s["TIMING_STATISTICS"].SetDefault(0).Get<int>();
  m_large_weight_fraction=s["TIMING_STATISTICS_LARGE_WEIGHT_FRACTION"].SetDefault(0.001).Get<double>();
}

Timing_Statistics::Stopwatch::Stopwatch()
{
  if (Instance().Active()) m_last=Now();
}

void Timing_Statistics::Stopwatch::Lap(const Type type)
{
  Timing_Statistics &timing(Instance());
  if (!timing.Active()) return;
  const Clock::time_point now(Now());
  Book(timing.Current(),type,now-m_last);
  m_last=now;
}

Timing_Statistics::Subprocess &Timing_Statistics::Current()
{
  if (!p_current) p_current=&m_subprocesses[""];
  return *p_current;
}

void Timing_Statistics::Book(Subprocess &sub, const Type type,
                             const Clock::duration &time)
{
  const double t(std::chrono::duration_cast<std::chrono::nanoseconds>(time).count()/1000000000.);
  ++sub.n[size_t(type)];
  sub.sum[size_t(type)]+=t;
  sub.sum2[size_t(type)]+=t*t;
}

void Timing_Statistics::StartEventGeneration()
{
  // the timing statistics describe the event generation only, the points of
  // the integration are not recorded
  m_active=Enabled();
}

void Timing_Statistics::StartEvent()
{
  m_pending=false;
}

void Timing_Statistics::Kept(const bool repetition)
{
  Subprocess &sub(Current());
  // kept events per subprocess, needed for the Neff/evt of every run, and of
  // those the repetitions of a hard event for its overweight
  ++sub.n[size_t(Type::kept)];
  if (repetition) ++sub.n[size_t(Type::keptcopy)];
  if (!m_active) return;
  const Clock::time_point now(Now());
  Book(sub,Type::overhead_after_kept,now-m_after);
  // The next repetition of this hard event for its overweight starts here,
  // its time after acceptance must not contain that of this event.
  m_after=now;
}

void Timing_Statistics::StartTrial(const Process_Integrator *pi)
{
  if (m_active && m_pending) {
    // the last accepted event or repetition was vetoed
    m_pending=false;
    Book(Current(),Type::overhead_after,Now()-m_after);
  }
  auto slot(m_slots.find(pi));
  if (slot==m_slots.end())
    slot=m_slots.emplace(pi,&m_subprocesses[pi->Process()->Name()]).first;
  p_current=slot->second;
  if (m_active) ++p_current->n[size_t(Type::trial)];
}

void Timing_Statistics::RejectTrial(const Stopwatch &trial, const bool pass)
{
  if (!m_active) return;
  const Clock::duration time(Now()-trial.m_last);
  Book(Current(),Type::total,time);
  if (pass) Book(Current(),Type::pass,time);
}

void Timing_Statistics::AcceptTrial(const Stopwatch &trial, const bool overweight)
{
  Subprocess &sub(Current());
  // accepted events per subprocess, needed for the Neff/evt of every run
  ++sub.n[size_t(Type::gen)];
  if (!m_active) return;
  // End of the trial. From here on, including the re-run of the pilot run, the
  // time is the time after acceptance, until the event is kept or vetoed.
  m_after=Now();
  m_pending=true;
  if (overweight) ++sub.n[size_t(Type::overw)];
  Book(sub,Type::total,m_after-trial.m_last);
  Book(sub,Type::pass,m_after-trial.m_last);
}

void Timing_Statistics::Count(const Type type)
{
  if (m_active) ++Current().n[size_t(type)];
}

void Timing_Statistics::StartRerun()
{
  m_rerun_active=m_active;
  if (m_active) m_rerun=Now();
  m_active=false;
}

void Timing_Statistics::EndRerun()
{
  m_active=m_rerun_active;
  if (m_active) Book(Current(),Type::rerun,Now()-m_rerun);
}

void Timing_Statistics::StartRepetition()
{
  ++Current().n[size_t(Type::copy)];
  // its time after acceptance is booked like that of an accepted event, as
  // vetoed if it is vetoed
  if (m_active) m_pending=true;
}

void Timing_Statistics::SetCrossSection(const std::string &name, const double xs)
{
  m_xsec_map[name]=xs;
}

void Timing_Statistics::SetChosen(const std::string &name, const double efficiency,
                                  const double alpha)
{
  m_chosen_efficiency_map[name]=efficiency;
  m_chosen_alpha_map[name]=alpha;
}

void Timing_Statistics::SetCappedFraction(const std::string &name,
                                          const Weight_Distribution &wd,
                                          const int i, const double max)
{
  //save lost xsec fraction due to overweight_threshold
  double capped_fraction = 0;
  for (int k=i+1;k<wd.Last()+1;k++) {
    double w = wd.BinWeight(k);
    if (w>max*wd.OverweightThreshold()) {
      capped_fraction += wd.All().Value(k) * (w-max*wd.OverweightThreshold());
    }
  }
  m_chosen_capped_fraction_map[name] = capped_fraction/wd.AbsSum();
}

void Timing_Statistics::ScanMaxEpsilon(const Process_Integrator &pi,
                                       const Weight_Distribution &wd)
{
  const std::string &name(pi.Process()->ResultsName());
  const Histogram *p_whisto(&wd.All());
  const Histogram *p_whisto_pos(&wd.Pos()), *p_whisto_neg(&wd.Neg());
  const double whisto_abs_sum(wd.AbsSum()), whisto_sum(wd.Sum()), whisto_sum2(wd.Sum2());
  const int whisto_fills(wd.Fills());
  const std::vector<double> &epsilon_max_values(m_epsilon_values);
  const size_t neps(epsilon_max_values.size());
  // Scan of the epsilon values for the large-weight region, i.e. the fraction
  // target_fraction of the largest weights of |w|. For every epsilon: the
  // reduced maximum, the statistical dilution alpha of the whole sample and
  // alpha times beta of the large-weight region, the efficiency including
  // the cut efficiency, and the generated events per accepted one for
  // unweighted events. An entry without a value stays -1.
  struct Eps_Scan { std::vector<double> wmax, alpha, alpha_beta_fraction, efficiency, genperacc; };
  auto scan = [&](const double target_fraction) {
    Eps_Scan res{std::vector<double>(neps,-1), std::vector<double>(neps,-1),
                 std::vector<double>(neps,-1), std::vector<double>(neps,-1),
                 std::vector<double>(neps,-1)};
    // sums over the weights below the large-weight region
    double whisto_sum2_fraction = 0;
    double whisto_sum_fraction = 0;
    double whisto_sum_abs_fraction = 0;
    int whisto_fills_fraction = 0;
    double fraction_pxs = (1-target_fraction)*whisto_abs_sum;
    // Only filled bins contribute, so start at the first filled bin, which can be
    // the underflow bin i=0 (one below would read bin -1).
    for (int i=wd.First();i<wd.Last()+1;i++) {
      if (p_whisto->Value(i)!=0) {
        double w = wd.BinWeight(i);
        //it is rounded up the target_fraction -> for fraction of 0 largest weight is picked
        if ((whisto_sum_abs_fraction+(p_whisto_pos->Value(i)+p_whisto_neg->Value(i))*w)>fraction_pxs) break;
        whisto_sum_abs_fraction += p_whisto->Value(i)*w;
        whisto_sum_fraction += p_whisto_pos->Value(i)*w;
        whisto_sum_fraction -= p_whisto_neg->Value(i)*w;
        whisto_sum2_fraction += p_whisto->Value(i)*w*w;
        whisto_fills_fraction += p_whisto->Value(i);
      }
    }

    //sums of |w|, w^2 and number of the points below the current bin edge
    double this_cutxs_abs = 0.;
    double this_cutxs2 = 0.;
    double this_cnt = 0.;
    int curr_epsilon_max_manual_index = neps-1;
    double next_pxs_manual = whisto_abs_sum*(1-exp(log(10.)*epsilon_max_values[curr_epsilon_max_manual_index]));
    for (int i=wd.First()-1;i<wd.Last()+1;i++) {
      // i=-1 if the underflow bin is filled
      if (i>=0) {
        double w = wd.BinWeight(i);
        this_cutxs_abs+= p_whisto->Value(i) * w;
        this_cutxs2+= p_whisto->Value(i) * pow(w,2);
        this_cnt+= p_whisto->Value(i);
      }
      //manual epsilon_max definition: all epsilon values whose threshold is
      //crossed in this bin, as for the chosen one in Process_Integrator::GetMaxEps()
      while (curr_epsilon_max_manual_index>=0 and this_cutxs_abs>=next_pxs_manual) {
        // reduced maximum, efficiency and statistical dilution alpha (eq. (39) in
        // arXiv:2506.06203v2) as for the chosen epsilon
        const Weight_Distribution::Unweighting u(wd.Unweight(i,this_cutxs_abs,this_cutxs2,this_cnt));
        const double fin_w_max(u.max), mean_efficiency(u.efficiency);
        msg_Debugging() << "Manual: " << epsilon_max_values[curr_epsilon_max_manual_index] << ": "<< mean_efficiency << std::endl;
        double mean_efficiency_fraction = (whisto_fills-this_cnt+(this_cutxs_abs-whisto_sum_abs_fraction)/fin_w_max)/(whisto_fills-whisto_fills_fraction);
        if (this_cutxs_abs<whisto_sum_abs_fraction) {//when over fraction, then efficiency is 1
          mean_efficiency_fraction = 1;
        }
        //in one variable, because "whisto_sum-whisto_sum_fraction" can be 0, but cancels in product
        double alpha_times_beta_fraction = 0;
        if (this_cutxs_abs>whisto_sum_abs_fraction) {//alpha is const else, because wmax is smaller than lowest w of fraction -> wmax cancels -> const
          alpha_times_beta_fraction = (whisto_sum*whisto_sum/whisto_fills*(whisto_abs_sum-whisto_sum_abs_fraction)/mean_efficiency)/(whisto_abs_sum)*(1.0/fin_w_max)/(((-whisto_sum_abs_fraction+this_cutxs_abs+(whisto_sum2-this_cutxs2)/fin_w_max)));
        } else {
          alpha_times_beta_fraction = (whisto_sum*whisto_sum/whisto_fills*(whisto_abs_sum-whisto_sum_abs_fraction)/mean_efficiency)/(whisto_abs_sum/mean_efficiency_fraction)*1.0/((whisto_sum2-whisto_sum2_fraction));
        }
        res.wmax[curr_epsilon_max_manual_index] = fin_w_max;
        res.alpha[curr_epsilon_max_manual_index] = u.alpha;
        res.alpha_beta_fraction[curr_epsilon_max_manual_index] = alpha_times_beta_fraction;
        res.efficiency[curr_epsilon_max_manual_index] = mean_efficiency*whisto_fills/p_whisto->Fills();//directly multiply with cut efficiency for manual definition
        res.genperacc[curr_epsilon_max_manual_index] = wd.GeneratedPerAccepted(fin_w_max,i,this_cutxs_abs,mean_efficiency*whisto_fills);
        if (--curr_epsilon_max_manual_index>=0)
          next_pxs_manual = whisto_abs_sum*(1-exp(log(10.)*epsilon_max_values[curr_epsilon_max_manual_index]));
      }
    }
    return res;
  };

  // An empty weight histogram has returned in GetMaxEps(), so every entry is set
  // below or stays -1 (no value for this epsilon).
  const Eps_Scan eps_scan(scan(m_large_weight_fraction));
  vector<double> wmax_manual_list(eps_scan.wmax);
  vector<double> alpha_manual_list(eps_scan.alpha);
  vector<double> alpha_manual_fraction_list(eps_scan.alpha_beta_fraction);
  vector<double> efficiency_manual_list(eps_scan.efficiency);
  // full unweighting and weighted events: no repetitions
  vector<double> genperacc_manual_list(eps_scan.genperacc);
  genperacc_manual_list.resize(neps+2,1.);
  wmax_manual_list.resize(neps+2,-1);
  alpha_manual_list.resize(neps+2,-1);
  efficiency_manual_list.resize(neps+2,-1);
  //unweighted with eps=0.0
  std::vector<double> unweighted_effiandeffevperev = pi.TotalEffiAndEffEvPerEv(1);
  efficiency_manual_list[neps]=unweighted_effiandeffevperev[0]*whisto_fills/p_whisto->Fills();//directly multiply with cut efficiency for manual definition
  alpha_manual_list[neps]=unweighted_effiandeffevperev[1];
  wmax_manual_list[neps]=wd.WeightMax();
  //weighted: every point that passes the cuts is an event, so that the efficiency
  //is the cut efficiency, as in the entries above, and alpha that of the weights
  efficiency_manual_list[neps+1]=double(whisto_fills)/p_whisto->Fills();
  std::vector<double> weighted_effiandeffevperev = pi.TotalEffiAndEffEvPerEv(0);
  alpha_manual_list[neps+1]=weighted_effiandeffevperev[1];
  wmax_manual_list[neps+1]=1;
  //save number of generatead events after cuts
  m_fills_map[name] = whisto_fills;//only to check average fill of weight-histo for Warning-printing
  m_efficiency_manual_map[name] = efficiency_manual_list;
  m_alpha_manual_map[name] = alpha_manual_list;
  m_alpha_manual_fraction_map[name] = alpha_manual_fraction_list;
  m_wmax_manual_map[name] = wmax_manual_list;
  m_genperacc_manual_map[name] = genperacc_manual_list;

  //the same scan for all fractions
  const std::vector<double> &fraction_values(m_epsilon_values);
  vector<vector<double>> alpha_manual_fscan_list(fraction_values.size());
  vector<vector<double>> efficiency_manual_fscan_list(fraction_values.size());
  for (size_t fi=0;fi<fraction_values.size();fi++) {
    const Eps_Scan fraction_scan(scan(exp(log(10.)*fraction_values[fi])));
    alpha_manual_fscan_list[fi] = fraction_scan.alpha_beta_fraction;
    efficiency_manual_fscan_list[fi] = fraction_scan.efficiency;
  }
  m_efficiency_manual_fscan_map[name] = efficiency_manual_fscan_list;
  m_alpha_manual_fscan_map[name] = alpha_manual_fscan_list;
}

std::map<std::string, double> Timing_Statistics::TimeMapAll() const
{
  std::map<std::string, double> time_map;
  for (const auto &[name, sub] : m_subprocesses)
    for (size_t t(0);t<size_t(Type::size);++t) {
      time_map[std::string("sum_")+s_typenames[t]+"_"+name]=sub.sum[t];
      time_map[std::string("sum2_")+s_typenames[t]+"_"+name]=sub.sum2[t];
    }
  return time_map;
}

std::map<std::string, long int> Timing_Statistics::NumberMapAll() const
{
  std::map<std::string, long int> number_map;
  for (const auto &[name, sub] : m_subprocesses)
    for (size_t t(0);t<size_t(Type::size);++t)
      number_map[std::string("n_")+s_typenames[t]+"_"+name]=sub.n[t];
  return number_map;
}

void Timing_Statistics::Report(const std::function<void()> &finish)
{
  //calculate sum_map etc.
  std::map<std::string, double>  time_map = TimeMapAll();
  std::map<std::string, long int> number_map = NumberMapAll();
  std::map<std::string, double>  xsec_map = m_xsec_map;
  // All subprocesses, including those without a trial in the event generation,
  // which have no entries in time_map and number_map.
  std::vector<std::string> subprocesses;
  for (auto const& [key, val] : xsec_map) subprocesses.push_back(key);
  std::vector<std::string> num_types = {"total", "trial", "PS", "ME", "gen", "overw", "maxoverw", "kept"};//
  std::map<std::string, double>  sum_map;
  std::map<std::string, double>  sum_mult_map;//defined as double to simplify following divisions
  for (const std::string& num_type : num_types) {
    sum_map[num_type] = 0;
  }
  for (const std::string& sub_name : subprocesses) {
    for (const std::string& num_type : num_types) {
      sum_map[num_type] += number_map["n_"+num_type+"_"+sub_name];
      std::string mult=sub_name.substr(0,4);
      if (sum_mult_map.find(mult+"_"+num_type) == sum_mult_map.end()) {
        sum_mult_map[mult+"_"+num_type] = 0;
      }
      sum_mult_map[mult+"_"+num_type] += number_map["n_"+num_type+"_"+sub_name];
    }
  }
  
  //calculate sudakov_efficiency per process (needs sum_map etc.)
  double sepsum = 0;
  std::map<std::string, double> sudakov_efficiency;
  //uncorrelated between subprocesses
  std::map<std::string, double> sudakov_efficiency_up;
  std::map<std::string, double> sudakov_efficiency_down;
  //correlated between subprocesses
  std::map<std::string, double> sudakov_efficiency_up_corr;
  std::map<std::string, double> sudakov_efficiency_down_corr;
  double average_sudakov = sum_map["kept"]/sum_map["gen"];
  for (const std::string& sub_name : subprocesses) {
	if (number_map["n_kept_"+sub_name]>0) {//take sub specific
	  sudakov_efficiency[sub_name] = number_map["n_kept_"+sub_name]*1.0/number_map["n_gen_"+sub_name];
	  double d_sud = sqrt(number_map["n_kept_"+sub_name]*1.0*max(1L, number_map["n_gen_"+sub_name]-number_map["n_kept_"+sub_name])/pow(number_map["n_gen_"+sub_name],3));
	  sudakov_efficiency_up[sub_name] = sudakov_efficiency[sub_name]+d_sud;
	  sudakov_efficiency_down[sub_name] = sudakov_efficiency[sub_name]-d_sud;
	  sudakov_efficiency_up_corr[sub_name] = sudakov_efficiency[sub_name];
	  sudakov_efficiency_down_corr[sub_name] = sudakov_efficiency[sub_name];
	} else {
	  std::string mult=sub_name.substr(0,4);
	  if (sum_mult_map.find(mult+"_kept") != sum_mult_map.end() and sum_mult_map[mult+"_kept"]>0) {
	    //assume average of same multi: e.g. 2_6
	    sudakov_efficiency[sub_name] = sum_mult_map[mult+"_kept"]/sum_mult_map[mult+"_gen"];
	  } else {
	    //assume average of whole sample
	    sudakov_efficiency[sub_name] = sum_map["kept"]/sum_map["gen"];
	  }
	  sudakov_efficiency_up[sub_name] = sudakov_efficiency[sub_name];
	  sudakov_efficiency_down[sub_name] = sudakov_efficiency[sub_name];
	  sudakov_efficiency_up_corr[sub_name] = 1.0;
	  sudakov_efficiency_down_corr[sub_name] = 0.0;
	}
	//msg_Info() << sub_name << ": sudakov_efficiency[sub_name]: " << sudakov_efficiency[sub_name] << std::endl;
	double this_sepsum = time_map["sum_overhead_after_"+sub_name]+time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_total_"+sub_name];
	sepsum += this_sepsum;
  }
  // For unweighted events a hard event with an overweight is repeated with unit
  // weight, see Signal_Processes::Treat(), so that sudakov_efficiency counts the
  // kept events including these repetitions per accepted hard event. The cross
  // section and the statistical dilution alpha refer to the hard events, each
  // counted once with its overweight, i.e. to the kept events without the
  // repetitions. For other event types both are the same.
  std::map<std::string, double> sudakov_efficiency_hard;
  double sum_keptcopy = 0;
  for (const std::string& sub_name : subprocesses)
    sum_keptcopy += number_map["n_keptcopy_"+sub_name];
  const double hard_fraction_all = sum_map["kept"]>0 ? 1.-sum_keptcopy/sum_map["kept"] : 1.;
  for (const std::string& sub_name : subprocesses) {
    const double nkept(number_map["n_kept_"+sub_name]);
    const double hard_fraction = nkept>0 ? 1.-number_map["n_keptcopy_"+sub_name]/nkept : hard_fraction_all;
    sudakov_efficiency_hard[sub_name] = sudakov_efficiency[sub_name]*hard_fraction;
  }
  // The tables below scan Max_Epsilon, for which the number of repetitions per
  // accepted event differs from that of this run. They therefore take the kept
  // events and the time after acceptance per generated event (accepted event or
  // repetition), and the generated events per accepted one of the scanned
  // epsilon from the weight histogram, see Process_Integrator::GetMaxEps(). A
  // repetition is assumed to take as long (apart from the re-run of the pilot
  // run, see time_after_acceptance below) and to be vetoed as often as the
  // accepted event, whose hard process it shares, so that the kept events per
  // generated event also give the kept hard events per accepted one.
  std::map<std::string, double> gen_per_acc_run;
  std::map<std::string, double> sudakov_gen, sudakov_gen_up, sudakov_gen_down;
  std::map<std::string, double> sudakov_gen_up_corr, sudakov_gen_down_corr;
  for (const std::string& sub_name : subprocesses) {
    const double ngen(number_map["n_gen_"+sub_name]);
    const double c(ngen>0 ? (ngen+number_map["n_copy_"+sub_name])/ngen : 1.);
    gen_per_acc_run[sub_name] = c;
    sudakov_gen[sub_name] = sudakov_efficiency[sub_name]/c;
    sudakov_gen_up[sub_name] = sudakov_efficiency_up[sub_name]/c;
    sudakov_gen_down[sub_name] = sudakov_efficiency_down[sub_name]/c;
    sudakov_gen_up_corr[sub_name] = sudakov_efficiency_up_corr[sub_name]/c;
    sudakov_gen_down_corr[sub_name] = sudakov_efficiency_down_corr[sub_name]/c;
  }
  Settings& s = Settings::GetMainSettings();
  double timing_statistics_large_weight_fraction=s["TIMING_STATISTICS_LARGE_WEIGHT_FRACTION"].SetDefault(0.001).Get<double>();
  double timing_statistics_det_sim=s["TIMING_STATISTICS_DET_SIM_IN_S"].SetDefault(0.0).Get<double>();
  int timing_statistics=s["TIMING_STATISTICS"].SetDefault(0).Get<int>();
  // Power p of alpha in the selection weight, see PHASIC::Process_Integrator::SelectionAlphaExponent()
  const size_t swmode=s["SELECTION_WEIGHT_MODE"].SetDefault(0).Get<int>();
  auto alpha_power = [swmode](const double alpha) {
    return PHASIC::Process_Integrator::SelectionAlphaPower(alpha,swmode);
  };
  // Kish denominator term sigma_i^2/(alpha_i*N_i) of one subprocess with N_i
  // accepted events, such that the dilution of the combined sample is
  // (sum sigma)^2/(sum N_i * sum sigma_i^2/(alpha_i*N_i)) for any allocation N_i.
  auto kish_term = [](const double xsec, const double alpha, const double nacc) {
    return (alpha>0. && nacc>0.) ? xsec*xsec/(alpha*nacc) : 0.;
  };
  int m_min_tot_unc_per_day=s["TIMING_STATISTICS_MIN_UNC_PER_DAY"].SetDefault(0).Get<int>();//min tot unc per computing time
  //only for the prediction tables below: the sampling itself follows SELECTION_WEIGHT_MODE
  int generation_mode=ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE"));
  //calculate chosen effevperev (needs sudakov_efficiency)
  std::map<std::string, double> chosen_alpha_map = m_chosen_alpha_map;
  std::map<std::string, double> chosen_efficiency_map = m_chosen_efficiency_map;
  double sum_p_unw = 0;
  double sum_p_eff_sign = 0;
  double sum_kish = 0;
  for (auto const& [key, val] : chosen_alpha_map) {
    std::string sub_name = key;
    // A subprocess without accepted events, e.g. without any point that passed
    // the cuts, has efficiency and alpha 0 (weighted events: -1) and contributes
    // nothing to the sample; its values would give 0/0 below.
    if (!(chosen_efficiency_map[sub_name]>0.) || !(chosen_alpha_map[sub_name]>0.)) continue;
    // Accepted events of the subprocess. For (partially) unweighted events they
    // follow from its selection weight, eff*|sigma|/(eff*alpha^p). For weighted
    // events the selection weight is that of the weighted sample including the
    // points that fail the cuts, which are not generated, so the generated events
    // of this run are taken; a subprocess without any does not contribute.
    double nacc;
    if (generation_mode==0) {
      nacc = number_map["n_gen_"+sub_name];
      if (nacc==0.) continue;
    } else {
      double curr_xsec = dabs(xsec_map[sub_name])/chosen_efficiency_map[sub_name]/alpha_power(chosen_alpha_map[sub_name]);
      nacc = chosen_efficiency_map[sub_name]*curr_xsec;
    }
    // generated events including the repetitions, effective events of the
    // hard events
    sum_p_unw += nacc*sudakov_efficiency[sub_name];
    sum_p_eff_sign += xsec_map[sub_name]*sudakov_efficiency_hard[sub_name];
    sum_kish += kish_term(xsec_map[sub_name],chosen_alpha_map[sub_name],nacc)*sudakov_efficiency_hard[sub_name];
  }
  double chosen_effevperev = pow(sum_p_eff_sign,2)/(sum_p_unw*sum_kish);
  msg_Info()<<"with "<< chosen_effevperev << " Neff/evt           "<<std::endl;
  finish();

  if (not timing_statistics) return;
  double m_ovwth = s["OVERWEIGHT_THRESHOLD"].SetDefault(1e12).Get<double>();

  //loop over subprocesses
  std::vector<std::string> sum_types = {"total", "PS", "ME", "overhead", "overhead_after_kept", "overhead_after", "rerun", "pass"};//
  std::map<std::string, double> total;
  for (const std::string& sum_type : sum_types) {
    total[sum_type] = 0;
  }
  for (const std::string& sub_name : subprocesses) {
      // subprocesses without a trial in the event generation have no timing
      if (timing_statistics>4 && number_map["n_trial_"+sub_name]>0) {
	  msg_Info() << std::endl;
	  for (const std::string& num_type : num_types) {
	    msg_Info()<<sub_name<<" : "<< "n_"<<num_type<<": "<<number_map["n_"+num_type+"_"+sub_name]<<std::endl;
	  }
	  for (const std::string& sum_type : sum_types) {
	    if (sum_type=="overhead") continue;
	    printTime(sub_name, sum_type, time_map, number_map);
	  }
	  printEffi(sub_name, "Cut", "ME", "trial", number_map);
	  printEffi(sub_name, "Unweighting", "gen", "ME", number_map);
	  printEffi(sub_name, "Total", "gen", "trial", number_map);
	  double w = time_map["sum_PS_"+sub_name]+time_map["sum_ME_"+sub_name];
	  msg_Info()<< sub_name << " : "<<"Summed time: "<<w<<" s"<<std::endl;
      }
      total["PS"] += time_map["sum_PS_"+sub_name];
      total["ME"] += time_map["sum_ME_"+sub_name];
      total["overhead"] += time_map["sum_total_"+sub_name]-(time_map["sum_PS_"+sub_name]+time_map["sum_ME_"+sub_name]);

      double w = time_map["sum_total_"+sub_name];
      if (timing_statistics>4) {
        msg_Info()<< sub_name << " : "<<"Total time: "<<w<<" s"<<std::endl;
      }
      total["total"] += w;
      w = time_map["sum_overhead_after_kept_"+sub_name];
      if (timing_statistics>4) {
        msg_Info()<< sub_name << " : "<<"overhead_after_kept time: "<<w<<" s"<<std::endl;
      }
      total["overhead_after_kept"] += w;
      w = time_map["sum_overhead_after_"+sub_name];
      if (timing_statistics>4) {
        msg_Info()<< sub_name << " : "<<"overhead_after time: "<<w<<" s"<<std::endl;
      }
      total["overhead_after"] += w;
      // contained in overhead_after_kept and overhead_after
      total["rerun"] += time_map["sum_rerun_"+sub_name];
      // contained in total
      total["pass"] += time_map["sum_pass_"+sub_name];
  }
  if (timing_statistics>4) {
    for (const std::string& num_type : num_types) {
      msg_Info()<<"sum_map['n_"<<num_type<<"']: "<<sum_map[num_type]<<std::endl;
    }
    for ( const auto &myPair : sum_mult_map) {
      msg_Info()<<"sum_mult_map['"<<myPair.first<<"']: "<<sum_mult_map[myPair.first]<<std::endl;
    }
    for (const std::string& sum_type : sum_types) {
      msg_Info()<<"total['"<<sum_type<<"']: "<<total[sum_type]<<std::endl;
    }
  }
  //following for chosen epsilon max
  //Total time and splitting into unweighting steps
  total["overhead_during"] = total["overhead"];
  total["overhead_after"] = total["overhead_after_kept"]+total["overhead_after"];
  if (timing_statistics>4) {
    msg_Info()<<std::endl<<"Total time: "<< total["total"]+total["overhead_after"] <<" s"<<std::endl;
    msg_Info()<<" (It sums up to 100%. 'overhead' is split into 'during' and 'after' unweighting.)"<<std::endl;
    for (const std::string& sum_type : {"PS", "ME", "overhead_during", "overhead_after"}) {
      msg_Info()<<" "<< std::right<<std::setw(2) << round(total[sum_type]/(total["total"]+total["overhead_after"])*100)<<"% in '"<<sum_type<<"'"<<std::endl;
    }
  }
  if (timing_statistics>4) msg_Info() << "Total separate sum: " << sepsum << "s resulting in: " << 60*60*24/sepsum*sum_map["kept"] << "evts/day" << std::endl;

  //epsilon_max scan manual definition
  std::vector<double> epsilon_values = m_epsilon_values;
  std::vector<double> fraction_values = m_epsilon_values;
  std::map<std::string, std::vector<double>> alpha_manual_map = m_alpha_manual_map;
  std::map<std::string, std::vector<double>> wmax_manual_map = m_wmax_manual_map;
  std::map<std::string, std::vector<double>> efficiency_manual_map = m_efficiency_manual_map;
  // rows: the scanned Max_Epsilon values, full unweighting (size), weighted
  // (size+1) and partially unweighted with Max_Epsilon=1 (size+2)
  std::vector<double> mean_manual_alpha(epsilon_values.size()+3, -1);
  std::vector<double> mean_manual_events(epsilon_values.size()+3, -1);
  std::vector<double> mean_manual_eff_events(epsilon_values.size()+3, -1);
  std::vector<double> mean_manual_events_up(epsilon_values.size()+1, -1);
  std::vector<double> mean_manual_events_down(epsilon_values.size()+1, -1);
  std::vector<std::map<std::string, double>> mean_manual_events_time(epsilon_values.size()+1);
  int optimal_manual_i = 0;
  double optimal_manual_sum_t_trial = 0;
  // Sum of |xsec| of all subprocesses, for the shares in the tables below. Taken
  // from xsec_map and not accumulated in the scan, whose maps stay empty for
  // Max_Epsilon=0 (GetMaxEps() is not called then).
  double plain_xsec_sum = 0;
  for (auto const& [key, val] : xsec_map) plain_xsec_sum += dabs(val);
  // Time per trial of a subprocess and per generated event after acceptance
  // (shower etc.), measured in this run. A subprocess without trials or without
  // accepted events gets 0, which underestimates its time; see the warning on
  // untimed subprocesses below.
  auto time_per_trial = [&](const std::string &sub_name) {
    const long int ntot = number_map["n_total_"+sub_name];
    return ntot>0 ? time_map["sum_total_"+sub_name]/ntot : 0.;
  };
  // Time after acceptance per accepted event with gen generated events (the
  // accepted one and its repetitions), including the detector simulation of the
  // kept ones. The re-run of the pilot run (see
  // Matrix_Element_Handler::GenerateOneTrialEvent()) is contained in the time
  // after acceptance, but occurs once per accepted event and not for the
  // repetitions, so it is taken out of the time per generated event.
  auto time_after_acceptance = [&](const std::string &sub_name, const double gen) {
    const long int nacc = number_map["n_gen_"+sub_name];
    const long int ngen = nacc+number_map["n_copy_"+sub_name];
    if (ngen==0) return 0.;
    const double rerun = time_map["sum_rerun_"+sub_name];
    const double after = (time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_overhead_after_"+sub_name]-rerun)/ngen;
    return rerun/nacc+gen*(after+sudakov_gen[sub_name]*timing_statistics_det_sim);
  };
  // Generated events per accepted one for the scanned epsilon i, see above;
  // repetitions only for unweighted events.
  std::map<std::string, std::vector<double>> genperacc_manual_map = m_genperacc_manual_map;
  auto gen_per_acc = [&](const std::string &sub_name, const size_t i) {
    if (generation_mode!=1) return 1.;
    const std::vector<double> &c(genperacc_manual_map[sub_name]);
    return (i<c.size() && c[i]>0.) ? c[i] : 1.;
  };
  // time per trial including the time after acceptance of the eff accepted
  // events per trial with gen generated events per accepted one
  auto trial_time = [&](const std::string &sub_name, const double eff, const double gen) {
    return time_per_trial(sub_name)+eff*time_after_acceptance(sub_name,gen);
  };
  for(int i=0; i < epsilon_values.size()+1; i++){
    double sum_t_trial = 0;
    double sum_p_unw = 0;
    double sum_p_eff_sign = 0;
    double sum_kish = 0;
    double sum_p_unw_up = 0;
    double sum_p_unw_down = 0;
    double sum_p_unw_up_corr = 0;
    double sum_p_unw_down_corr = 0;
    double sum_t_ME = 0;
    double sum_t_PS = 0;
    double sum_t_ov_during = 0;
    double sum_t_ov_after = 0;
    //double current_epsilon = epsilon_values[i];
    for (auto const& [key, val] : alpha_manual_map) {
	std::string sub_name = key;
	//std::cout << sub_name << std::endl;
	//need to weight with sampling probability. Why not sudakov? - bacause happens afterwards - but still more events needed for optimal eff events? no
	// Skip a subprocess without a valid value for this epsilon (-1: not set), which
	// would turn the sums into NaN or, for SELECTION_WEIGHT_MODE 1 and 2, into
	// finite but wrong numbers.
	if (!(efficiency_manual_map[sub_name][i]>0.) || !(alpha_manual_map[sub_name][i]>0.)) {
	  msg_Info() << "WARNING: for " << sub_name << " there is no alpha value for i=" << i << " corresponding to eps=" << exp(log(10)*epsilon_values[i]) << ", skipped." << std::endl;
	  continue;
	}
	double curr_xsec = dabs(xsec_map[sub_name])/efficiency_manual_map[sub_name][i]/alpha_power(alpha_manual_map[sub_name][i]);
	// generated events per accepted one for this epsilon
	const double gen = gen_per_acc(sub_name,i);
	double overhead_after = time_after_acceptance(sub_name,gen);
	double t_trial = trial_time(sub_name,efficiency_manual_map[sub_name][i],gen);
	if (m_min_tot_unc_per_day) curr_xsec = t_trial?dabs(xsec_map[sub_name])/sqrt(t_trial*efficiency_manual_map[sub_name][i]*alpha_manual_map[sub_name][i]):0;///(sudakov_efficiency[sub_name]?sudakov_efficiency[sub_name]:average_sudakov)
	sum_t_trial += t_trial*curr_xsec;
	// kept events including the repetitions, effective events of the hard events
	sum_p_unw += efficiency_manual_map[sub_name][i]*gen*sudakov_gen[sub_name]*curr_xsec;
	sum_p_eff_sign += xsec_map[sub_name]*sudakov_gen[sub_name];
	sum_kish += kish_term(xsec_map[sub_name],alpha_manual_map[sub_name][i],efficiency_manual_map[sub_name][i]*curr_xsec)*sudakov_gen[sub_name];

	//100% correlated for very few kept events
	sum_p_unw_up_corr += efficiency_manual_map[sub_name][i]*gen*sudakov_gen_up_corr[sub_name]*curr_xsec;
	sum_p_unw_down_corr += efficiency_manual_map[sub_name][i]*gen*sudakov_gen_down_corr[sub_name]*curr_xsec;
	//0% correlated for a lot kept events - could even make 2 sudakov uncertainties for this ...
	sum_p_unw_up += pow(efficiency_manual_map[sub_name][i]*gen*(sudakov_gen_up[sub_name]-sudakov_gen[sub_name])*curr_xsec,2);
	sum_p_unw_down += pow(efficiency_manual_map[sub_name][i]*gen*(sudakov_gen[sub_name]-sudakov_gen_down[sub_name])*curr_xsec,2);

	//split by evegen step
	if (number_map["n_total_"+sub_name]>0) {
	  sum_t_ME += time_map["sum_ME_"+sub_name]/number_map["n_total_"+sub_name]*curr_xsec;
	  sum_t_PS += time_map["sum_PS_"+sub_name]/number_map["n_total_"+sub_name]*curr_xsec;
	  sum_t_ov_during += (time_map["sum_total_"+sub_name]-time_map["sum_ME_"+sub_name]-time_map["sum_PS_"+sub_name])/number_map["n_total_"+sub_name]*curr_xsec;
	}
	sum_t_ov_after += overhead_after*efficiency_manual_map[sub_name][i]*curr_xsec;
    }
    mean_manual_alpha[i] = pow(sum_p_eff_sign,2)/(sum_p_unw*sum_kish);
    mean_manual_events[i] = 60*60*24/(sum_t_trial/sum_p_unw);
    mean_manual_eff_events[i] = mean_manual_events[i]*mean_manual_alpha[i];
    mean_manual_events_up[i] = 60*60*24/(sum_t_trial/(sum_p_unw_up_corr+sqrt(sum_p_unw_up)));
    mean_manual_events_down[i] = 60*60*24/(sum_t_trial/(sum_p_unw_down_corr-sqrt(sum_p_unw_down)));
    mean_manual_events_time[i]["ME"] = sum_t_ME;
    mean_manual_events_time[i]["PS"] = sum_t_PS;
    mean_manual_events_time[i]["ov_during"] = sum_t_ov_during;
    mean_manual_events_time[i]["ov_after"] = sum_t_ov_after;
    mean_manual_events_time[i]["sum"] = sum_t_trial;
    //"=" to update optimal_manual_sum_eff, if i=0 is optimal
    if (mean_manual_eff_events[i]>=mean_manual_eff_events[optimal_manual_i] && epsilon_values.size()>i) {
	optimal_manual_i=i;
	optimal_manual_sum_t_trial=sum_t_trial;
    }
  }
  //for weighted: take care of potentially non-optimal selection weight
  // Every selected trial that passes the cuts is an event. With the cut efficiency
  // eps of the subprocess, a trial takes the time of a trial of this run plus eps
  // times the time after acceptance, as in the scans above, and yields eps events.
  // The selection weight per trial is |sigma|/alpha_all^p with alpha_all =
  // eps*alpha the dilution of all trial points, see
  // Process_Integrator::SelectionWeight(). Weighted events have no pilot run: a
  // trial that passes the cuts is calculated once, with the variations. With the
  // pilot run, this is the re-run, which the time after acceptance contains, so
  // that the pilot trial that passed the cuts (sum_pass) is taken out.
  // Partially unweighted events with Max_Epsilon=1 are the same sample, every
  // trial that passes the cuts is accepted, but with the pilot run: their time
  // per trial is that of this run (trial_time), which includes it.
  auto weighted_trial_time = [&](const std::string &sub_name, const double eps) {
    const long int npass = number_map["n_pass_"+sub_name];
    const double pilot = (number_map["n_rerun_"+sub_name]>0 && npass>0) ?
      time_map["sum_pass_"+sub_name]/npass : 0.;
    return time_per_trial(sub_name)+eps*(time_after_acceptance(sub_name,1.)-pilot);
  };
  int i = epsilon_values.size()+1;
  double sum_t_trial = 0;
  double sum_t_trial_pilot = 0;
  sum_p_unw = 0;
  double sum_xsec = 0;
  double sum_complex = 0;
  double sum_effiselw = 0;
  for (auto const& [key, val] : alpha_manual_map) {
    std::string sub_name = key;
    //std::cout << sub_name << std::endl;
    const double eps = efficiency_manual_map[sub_name][i];
    const double alpha = alpha_manual_map[sub_name][i];
    if (!(eps>0.) || !(alpha>0.)) continue;
    const double selw = dabs(xsec_map[sub_name])/alpha_power(eps*alpha);
    const double nacc = eps*selw;
    // weighted events are not repeated for an overweight
    sum_t_trial += weighted_trial_time(sub_name,eps)*selw;
    sum_t_trial_pilot += trial_time(sub_name,eps,1.)*selw;
    sum_p_unw += sudakov_gen[sub_name]*nacc;

    sum_xsec += xsec_map[sub_name]*sudakov_gen[sub_name];
    sum_effiselw += sudakov_gen[sub_name]*nacc;
    //std::cout << " alpha_manual_map[sub_name][i]:" << alpha_manual_map[sub_name][i] << std::endl;
    //std::cout << " sudakov_efficiency[sub_name]:" << sudakov_efficiency[sub_name] << std::endl;
    //std::cout << " selw:" << selw << std::endl;
    //std::cout << " xsec_map[sub_name]:" << xsec_map[sub_name] << std::endl;
    //std::cout << " sum_complex+:" << pow(xsec_map[sub_name]*sudakov_efficiency[sub_name],2)/(alpha_manual_map[sub_name][i]*sudakov_efficiency[sub_name]*selw) << std::endl;
    sum_complex += pow(xsec_map[sub_name],2)*sudakov_gen[sub_name]/(alpha*nacc);
    //std::cout << " sum_complex:" << sum_complex << std::endl;
  }
  mean_manual_alpha[i] = pow(sum_xsec,2)/(sum_effiselw*sum_complex);//von Zettel
  mean_manual_events[i] = 60*60*24/(sum_t_trial/sum_p_unw);
  mean_manual_eff_events[i] = mean_manual_events[i]*mean_manual_alpha[i];
  mean_manual_alpha[i+1] = mean_manual_alpha[i];
  mean_manual_events[i+1] = 60*60*24/(sum_t_trial_pilot/sum_p_unw);
  mean_manual_eff_events[i+1] = mean_manual_events[i+1]*mean_manual_alpha[i+1];

  //same as above, but for fraction 0.001
  //epsilon_max scan manual definition
  std::map<std::string, std::vector<double>> alpha_manual_fraction_map = m_alpha_manual_fraction_map;
  std::vector<double> mean_manual_fraction_alpha(epsilon_values.size()+1, -1);
  std::vector<double> mean_manual_fraction_events(epsilon_values.size()+1, -1);
  std::vector<double> mean_manual_fraction_eff_events(epsilon_values.size()+1, -1);
  int optimal_manual_fraction_i = 0;
  double optimal_manual_fraction_sum_t_trial = 0;
  for(int i=0; i < epsilon_values.size(); i++){
    double sum_complex = 0;
    double sum_t_trial = 0;
    double sum_p_unw = 0;
    double sum_p_eff_sign = 0;
    for (auto const& [key, val] : alpha_manual_fraction_map) {
	std::string sub_name = key;
	//msg_Info() << sub_name << std::endl;
	//need to weight with sampling probability. Why not sudakov? - bacause happens afterwards - but still more events needed for optimal eff events? no
	if (!(efficiency_manual_map[sub_name][i]>0.) || !(alpha_manual_map[sub_name][i]>0.) ||
	    !(alpha_manual_fraction_map[sub_name][i]>0.)) {
	  msg_Info() << "WARNING: for " << sub_name << " there is no alpha value for i=" << i << " corresponding to eps=" << exp(log(10)*epsilon_values[i]) << ", skipped." << std::endl;
	  continue;
	}
	double curr_xsec = dabs(xsec_map[sub_name])/efficiency_manual_map[sub_name][i]/alpha_power(alpha_manual_map[sub_name][i]);
	const double gen = gen_per_acc(sub_name,i);
	double t_trial = trial_time(sub_name,efficiency_manual_map[sub_name][i],gen);
	if (m_min_tot_unc_per_day) curr_xsec = t_trial?dabs(xsec_map[sub_name])/sqrt(t_trial*efficiency_manual_map[sub_name][i]*alpha_manual_map[sub_name][i]):0;///(sudakov_efficiency[sub_name]?sudakov_efficiency[sub_name]:average_sudakov)
      sum_t_trial += t_trial*curr_xsec;
      sum_p_unw += efficiency_manual_map[sub_name][i]*gen*sudakov_gen[sub_name]*curr_xsec;
      sum_p_eff_sign += xsec_map[sub_name]*sudakov_gen[sub_name];
	//sum_complex += pow(xsec_map[sub_name]*sudakov_efficiency[sub_name],2)/(alpha_manual_fraction_map[sub_name][i]*dabs(xsec_map[sub_name])/alpha_manual_map[sub_name][i]*sudakov_efficiency[sub_name]);
	sum_complex += pow(xsec_map[sub_name],2)*sudakov_gen[sub_name]/(alpha_manual_fraction_map[sub_name][i]*efficiency_manual_map[sub_name][i]*curr_xsec);
	//sum_complex += dabs(xsec_map[sub_name])*sudakov_efficiency[sub_name]/(alpha_manual_fraction_map[sub_name][i]/alpha_manual_map[sub_name][i]);
    }
    mean_manual_fraction_alpha[i] = pow(sum_p_eff_sign,2)/(sum_p_unw*sum_complex);
    mean_manual_fraction_events[i] = 60*60*24/(sum_t_trial/sum_p_unw);
    mean_manual_fraction_eff_events[i] = mean_manual_fraction_events[i]*mean_manual_fraction_alpha[i];
    //"=" to update optimal_manual_sum_eff, if i=0 is optimal
    if (mean_manual_fraction_eff_events[i]>=mean_manual_fraction_eff_events[optimal_manual_fraction_i]) {
	optimal_manual_fraction_i=i;
	optimal_manual_fraction_sum_t_trial=sum_t_trial;
    }
  }





  //same as above, but for all fractions
  //fraction scan manual definition
  std::map<std::string, std::vector<std::vector<double>>> alpha_manual_fscan_map = m_alpha_manual_fscan_map;
  std::map<std::string, std::vector<std::vector<double>>> efficiency_manual_fscan_map = m_efficiency_manual_fscan_map;
  std::vector<double> mean_manual_alpha_fscan(fraction_values.size(), -1);
  std::vector<double> mean_manual_events_fscan(fraction_values.size(), -1);
  std::vector<double> mean_manual_eff_events_fscan(fraction_values.size(), -1);
  for(int fi=0; fi < fraction_values.size(); fi++){
    double sum_t_trial = 0;
    double sum_p_unw = 0;
    double sum_p_eff_sign = 0;
    double sum_kish = 0;
    for (auto const& [key, val] : alpha_manual_fscan_map) {
	std::string sub_name = key;
	//msg_Info() << "  " << sub_name << std::endl;
	double opt_t_trial = 0;
	double opt_p_unw = 0;
	double opt_p_eff = -1;
	double opt_p_eff_sign = 0;
	double opt_kish = 0;
	for(int i=0; i < epsilon_values.size(); i++){
	  //msg_Info() << "   " << i << std::endl;
	  if (!(efficiency_manual_fscan_map[sub_name][fi][i]>0.) || !(alpha_manual_fscan_map[sub_name][fi][i]>0.)) {
	    msg_Info() << "WARNING: for " << sub_name << " there is no alpha value for fraction=" << exp(log(10)*fraction_values[fi]) << " and eps=" << exp(log(10)*epsilon_values[i]) << ", skipped." << std::endl;
	    continue;
	  }
	  double curr_xsec = dabs(xsec_map[sub_name])/efficiency_manual_fscan_map[sub_name][fi][i]/alpha_power(alpha_manual_fscan_map[sub_name][fi][i]);
	  // the maximum, and so the repetitions, depend on epsilon only
	  const double gen = gen_per_acc(sub_name,i);
	  double t_trial = trial_time(sub_name,efficiency_manual_fscan_map[sub_name][fi][i],gen);
	  if (m_min_tot_unc_per_day) curr_xsec = t_trial?dabs(xsec_map[sub_name])/sqrt(t_trial*efficiency_manual_fscan_map[sub_name][fi][i]*alpha_manual_fscan_map[sub_name][fi][i]):0;///(sudakov_efficiency[sub_name]?sudakov_efficiency[sub_name]:average_sudakov)
	  double this_t_trial = t_trial*curr_xsec;
	  double this_p_unw = efficiency_manual_fscan_map[sub_name][fi][i]*gen*sudakov_gen[sub_name]*curr_xsec;
	  double this_p_eff = dabs(xsec_map[sub_name])*sudakov_gen[sub_name];
	  double this_p_eff_sign = xsec_map[sub_name]*sudakov_gen[sub_name];
	  if (opt_p_eff==-1 or opt_t_trial/opt_p_eff/(this_t_trial/this_p_eff)>1.0001 or (opt_t_trial/opt_p_eff/(this_t_trial/this_p_eff)>0.9999 and opt_p_eff/opt_p_unw/(this_p_eff/this_p_unw)<1)) {
	    opt_t_trial = this_t_trial;
	    opt_p_unw = this_p_unw;
	    opt_p_eff = this_p_eff;
	    opt_p_eff_sign = this_p_eff_sign;
	    opt_kish = kish_term(xsec_map[sub_name],alpha_manual_fscan_map[sub_name][fi][i],efficiency_manual_fscan_map[sub_name][fi][i]*curr_xsec)*sudakov_gen[sub_name];
	  }
	}
	sum_t_trial += opt_t_trial;
	sum_p_unw += opt_p_unw;
	sum_p_eff_sign += opt_p_eff_sign;
	sum_kish += opt_kish;
    }
    mean_manual_alpha_fscan[fi] = pow(sum_p_eff_sign,2)/(sum_p_unw*sum_kish);
    mean_manual_events_fscan[fi] = 60*60*24/(sum_t_trial/sum_p_unw);
    mean_manual_eff_events_fscan[fi] = mean_manual_events_fscan[fi]*mean_manual_alpha_fscan[fi];
  }








  
  //make nice final table IV (for chosen emax)
  if (timing_statistics>3) {
    std::cout << "┌────────────────────────────────┬────────────────┬─────────────────────────────────────────────────────────┐" << std::endl;
    std::cout << "│     sampling contribution      │                │                                                         │" << std::endl;
    std::cout << "│ xsec*enh  efficiency  stat.dil │ time |xsec*enh|│ subprocess                                              │" << std::endl;
    std::cout << "├────────────────────────────────┼────────────────┼─────────────────────────────────────────────────────────┤" << std::endl;
    const auto default_precision{std::cout.precision()};
    std::cout << std::setprecision(3);
    for (const std::string& sub_name : subprocesses) {
      std::cout<<"│  "<< std::left<<std::setw(9) << xsec_map[sub_name]<<" ";
      std::cout<<" "<< std::left<<std::setw(9) << chosen_efficiency_map[sub_name]<<"";
      std::cout<<" "<< std::left<<std::setw(8) << chosen_alpha_map[sub_name]<<" ";
      double this_sepsum = time_map["sum_overhead_after_"+sub_name]+time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_total_"+sub_name];
      std::cout <<"│ " <<std::right<<std::setw(2) << round(this_sepsum/(total["overhead_after"]+total["total"])*100)<<"%  ";
      std::cout <<"  " <<std::right<<std::setw(4) << round(dabs(xsec_map[sub_name])/plain_xsec_sum*1000.)/10.<<"%   ";
      std::cout <<"│ " <<std::left<<std::setw(55)<< sub_name << " │" << std::endl;
    }
    std::cout << std::setprecision(default_precision);
    std::cout << "└────────────────────────────────┴────────────────┴─────────────────────────────────────────────────────────┘" << std::endl << std::endl;
  }

  //make nice final table III (for chosen emax)
  if (timing_statistics>2) {
    std::cout << "┌─────────────────────────┬────────────────┬─────────────────────────────────────────────────────────┐" << std::endl;
    std::cout << "│    time contribution    │                │                                                         │" << std::endl;
    std::cout << "│ PS   ME   ov.h.  shower │ time |xsec*enh|│ subprocess                                              │" << std::endl;
    std::cout << "├─────────────────────────┼────────────────┼─────────────────────────────────────────────────────────┤" << std::endl;
    for (const std::string& sub_name : subprocesses) {
      if (number_map["n_total_"+sub_name]>0) {
        double t_ME = time_map["sum_ME_"+sub_name]/number_map["n_total_"+sub_name];
        double t_PS = time_map["sum_PS_"+sub_name]/number_map["n_total_"+sub_name];
        double t_ov_during = (time_map["sum_total_"+sub_name]-time_map["sum_ME_"+sub_name]-time_map["sum_PS_"+sub_name]-time_map["sum_PDF_"+sub_name])/number_map["n_total_"+sub_name];
        double t_ov_after = (time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_overhead_after_"+sub_name])/number_map["n_total_"+sub_name];
        double sum = t_ME+t_PS+t_ov_during+t_ov_after;
        msg_Info()<<"│ "<< std::right<<std::setw(2) << round(t_PS/sum*100)<<"% ";
        msg_Info()<<" "<< std::right<<std::setw(2) << round(t_ME/sum*100)<<"% ";
        msg_Info()<<"  "<< std::right<<std::setw(2) << round(t_ov_during/sum*100)<<"% ";
        msg_Info()<<"   "<< std::right<<std::setw(2) << round(t_ov_after/sum*100)<<"%   ";
      } else {
        std::cout << "│ --   --    --     --   ";
      }
      double this_sepsum = time_map["sum_overhead_after_"+sub_name]+time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_total_"+sub_name];
      std::cout <<"│ " <<std::right<<std::setw(2) << round(this_sepsum/(total["overhead_after"]+total["total"])*100)<<"%  ";
      std::cout <<"  " <<std::right<<std::setw(4) << round(dabs(xsec_map[sub_name])/plain_xsec_sum*1000.)/10.<<"%   ";
      std::cout <<"│ " <<std::left<<std::setw(55)<< sub_name << " │" << std::endl;
    }
    std::cout << "└─────────────────────────┴────────────────┴─────────────────────────────────────────────────────────┘" << std::endl;
  }

  //make nice final table II
  if (timing_statistics>1) {
    std::cout << "┌─────────────┬─────────────────────────────────────────────┐" << std::endl;
    std::cout << "│ Max_Epsilon │ PS time   ME time   overhead    shower etc  │" << std::endl;
    std::cout << "├─────────────┼─────────────────────────────────────────────┤" << std::endl;
    for(int i=0; i < epsilon_values.size(); i++){
      std::cout << "│ 1e" <<std::left<<std::setw(9)<< epsilon_values[i] << " │ ";
      msg_Info()<<" "<< std::right<<std::setw(2) << round(mean_manual_events_time[i]["PS"]/mean_manual_events_time[i]["sum"]*100)<<"%      ";
      msg_Info()<<" "<< std::right<<std::setw(2) << round(mean_manual_events_time[i]["ME"]/mean_manual_events_time[i]["sum"]*100)<<"%      ";
      msg_Info()<<" "<< std::right<<std::setw(2) << round(mean_manual_events_time[i]["ov_during"]/mean_manual_events_time[i]["sum"]*100)<<"%       ";
      msg_Info()<<" "<< std::right<<std::setw(2) << round(mean_manual_events_time[i]["ov_after"]/mean_manual_events_time[i]["sum"]*100)<<"%         │" << std::endl;
    }
    std::cout << "└─────────────┴─────────────────────────────────────────────┘" << std::endl;
  }

  //make nice final table I.5
  std::cout << "┌─────────────────────────────┬─────────────────────────────────┐" << std::endl;
  std::cout << "│                             │      Corresponding region       │" << std::endl;
  std::cout << "│ Fraction       Events/day   │ Eff. events/day     Sample size │" << std::endl;
  std::cout << "├─────────────────────────────┼─────────────────────────────────┤" << std::endl;
  const auto default_precision{std::cout.precision()};
  std::cout << std::setprecision(4);
  for(int i=0; i < fraction_values.size(); i++){
    std::cout << "│ 1e" <<std::left<<std::setw(9)<< fraction_values[i] << "    " <<std::setw(12)<< mean_manual_events_fscan[i] << " │ ";
    std::cout <<std::setw(15)<< mean_manual_eff_events_fscan[i] << "     " <<std::left<<std::setw(11) << 1./mean_manual_alpha_fscan[i] << " │ " << std::endl;
  }
  std::cout << "└─────────────────────────────┴─────────────────────────────────┘" << std::endl;
  
  //make nice final table I
  std::cout << "┌─────────────────────────────┬─────────────────────────────────┬─────────────────────────────────┐" << std::endl;
  std::cout << "│                             │          Average region         │       Large weight region       │" << std::endl;
  std::cout << "│                             │           (fraction: 1)         │        (fraction: "<< timing_statistics_large_weight_fraction <<")        │" << std::endl;
  std::cout << "│ Max_Epsilon    events/day   │ eff. events/day     sample size │ eff. events/day     sample size │" << std::endl;
  std::cout << "├─────────────────────────────┼─────────────────────────────────┼─────────────────────────────────┤" << std::endl;
  std::cout << "│ 0.0 (unw.) " <<std::left<< "    " <<std::setw(12)<< mean_manual_events[epsilon_values.size()] << " │ ";
  std::cout <<std::setw(15)<< mean_manual_eff_events[epsilon_values.size()] << "     " <<std::left<<std::setw(11) << 1./mean_manual_alpha[epsilon_values.size()] << " │ ";
  std::cout << "-                   -          " << " │" << std::endl;
  //std::cout <<std::setw(15)<< mean_manual_fraction_eff_events[epsilon_values.size()] << "     " <<std::left<<std::setw(11) << 1./mean_manual_fraction_alpha[epsilon_values.size()] << " │" << std::endl;
  for(int i=0; i < epsilon_values.size(); i++){
    std::cout << "│ 1e" <<std::left<<std::setw(9)<< epsilon_values[i] << "    " <<std::setw(12)<< mean_manual_events[i] << " │ ";
    if (i==optimal_manual_i) {
      std::cout <<std::setw(15)<< mean_manual_eff_events[i] << "<--  " <<std::left<<std::setw(11) << 1./mean_manual_alpha[i] << " │ ";
    } else {
      std::cout <<std::setw(15)<< mean_manual_eff_events[i] << "     " <<std::left<<std::setw(11) << 1./mean_manual_alpha[i] << " │ ";
    }
    if (i==optimal_manual_fraction_i) {
      std::cout <<std::setw(15)<< mean_manual_fraction_eff_events[i] << "<--  " <<std::left<<std::setw(11) << 1./mean_manual_fraction_alpha[i] << " │" << std::endl;
    } else {
      std::cout <<std::setw(15)<< mean_manual_fraction_eff_events[i] << "     " <<std::left<<std::setw(11) << 1./mean_manual_fraction_alpha[i] << " │" << std::endl;
    }
  }
  // For partially unweighted events this is the row 1e0, for weighted ones (no
  // pilot run measured) the row weighted. For unweighted ones the row 1e0
  // contains the repetitions of overweight events.
  if (generation_mode==1) {
    std::cout << "│ 1.0 (p.unw.)" <<std::left<< "   " <<std::setw(12)<< mean_manual_events[epsilon_values.size()+2] << " │ ";
    std::cout <<std::setw(15)<< mean_manual_eff_events[epsilon_values.size()+2] << "     " <<std::left<<std::setw(11) << 1./mean_manual_alpha[epsilon_values.size()+2] << " │ ";
    std::cout << "-                   -          " << " │" << std::endl;
  }
  std::cout << "│ weighted   " <<std::left<< "    " <<std::setw(12)<< mean_manual_events[epsilon_values.size()+1] << " │ ";
  std::cout <<std::setw(15)<< mean_manual_eff_events[epsilon_values.size()+1] << "     " <<std::left<<std::setw(11) << 1./mean_manual_alpha[epsilon_values.size()+1] << " │ ";
  std::cout << "-                   -          " << " │" << std::endl;
  std::cout << "└─────────────────────────────┴─────────────────────────────────┴─────────────────────────────────┘" << std::endl;
  //show relative sudakov uncertainty for average region optimal point
  std::cout << "Relative Sudakov uncertainty for average region optimal point: +"<< (mean_manual_events_up[optimal_manual_i]/mean_manual_events[optimal_manual_i]-1)*100. <<"% -"<< (1-mean_manual_events_down[optimal_manual_i]/mean_manual_events[optimal_manual_i])*100. <<"%" << std::endl;
  std::cout << std::setprecision(default_precision);
  if (timing_statistics_det_sim!=0) {
    std::cout << "With assumed detector simulation time of "<< timing_statistics_det_sim<<" s per kept event." << std::endl;
  } else {
    std::cout << "Note: No detector simulation time was considered for the above table(s). Use TIMING_STATISTICS_DET_SIM_IN_S to set it." << std::endl;
  }

  //time saving potential: compare the allocation of SELECTION_WEIGHT_MODE with
  //the time-optimal one of TIMING_STATISTICS_MIN_UNC_PER_DAY, which selects
  //subprocess i with |sigma_i|/sqrt(t_i*eff_i*alpha_i) (t_i: time per trial,
  //eff_i: accepted events per trial) and thereby minimises the total
  //cross-section uncertainty per computing time. Both use the timing measured
  //in this run and the efficiencies of the epsilon scan (including the cut
  //efficiency), and assume that the time per trial of a subprocess does not
  //depend on how often it is selected.
  {
    // sud: kept events per accepted one, sud_hard: the same without the
    // repetitions for an overweight, for the epsilon of the row (see sudakov_gen)
    struct Sub_Point { double xsec, eff, alpha, t_trial, sud, sud_hard; };
    auto eff_events_per_day = [&](const std::vector<Sub_Point> &subs, const bool timeopt) {
      double sum_t(0.), sum_n(0.), sum_sign(0.), sum_k(0.);
      for (const Sub_Point &sp : subs) {
        const double curr = timeopt ?
          (sp.t_trial>0. ? dabs(sp.xsec)/sqrt(sp.t_trial*sp.eff*sp.alpha) : 0.) :
          dabs(sp.xsec)/sp.eff/alpha_power(sp.alpha);
        sum_t += sp.t_trial*curr;
        sum_n += sp.eff*sp.sud*curr;
        sum_sign += sp.xsec*sp.sud_hard;
        sum_k += kish_term(sp.xsec,sp.alpha,sp.eff*curr)*sp.sud_hard;
      }
      if (!(sum_t>0.) || !(sum_n>0.) || !(sum_k>0.)) return 0.;
      return 60*60*24*sum_n/sum_t*sum_sign*sum_sign/(sum_n*sum_k);
    };
    // subprocesses without trials or without generated events have no or an
    // incomplete time per trial, which biases the time-optimal allocation
    std::vector<std::string> untimed;
    for (auto const& [key, val] : alpha_manual_map)
      if (number_map["n_total_"+key]==0 || number_map["n_gen_"+key]==0)
        untimed.push_back(key);
    // rows: the scanned Max_Epsilon values and full unweighting (i=size)
    const size_t nrows(epsilon_values.size()+1);
    std::vector<double> neff_sw(nrows,0.), neff_time(nrows,0.);
    for (size_t i(0);i<nrows;++i) {
      std::vector<Sub_Point> subs;
      for (auto const& [key, val] : alpha_manual_map) {
        const double eff(efficiency_manual_map[key][i]), alpha(alpha_manual_map[key][i]);
        // no value for this epsilon (see the warnings above)
        if (!(eff>0.) || !(alpha>0.)) continue;
        const double gen(gen_per_acc(key,i));
        subs.push_back({xsec_map[key],eff,alpha,trial_time(key,eff,gen),gen*sudakov_gen[key],
                        sudakov_gen[key]});
      }
      neff_sw[i] = eff_events_per_day(subs,false);
      neff_time[i] = eff_events_per_day(subs,true);
    }
    const size_t best_sw = std::max_element(neff_sw.begin(),neff_sw.end())-neff_sw.begin();
    const size_t best_time = std::max_element(neff_time.begin(),neff_time.end())-neff_time.begin();
    auto row_name = [&](const size_t i) {
      return i<epsilon_values.size() ? "1e"+ToString(epsilon_values[i]) : std::string("0.0 (unw.)");
    };
    // columns: 24 + 17 + 8 characters after "│ Max_Epsilon │ "
    std::cout << std::setprecision(4) << std::left;
    std::cout << "┌─────────────┬──────────────────────────────────────────────────┐" << std::endl;
    std::cout << "│             │ " << std::setw(49) << "eff. events/day, selection of subprocesses" << "│" << std::endl;
    std::cout << "│ Max_Epsilon │ " << std::setw(24) << "SELECTION_WEIGHT_MODE "+ToString(swmode)
              << std::setw(17) << "time-optimal" << std::setw(8) << "ratio" << "│" << std::endl;
    std::cout << "├─────────────┼──────────────────────────────────────────────────┤" << std::endl;
    for (size_t i(0);i<nrows;++i) {
      std::cout << "│ " << std::setw(11) << row_name(i) << " │ ";
      std::cout << std::setw(19) << neff_sw[i] << (i==best_sw ? "<--  " : "     ");
      std::cout << std::setw(13) << neff_time[i] << (i==best_time ? "<-- " : "    ");
      std::cout << std::setw(8) << (neff_sw[i]>0. ? neff_time[i]/neff_sw[i] : 0.) << "│" << std::endl;
    }
    std::cout << "└─────────────┴──────────────────────────────────────────────────┘" << std::endl;
    if (neff_sw[best_sw]>0. && neff_time[best_time]>0.)
      std::cout << "Time saving potential of the time-optimal selection (both at their optimal"
                << " Max_Epsilon): " << (1.-neff_sw[best_sw]/neff_time[best_time])*100.
                << "% less computing time for the same number of effective events." << std::endl;
    std::cout << std::setprecision(default_precision) << std::right;
    if (!untimed.empty()) {
      std::cout << "Warning: " << untimed.size() << " subprocess(es) have no trials or no generated"
                << " events in this run, so their time per trial is missing or incomplete and the"
                << " time-optimal allocation above is biased. Generate more events." << std::endl;
      for (const std::string &name : untimed) msg_Tracking() << "  " << name << std::endl;
    }
  }

  // Average over the subprocesses with a weight histogram, i.e. those for which
  // GetMaxEps() ran; for Max_Epsilon=0 there are none and no statement is made.
  std::map<std::string, int> whistofill_map = m_fills_map;
  double average_whistofill = 0;
  double fills_xsecsum = 0;
  for (auto const& [sub_name, fills] : whistofill_map) {
    average_whistofill += dabs(xsec_map[sub_name])*fills;
    fills_xsecsum += dabs(xsec_map[sub_name]);
  }
  if (fills_xsecsum>0.) {
    average_whistofill = average_whistofill/fills_xsecsum;
    if (average_whistofill<pow(10,5)) {
      std::cout << "Warning: There are on average only "<< average_whistofill << " events generated in the integration phase per process. For a reliable prediction of the impact of Max_Epslion choises more than 1e6 are recommended. Please increase PSI: MAX_OPT during integration." << std::endl;
    } else if (average_whistofill<pow(10,6)) {
      std::cout << "Info: There are on average only "<< average_whistofill << " events generated in the integration phase per process. For a reliable prediction of the impact of Max_Epslion choises more than 1e6 are recommended. One could increase PSI: MAX_OPT during integration." << std::endl;
    }
  }

  // capped cross section relative to the total |xsec| of all subprocesses
  std::map<std::string, double> chosen_capped_fraction_map = m_chosen_capped_fraction_map;
  double capfraction = 0;
  for (auto const& [sub_name, fraction] : chosen_capped_fraction_map)
    capfraction += dabs(xsec_map[sub_name])*fraction;
  capfraction = plain_xsec_sum>0. ? capfraction/plain_xsec_sum*100 : 0.;
  if (capfraction>0.01) {
    std::cout << "Warning: OVERWEIGHT_THRESHOLD is set to: " << m_ovwth << " which results in a loss of total cross section due to max weight capping of: " << std::round(capfraction*100)/100 << "%. Please consider chosing a smaller Max_Epsilon or a larger OVERWEIGHT_THRESHOLD." << std::endl;
  }
  if (generation_mode==0) {
    //there is no pilot run for weighted events and thus the corresponding overhead can not be estimated in this case. Try to convey this message in user friendly language:
    std::cout << "Warning: TIMING_STATISTICS was run in 'weighted' mode. The estimate of 'shower etc.' is only correctly determined from 'unweighted' or '(partially) unweighted' events." << std::endl;
  }
}
