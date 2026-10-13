#ifndef PHASIC_Main_Timing_Statistics_H
#define PHASIC_Main_Timing_Statistics_H

#include <array>
#include <chrono>
#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace PHASIC {

  class Process_Integrator;
  class Weight_Distribution;

  // Statistics of the event generation for TIMING_STATISTICS: the time of the
  // trials (phase space, matrix element) and after the acceptance of an event
  // (re-run of the pilot run, shower etc.) per subprocess, the epsilon scans of
  // the weight histograms, and the tables at the end of the run, which predict
  // the events per day for other Max_Epsilon values. The accepted and kept events
  // and the unweighting statistics of the chosen Max_Epsilon are recorded in every
  // run, for the Neff/evt of the summary.
  //
  // Without TIMING_STATISTICS no clock is read: the timing calls only check a flag.
  class Timing_Statistics {
  public:
    typedef std::chrono::high_resolution_clock Clock;

    // what is counted and, for the intervals, timed per subprocess
    enum class Type { trial, total, PS, ME, gen, overw, maxoverw, kept, keptcopy,
                      copy, overhead_after, overhead_after_kept, rerun, pass, size };

    // Times consecutive intervals of a trial of the current subprocess: started on
    // construction, Lap() books the time since the start or the previous lap.
    class Stopwatch {
      Clock::time_point m_last;
      friend class Timing_Statistics;
    public:
      Stopwatch();
      void Lap(const Type type);
    };

    static Timing_Statistics &Instance();

    // TIMING_STATISTICS>0: the epsilon scans are made and the tables printed
    inline bool Enabled() const { return m_level>0; }
    // the timing is recorded, from the first event on and not during the re-run
    inline bool Active() const { return m_active; }

    // event generation, see SHERPA::Sherpa::GenerateOneEvent()
    void StartEventGeneration();
    void StartEvent();
    // the event is kept (not vetoed): end of its time after acceptance
    void Kept(const bool repetition);

    // trials, see Phase_Space_Handler::Differential() and
    // SHERPA::Matrix_Element_Handler::GenerateOneTrialEvent(). StartTrial() closes
    // the time after acceptance of a vetoed event and makes pi the current
    // subprocess. pass: the trial gave a weight, i.e. passed the cuts.
    void StartTrial(const Process_Integrator *pi);
    void RejectTrial(const Stopwatch &trial, const bool pass);
    void AcceptTrial(const Stopwatch &trial, const bool overweight);
    void Count(const Type type);
    // re-run of an accepted event with all variations (pilot run): part of the
    // time after acceptance and timed on its own, its trial is not recorded
    void StartRerun();
    void EndRerun();
    // repetition of the last hard event for its overweight (unweighted events),
    // see SHERPA::Signal_Processes::Treat()
    void StartRepetition();

    // results of Process_Integrator::GetMaxEps()
    void SetCrossSection(const std::string &name, const double xs);
    void SetChosen(const std::string &name, const double efficiency, const double alpha);
    void ScanMaxEpsilon(const Process_Integrator &pi, const Weight_Distribution &wd);
    void SetCappedFraction(const std::string &name, const Weight_Distribution &wd,
                           const int i, const double max);

    // Summary line Neff/evt and, for TIMING_STATISTICS>0, the tables. finish() is
    // called in between (SHERPA::Event_Handler::Finish() prints the cross sections).
    void Report(const std::function<void()> &finish);

  private:

    struct Subprocess {
      std::array<long int,size_t(Type::size)> n{};
      std::array<double,size_t(Type::size)> sum{}, sum2{};
    };

    int m_level;
    double m_large_weight_fraction;
    bool m_active = false, m_rerun_active = false;
    // the time after acceptance of the last accepted event or repetition is
    // still to be booked as vetoed (overhead_after) or kept (overhead_after_kept)
    bool m_pending = false;
    Clock::time_point m_after, m_rerun;

    // per subprocess name, and a slot per integrator for the trials
    std::map<std::string, Subprocess> m_subprocesses;
    std::unordered_map<const Process_Integrator*, Subprocess*> m_slots;
    Subprocess *p_current = nullptr;

    // per process ResultsName, see Process_Integrator::GetMaxEps()
    const std::vector<double> m_epsilon_values = {-5.0, -4.9, -4.8, -4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1, -4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0};
    std::map<std::string, double> m_chosen_alpha_map;
    std::map<std::string, double> m_chosen_efficiency_map;
    std::map<std::string, double> m_chosen_capped_fraction_map;
    std::map<std::string, double> m_xsec_map;
    std::map<std::string, int> m_fills_map;
    std::map<std::string, std::vector<double>> m_wmax_manual_map;
    // generated events per accepted one for unweighted events, whose hard
    // events with an overweight are repeated, per scanned Max_Epsilon
    std::map<std::string, std::vector<double>> m_genperacc_manual_map;
    std::map<std::string, std::vector<double>> m_alpha_manual_map;
    std::map<std::string, std::vector<double>> m_alpha_manual_fraction_map;
    std::map<std::string, std::vector<std::vector<double>>> m_alpha_manual_fscan_map;
    std::map<std::string, std::vector<double>> m_efficiency_manual_map;
    std::map<std::string, std::vector<std::vector<double>>> m_efficiency_manual_fscan_map;

    Timing_Statistics();

    static inline Clock::time_point Now() { return Clock::now(); }
    Subprocess &Current();
    static void Book(Subprocess &sub, const Type type, const Clock::duration &time);

    // the counters and times as name -> value, for the report
    std::map<std::string, double> TimeMapAll() const;
    std::map<std::string, long int> NumberMapAll() const;
  };

}

#endif
