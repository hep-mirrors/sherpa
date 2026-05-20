#include "ATOOLS/Org/CXXFLAGS.H"

#include <mpi.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <functional>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <regex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <unistd.h>

#include <hdf5.h>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>

#include "ATOOLS/Org/Message.H"

#define FIX__BROKEN_EVENT_FILES

using namespace HighFive;

namespace LHEH5 {

  struct ProcInfo {
    int pid, nplo, npnlo;
    double unitwgt, xsec;
    inline ProcInfo(int _pid,int _nplo,int _npnlo,
		    double _unitwgt,double _xsec):
      pid(_pid), nplo(_nplo), npnlo(_npnlo),
      unitwgt(_unitwgt), xsec(_xsec) {}
  };//end of struct ProcInfo

  inline std::ostream &operator<<(std::ostream &s,const ProcInfo &p)
  { return s<<"[pid="<<p.pid<<",nplo="<<p.nplo<<",npnlo="<<p.npnlo
	    <<",unitwgt="<<p.unitwgt<<",xsec="<<p.xsec<<"]"; }

  struct Particle {
    int id, st, mo1, mo2, cl1, cl2;
    double px, py, pz, e, m;
    inline Particle(int _id,int _st,int _mo1,int _mo2,int _cl1,int _cl2,
		    double _px,double _py,double _pz,double _e,double _m):
      id(_id), st(_st), mo1(_mo1), mo2(_mo2), cl1(_cl1), cl2(_cl2),
      px(_px), py(_py), pz(_pz), e(_e), m(_m) {}
  };// end of struct Particle

  inline std::ostream &operator<<(std::ostream &s,const Particle &p)
  { return s<<"{id="<<p.id<<",st="<<p.st
	    <<",mo=["<<p.mo1<<","<<p.mo2<<"]"
	    <<",cl=["<<p.cl1<<","<<p.cl2<<"]"
	    <<",p=("<<p.e<<","<<p.px<<","<<p.py<<","<<p.pz<<")}"; }

  struct Event: public std::vector<Particle> {
    ProcInfo pinfo;
    size_t trials;
    std::vector<double> wgts;
    double mur, muf, muq, aqed, aqcd;
    std::vector<Particle> ctparts;
    int ijt, kt, i, j, k;
    double z1, z2, bbpsw, psw;
    inline Event(const ProcInfo &_pinfo,
		 size_t _trials,std::vector<double> _wgts,
		 double _mur,double _muf,double _muq,
		 double _aqed,double _aqcd):
      pinfo(_pinfo), trials(_trials), wgts(_wgts),
      mur(_mur), muf(_muf), muq(_muq), aqed(_aqed), aqcd(_aqcd),
      ijt(-1), kt(-1), i(-1), j(-1), k(-1),
      z1(0), z2(0), bbpsw(0), psw(0) {}
    inline void AddCTInfo(int _ijt, int _kt,int _i,int _j,int _k,
			  double _z1,double _z2,double _bbw,double _w)
    { ijt=_ijt; kt=_kt, i=_i; j=_j; k=_k;
      z1=_z1; z2=_z2; bbpsw=_bbw; psw=_w; }
  };// end of struct Event

  inline std::ostream &operator<<(std::ostream &s,const Event &e)
  { s<<"Event "<<e.pinfo<<" {\n"
     <<"  trials="<<e.trials<<",weights=("<<e.wgts[0];
    for (size_t i(1);i<e.wgts.size();++i) s<<","<<e.wgts[i];
    s<<")\n  mur="<<e.mur<<", muf="<<e.muf<<", muq="<<e.muq
     <<",aqed="<<e.aqed<<",aqcd="<<e.aqcd<<"\n";
    for (size_t i(0);i<e.size();++i) s<<"  "<<e[i]<<"\n";
    if (!e.ctparts.empty() || e.psw) {
      s<<"  ("<<e.ijt<<","<<e.kt<<")->("<<e.i<<","<<e.j
       <<","<<e.k<<"), z1="<<e.z1<<", z2="<<e.z2
       <<", bbpsw="<<e.bbpsw<<", psw="<<e.psw<<"\n";
      for (size_t i(0);i<e.ctparts.size();++i) s<<"  "<<e.ctparts[i]<<"\n";
    }
    return s<<"}"; }

  class LHEFile {
  private:

    std::vector<int> version;
    std::vector<std::vector<double> > evts, parts, pinfo;
    std::vector<std::vector<double> > ctevts, ctparts;
    std::vector<std::string> wgtnames;

    inline Particle GetParticle(size_t i) const
    {
      return Particle(parts[i][0],parts[i][1],parts[i][2],parts[i][3],
		      parts[i][4],parts[i][5],parts[i][6],parts[i][7],
		      parts[i][8],parts[i][9],parts[i][10]);
    }

    inline Particle GetCTParticle(size_t i) const
    {
      return Particle(-1,-1,-1,-1,-1,-1,ctparts[i][0],ctparts[i][1],
		      ctparts[i][2],ctparts[i][3],-1);
    }

  public:

    inline const std::vector<int> &Version() const { return version; }

    inline double TotalXS() const
    {
      double xs(0.);
      for (int i(0);i<pinfo.size();++i) xs+=pinfo[i][3];
      return xs;
    }

    inline double UnitWeight() const
    {
      double wgt(0.);
      for (int i(0);i<pinfo.size();++i) wgt+=pinfo[i][5];
      return wgt;
    }

    inline const std::vector<std::string> &
    WeightNames() const { return wgtnames; }

    inline size_t NProcesses() const { return pinfo.size(); }
    inline ProcInfo GetProcInfo(const size_t pid) const
    {
      return ProcInfo(pid,pinfo[pid][1],pinfo[pid][2],
		      pinfo[pid][5],pinfo[pid][3]);
    }

    inline size_t NEvents() const { return evts.size(); }
    inline Event GetEvent(size_t i) const
    {
      std::vector<double> wgts(evts[i].begin()+9,evts[i].end());
      Event e(GetProcInfo(evts[i][0]?evts[i][0]-1:0),evts[i][3],wgts,
	      evts[i][6],evts[i][5],evts[i][4],evts[i][7],evts[i][8]);
      double wgt(0.);
      for (std::vector<double>::const_iterator
	     it(wgts.begin());it!=wgts.end();++it) wgt+=std::abs(*it);
      if (!wgt) return e;
      for (int n(0);n<evts[i][1];++n)
	e.push_back(GetParticle(evts[i][2]-evts[0][2]+n));
      if (!ctevts.empty()) {
	e.AddCTInfo(ctevts[i][0],ctevts[i][1],ctevts[i][2],
		    ctevts[i][3],ctevts[i][4],ctevts[i][5],
		    ctevts[i][6],ctevts[i][7],ctevts[i][8]);
	if (ctevts[i][0]>=0 && ctevts[i][1]>=0)
	  for (int n(0);n<evts[i][1]+(ctevts[i][0]>=0?1:0);++n)
	    e.ctparts.push_back(GetCTParticle(evts[i][2]-evts[0][2]+n));
      }
      return e;
    }

    void ReadHeader(File &file)
    {
      auto xfer_props = DataTransferProps{};
      file.getDataSet("version").read(version,xfer_props);
      file.getDataSet("procInfo").read(pinfo,xfer_props);
      DataSet events(file.getDataSet("events"));
      auto attr_keys(events.listAttributeNames());
      Attribute a(events.getAttribute(attr_keys[0]));
      a.read(wgtnames);
      for (int i(0);i<9;++i) wgtnames.erase(wgtnames.begin());
    }
    void ReadEvents(File &file,size_t first_event,size_t n_events)
    {
      auto xfer_props = DataTransferProps{};
      DataSet events(file.getDataSet("events"));
      std::vector<size_t> eoffsets{first_event,0};
      std::vector<size_t> ecounts{n_events,9+wgtnames.size()};
      evts.resize(n_events,std::vector<double>(9+wgtnames.size()));
      events.select(eoffsets,ecounts).read(evts,xfer_props);
      DataSet particles(file.getDataSet("particles"));
      std::vector<size_t> poffsets{(size_t)evts.front()[2],0};
      size_t nmax(0);
      for (size_t i(0);i<pinfo.size();++i)
	nmax=std::max((size_t)std::max(pinfo[i][1],pinfo[i][2]+1),nmax);
      std::vector<size_t> pcounts{n_events*nmax,13};
      parts.resize(n_events*nmax,std::vector<double>(13));
      particles.select(poffsets,pcounts).read(parts,xfer_props);
      if (file.exist("ctevents")) {
	DataSet events(file.getDataSet("ctevents"));
	std::vector<size_t> eoffsets{first_event,0};
	std::vector<size_t> ecounts{n_events,9};
	ctevts.resize(n_events,std::vector<double>(9));
	events.select(eoffsets,ecounts).read(ctevts,xfer_props);
	DataSet particles(file.getDataSet("ctparticles"));
	std::vector<size_t> poffsets{(size_t)evts.front()[2],0};
	std::vector<size_t> pcounts{n_events*nmax,4};
	ctparts.resize(n_events*nmax,std::vector<double>(4));
	particles.select(poffsets,pcounts).read(ctparts,xfer_props);
      }
    }

  };// end of struct LHEFile

  // Enables the HDF5 core (in-memory) virtual file driver on a
  // FileAccessProps. backing_store=0 means nothing is written to disk.
  struct CoreFileAccess {
    size_t increment;
    hbool_t backing_store;
    void apply(hid_t fapl) const {
      if (H5Pset_fapl_core(fapl, increment, backing_store) < 0)
        throw std::runtime_error("H5Pset_fapl_core failed");
    }
  };

}// end of namespace LHEH5

#include "PHASIC++/Main/Event_Reader.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "PDF/Main/PDF_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"

#include <pepper/pepper.h>

using namespace PHASIC;
using namespace ATOOLS;

namespace LHEH5 {

  // Parameters of a Sherpa "Mass" selector entry,
  // i.e. `[Mass, kf1, kf2, min, max]`.
  struct Mass_Selector_Params {
    int kf1;
    int kf2;
    double min;
    double max;
  };

  // Scan the SELECTORS list in the Sherpa main settings and return the
  // parameters of every entry shaped like `[Mass, kf1, kf2, min, max]`.
  // Other selector kinds are ignored here; Pepper only consumes a subset.
  std::vector<Mass_Selector_Params> ReadMassSelectorParams()
  {
    std::vector<Mass_Selector_Params> result;
    auto items = Settings::GetMainSettings()["SELECTORS"].GetItems();
    for (auto& item : items) {
      if (!item.IsList()) continue;
      const auto parameters
          = item.SetDefault<std::string>({}).GetVector<std::string>();
      if (parameters.size() != 5) continue;
      if (parameters[0] != "Mass") continue;
      Mass_Selector_Params p;
      p.kf1 = item.Interprete<int>(parameters[1]);
      p.kf2 = item.Interprete<int>(parameters[2]);
      p.min = item.Interprete<double>(parameters[3]);
      p.max = item.Interprete<double>(parameters[4]);
      // Pepper only implements lepton pair mass cuts.
      if (!Flavour{p.kf1}.IsLepton() || !Flavour{p.kf2}.IsLepton()) {
        continue;
      }
      result.push_back(p);
    }
    return result;
  }

  // Serial FIFO worker shared across all Pepper_Reader instances.
  // Pepper drives a single backend (CPU or GPU); running fills from
  // sibling readers concurrently would multiply GPU memory pressure
  // and re-enter Pepper, which is not safe. A single worker thread
  // processes submitted jobs one at a time, so each reader sees its
  // own future complete, but only one fill is ever in flight.
  class Fill_Worker {
  public:
    Fill_Worker() : m_thread([this] { Run(); }) {}
    ~Fill_Worker()
    {
      {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_stop = true;
      }
      m_cv.notify_all();
      m_thread.join();
    }

    Fill_Worker(const Fill_Worker&) = delete;
    Fill_Worker& operator=(const Fill_Worker&) = delete;

    std::future<size_t> Submit(std::function<size_t()> job)
    {
      auto task = std::make_shared<std::packaged_task<size_t()>>(
          std::move(job));
      auto fut = task->get_future();
      {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_queue.push_back([task]() { (*task)(); });
      }
      m_cv.notify_one();
      return fut;
    }

  private:
    void Run()
    {
      while (true) {
        std::function<void()> job;
        {
          std::unique_lock<std::mutex> lock(m_mutex);
          m_cv.wait(lock, [this] { return m_stop || !m_queue.empty(); });
          if (m_stop && m_queue.empty()) return;
          job = std::move(m_queue.front());
          m_queue.pop_front();
        }
        job();
      }
    }

    std::mutex m_mutex;
    std::condition_variable m_cv;
    std::deque<std::function<void()>> m_queue;
    bool m_stop {false};
    std::thread m_thread;
  };

  // Translate Sherpa's EW input scheme into the equivalent Pepper
  // configuration and forward the relevant numeric parameters. Sherpa
  // exposes many schemes; Pepper only implements `Gmu` and
  // `alphamZsW`, and `UserDefined`, so other choices are rejected here rather
  // than silently mismatching the rest of the run.
  void SyncEWScheme(Pepper::Initialization_settings& settings)
  {
    const auto sherpa_scheme = ToType<MODEL::ew_scheme::code>(
        rpa->gen.Variable("EW_SCHEME"));
    Settings& sherpa_settings {Settings::GetMainSettings()};
    const Complex csin2thetaW {
        MODEL::s_model->ComplexConstant("csin2_thetaW")};
    const double alpha_em {MODEL::aqed->Default()};
    switch (sherpa_scheme) {
      case MODEL::ew_scheme::UserDefined: {
        // All EW parameters given explicitly by the user on the Sherpa
        // side. Pepper's `none` scheme also consumes alpha_em and
        // sin^2(theta_W) verbatim from settings, with no further
        // derivation, so the two are equivalent as long as we propagate
        // both values (handled unconditionally below).
        settings.set_ew_scheme("none");
        break;
      }
      case MODEL::ew_scheme::Gmu: {
        settings.set_ew_scheme("gmu");
        // Pepper's "gmu" always derives sin^2(theta_W) from the complex
        // W/Z masses, so it is equivalent to Sherpa's CMS width scheme.
        // The Gmu_cms_alpha_qed_convention only selects the alpha_QED
        // formula:
        //   - "abs"     matches GMU_CMS_AQED_CONVENTION=0,
        //   - "sherpa2" matches GMU_CMS_AQED_CONVENTION=4 (real-mass
        //     formula sqrt(2) GF / pi * MW^2 * (1 - (MW/MZ)^2)).
        // WIDTH_SCHEME=Fixed cannot be mapped exactly: Sherpa keeps
        // sin^2(theta_W) real there, while Pepper unconditionally uses
        // the complex-mass relation, so the propagated EW parameters
        // would silently disagree.
        const std::string width_scheme {
            sherpa_settings["WIDTH_SCHEME"].Get<std::string>()};
        if (width_scheme != "CMS") {
          THROW(not_implemented,
                "Pepper's Gmu scheme requires WIDTH_SCHEME=CMS; got "
                + width_scheme + ".");
        }
        const size_t conv {
            sherpa_settings["GMU_CMS_AQED_CONVENTION"].Get<size_t>()};
        switch (conv) {
          case 0:
            settings.set_gmu_cms_alpha_qed_convention("abs");
            break;
          case 4:
            settings.set_gmu_cms_alpha_qed_convention("sherpa2");
            break;
          default:
            THROW(not_implemented,
                  "Pepper supports GMU_CMS_AQED_CONVENTION=0 (abs) and 4 "
                  "(sherpa2); got " + ToString(conv) + ".");
        }
        settings.set_GF(sherpa_settings["GF"].Get<double>());
        break;
      }
      case MODEL::ew_scheme::alphamZsW: {
        settings.set_ew_scheme("alpha_mz_sinthetaw");
        break;
      }
      default:
        THROW(not_implemented,
              "Pepper supports only EW_SCHEME=UserDefined (0), Gmu (3), "
              "or alphamZsW (4); got " + ToString(sherpa_scheme) + ".");
    }
    settings.set_alpha_em(alpha_em);
    settings.set_sin2_theta_w(csin2thetaW.real());
  }

  // Push the heavy-particle pole masses and decay widths that Pepper
  // consumes (top, EW gauge bosons, Higgs, b, c, tau) from Sherpa's
  // particle data table into Pepper's initialization settings, so that
  // the two sides agree on the parameter point without the user having
  // to mirror every value in the Pepper run card. File-based Pepper
  // settings still win.
  void SyncParticleProperties(Pepper::Initialization_settings& settings)
  {
    static constexpr int kfs[] = {kf_t, kf_Z, kf_Wplus};
    for (const int kf : kfs) {
      const ATOOLS::Flavour fl{static_cast<long int>(kf)};
      settings.set_particle_mass(kf, fl.Mass());
      settings.set_particle_width(kf, fl.Width());
    }
  }

  // Try to map Sherpa's SCALES string onto one of Pepper's hard-coded
  // μ² scale setters. Pepper only ships a handful (m_Z^2, m_W^2,
  // m_t^2, H_Tp^2, H_Tp^2/2, H_T^2/2, H_TM^2/2), so we recognise the
  // two simplest Sherpa forms: VAR{sqr(<mass>)}{sqr(<mass>)} for the
  // fixed-mass cases, and VAR{H_X[/2]}{H_X[/2]} for the H_T-style
  // ones where the divisor matches Pepper's hard-coded /2. Sherpa's
  // H_T scales can be normalised with other factors (e.g. /4); those
  // fall through to the warning case below. Returns the Pepper spec
  // string on success, or an empty string when no clean match was
  // found.
  std::string MapSherpaScalesToPepperSpec(const std::string& sherpa_scales)
  {
    // Sherpa's bare-named core scale setters. "Default" replicates
    // Sherpa's Default_Core_Scale on the Pepper side, currently only
    // implemented for hh collisions.
    static const std::regex default_pattern{R"(^\s*Default\s*$)"};
    if (std::regex_match(sherpa_scales, default_pattern))
      return "Sherpa_default_hh";

    // Accept either VAR{X} or VAR{X}{X}: Pepper carries a single μ²,
    // so a one-arg form is unambiguous, and a two-arg form is only
    // synced when μ_F and μ_R agree (otherwise we'd silently collapse
    // a genuine factorisation/renormalisation split).
    static const std::regex var_pattern{
      R"(^\s*VAR\s*\{\s*([^{}]*?)\s*\}(?:\s*\{\s*([^{}]*?)\s*\})?)"};
    std::smatch m;
    if (!std::regex_search(sherpa_scales, m, var_pattern)) return {};
    const std::string arg{m[1].str()};
    if (m[2].matched && arg != m[2].str()) return {};

    // Fixed mass squared: sqr(N). Compare against the Sherpa pole
    // masses for the bosons/top that Pepper hard-codes.
    static const std::regex sqr_pattern{
      R"(^sqr\(\s*([0-9.eE+-]+)\s*\)$)"};
    std::smatch sm;
    if (std::regex_match(arg, sm, sqr_pattern)) {
      const double mu2 = std::stod(sm[1].str());
      auto close = [&](double mass) {
        return std::abs(mu2 - mass) < 1e-2 * mass;
      };
      if (close(ATOOLS::Flavour(kf_Z).Mass()))     return "m_Z^2";
      if (close(ATOOLS::Flavour(kf_Wplus).Mass())) return "m_W^2";
      if (close(ATOOLS::Flavour(kf_t).Mass()))     return "m_t^2";
      return {};
    }

    // Dynamic H_T-like scales: VAR{H_X[/2]}. Pepper offers H_Tp^2,
    // H_Tp^2/2, H_T^2/2 and H_TM^2/2 — no half-less variants of H_T
    // or H_TM.
    static const std::regex ht_pattern{
      R"(^(H_Tp2|H_T2|H_TM2)\s*(/\s*2)?$)"};
    std::smatch hm;
    if (std::regex_match(arg, hm, ht_pattern)) {
      const std::string var{hm[1].str()};
      const bool halved{hm[2].matched};
      if (var == "H_Tp2") return halved ? "H_Tp^2/2" : "H_Tp^2";
      if (var == "H_T2"  && halved) return "H_T^2/2";
      if (var == "H_TM2" && halved) return "H_TM^2/2";
      return {};
    }

    return {};
  }

  // Resolve the SCALES string that actually drives the (unclustered)
  // hard process. With SCALES: METS, the shower decides emission
  // scales but the core process — which is all Pepper sees — uses
  // CORE_SCALE instead. Honour an inline `METS[C:VAR{...}]` override
  // first, then fall back to the global MEPS.CORE_SCALE setting.
  std::string ResolveSherpaCoreScale(const std::string& sherpa_scales)
  {
    static const std::regex mets_prefix{R"(^\s*METS\b)"};
    if (!std::regex_search(sherpa_scales, mets_prefix)) return sherpa_scales;
    static const std::regex inline_c{
      R"(\[[^\]]*\bC\s*:\s*([^\s,\]]+))"};
    std::smatch m;
    if (std::regex_search(sherpa_scales, m, inline_c)) return m[1].str();
    return ATOOLS::Settings::GetMainSettings()["MEPS"]["CORE_SCALE"]
      .GetScalarWithOtherDefault<std::string>("Default");
  }

  // Align Pepper's μ² scale setter with Sherpa's SCALES setting when
  // we can recognise it; otherwise warn. Pepper's scale choice does
  // not change physics, since Sherpa reweights each event to its own
  // scale on the fly, but a large mismatch reduces Pepper's
  // unweighting efficiency.
  void SyncScaleSetter(Pepper::Initialization_settings& settings)
  {
    const std::string sherpa_scales{
      ATOOLS::Settings::GetMainSettings()["SCALES"].Get<std::string>()};
    const std::string effective{ResolveSherpaCoreScale(sherpa_scales)};
    const std::string pepper_spec{MapSherpaScalesToPepperSpec(effective)};
    if (!pepper_spec.empty()) {
      settings.set_scale_setter(pepper_spec);
      msg_Info() << "Pepper_Interface: scale setter synced to '"
                 << pepper_spec << "' (from SCALES=\""
                 << sherpa_scales << "\"";
      if (effective != sherpa_scales)
        msg_Info() << ", core scale \"" << effective << "\"";
      msg_Info() << ").\n";
      return;
    }
    msg_Info() << om::brown << om::bold << "WARNING" << om::reset
               << ": Could not map Sherpa SCALES=\"" << sherpa_scales << "\"";
    if (effective != sherpa_scales)
      msg_Info() << " (core scale \"" << effective << "\")";
    msg_Info() << " onto one of Pepper's hard-coded scale setters "
                  "(m_Z^2, m_W^2, m_t^2, H_Tp^2, H_Tp^2/2, H_T^2/2, H_TM^2/2, "
		  "Sherpa_default_hh). This does NOT affect physics, "
                  "since Sherpa reweights each event to its own scale, "
                  "but Pepper's unweighting efficiency may be reduced. "
                  "Please contact the Pepper authors if you need a "
                  "closer match for your setup.\n";
  }

  // Map a signed Sherpa/PDG KF code onto the particle-name string that
  // Pepper's process-spec parser expects (see Pepper's
  // `string_from_particle`). Sherpa's Flavour::IDName disagrees with
  // Pepper on a handful of names (e.g. "G" vs "g", "W+" vs "W"), so we
  // do the conversion off the numeric code instead.
  std::string KfcodeToPepperName(long int kf)
  {
    const long int akf = std::abs(kf);
    if (akf == kf_gluon) return "g";
    if (akf == kf_Z)     return "Z";
    if (akf == kf_Wplus) return (kf > 0) ? "W" : "W-";
    std::string base;
    bool has_charge = false; // charged leptons get explicit +/- in Pepper
    switch (akf) {
      case kf_d:     base = "d";    break;
      case kf_u:     base = "u";    break;
      case kf_s:     base = "s";    break;
      case kf_c:     base = "c";    break;
      case kf_b:     base = "b";    break;
      case kf_t:     base = "t";    break;
      case kf_e:     base = "e";    has_charge = true; break;
      case kf_mu:    base = "mu";   has_charge = true; break;
      case kf_tau:   base = "tau";  has_charge = true; break;
      // Pepper's `particle_from_string` accepts the two-character
      // "v<flavour>" forms but errors on the longer "vmu"/"vtau" forms
      // that `string_from_particle` happens to emit, so we use the
      // round-tripping forms here.
      case kf_nue:   base = "ve";   break;
      case kf_numu:  base = "vm";   break;
      case kf_nutau: base = "vt";   break;
      default:
        THROW(not_implemented,
              "Cannot translate KF code " + ToString(kf)
              + " to a Pepper particle name.");
    }
    if (kf >= 0) {
      if (has_charge) base += "-";
      return base;
    }
    return has_charge ? base + "+" : base + "b";
  }

  // Build the Pepper process-spec string from the signed KF codes that
  // Process_Base forwards via Event_Reader_Key.
  //
  // Two cases:
  //  - Sherpa's "jet" container (kf=93) on both initial-state lines is
  //    Pepper's "pp" beam: we emit a compound name like "ppee", "ppev",
  //    "pptt", "ppjj", ..., which Pepper resolves to a process-data file
  //    that already sums over the partonic channels. Outgoing kf=93s are
  //    counted and turn into trailing 'j' characters.
  //  - Otherwise we fall back to a partonic-channel spec like
  //    "u ub -> Z g" assembled via `KfcodeToPepperName`.
  std::string BuildPepperProcessSpec(const std::vector<long int>& flavours,
                                     std::size_t nin)
  {
    if (flavours.size() < nin + 1 || nin == 0)
      THROW(invalid_input,
            "Pepper_Reader: cannot build a process spec from "
            + ToString(flavours.size()) + " flavours with "
            + ToString(nin) + " incoming.");

    if (nin == 2 && flavours[0] == kf_jet && flavours[1] == kf_jet) {
      std::size_t n_out_jets = 0;
      std::vector<long int> rest;
      for (std::size_t i = nin; i < flavours.size(); ++i) {
        if (flavours[i] == kf_jet) ++n_out_jets;
        else rest.push_back(flavours[i]);
      }
      std::sort(rest.begin(), rest.end());
      std::string token;
      if (rest.empty()) {
        token = "";
      } else if (rest == std::vector<long int>{-11, 11}) {
        token = "ee";
      } else if (rest == std::vector<long int>{-12, 12}) {
        token = "vv";
      } else if (rest == std::vector<long int>{-12, 11}
                 || rest == std::vector<long int>{-11, 12}) {
        token = "ev";
      } else if (rest == std::vector<long int>{-6, 6}) {
        token = "tt";
      } else {
        std::string outgoing;
        for (auto kf : rest) {
          if (!outgoing.empty()) outgoing += ' ';
          outgoing += ToString(kf);
        }
        THROW(not_implemented,
              "Cannot map Sherpa final-state flavours {" + outgoing
              + "} (alongside kf_jet=93 incoming) onto a Pepper compound "
              "process name. Supported compound finals are: (empty), "
              "e+ e-, ve veb, e ve, t tb.");
      }
      std::string spec {"pp" + token};
      spec.append(n_out_jets, 'j');
      return spec;
    }

    std::string spec;
    for (std::size_t i = 0; i < flavours.size(); ++i) {
      if (i > 0) spec += ' ';
      if (i == nin) spec += "-> ";
      spec += KfcodeToPepperName(flavours[i]);
    }
    return spec;
  }

  // Owns the once-per-process Pepper library bring-up and teardown.
  // A static weak_ptr deduplicates instances across Pepper_Reader
  // siblings (one per jet multiplicity); finalize fires when the
  // last reader holding a shared_ptr is destroyed, i.e. during
  // Sherpa's normal shutdown rather than at static-destructor time.
  class Pepper_Interface {
  public:
    Pepper_Interface()
    {
      // We do not currently plumb Sherpa's argc/argv through to
      // Pepper; a dummy pair is sufficient since MPI/Kokkos init
      // do not need them once Pepper's MPI bring-up is disabled.
      static int s_argc = 0;
      static char *s_argv[] = {nullptr};
      Pepper::Initialization_settings settings(s_argc, s_argv);
      settings.disable_mpi_initialization();

      // Keep Pepper's beam energy in lockstep with Sherpa's so the two
      // sides do not silently disagree on √s. Sherpa's Run_Parameter is
      // populated during framework initialisation, which precedes
      // Matrix_Element_Handler::InitializeProcesses (where readers are
      // constructed), so rpa->gen.Ecms() is available here.
      settings.set_e_cms(rpa->gen.Ecms());

      // Pin Pepper's LHAPDF choice to whatever Sherpa initialised for
      // the hard process (beam 0); the second beam uses the same set
      // in all configurations Pepper supports.
      if (auto* pdf = rpa->gen.PDF(0))
        settings.set_pdf(pdf->Set(), pdf->Member());

      // Mirror Sherpa's EW input scheme and the corresponding numeric
      // parameters. Pepper only implements `Gmu`, `alphamZsW`, and
      // `UserDefined`, so any other Sherpa scheme is a hard error.
      SyncEWScheme(settings);

      // Forward Sherpa's heavy-particle masses and widths so Pepper
      // evaluates matrix elements at the same parameter point.
      SyncParticleProperties(settings);

      // Best-effort alignment of Pepper's μ² scale setter with
      // Sherpa's SCALES choice; mismatches are warned about, not
      // fatal, because they only affect Pepper unweighting efficiency.
      SyncScaleSetter(settings);

      // Translate Sherpa "Mass" selectors into Pepper's lepton-pair
      // invariant-mass cut. Pepper exposes a single (m2_min, m2_max)
      // pair, so when multiple Mass selectors are present we narrow to
      // their intersection.
      const auto mass_selectors = ReadMassSelectorParams();
      double m2_min = 0.0;
      double m2_max = std::numeric_limits<double>::infinity();
      for (const auto& p : mass_selectors) {
        m2_min = std::max(m2_min, p.min * p.min);
        m2_max = std::min(m2_max, p.max * p.max);
        msg_Debugging() << "Pepper_Interface: Mass selector "
                        << "[kf1=" << p.kf1 << ",kf2=" << p.kf2
                        << ",min=" << p.min << ",max=" << p.max << "]\n";
      }
      if (!mass_selectors.empty()) {
        settings.set_ll_m2_min(m2_min);
        settings.set_ll_m2_max(m2_max);
      }

      Pepper::initialize(settings);
    }

    // The worker is destroyed (and joined) before Pepper::finalize()
    // runs, so no in-flight fill outlives Pepper itself.
    ~Pepper_Interface() { m_fill_worker.reset(); Pepper::finalize(); }

    Pepper_Interface(const Pepper_Interface &) = delete;
    Pepper_Interface &operator=(const Pepper_Interface &) = delete;

    Fill_Worker &fill_worker() { return *m_fill_worker; }

    static std::shared_ptr<Pepper_Interface> Acquire()
    {
      static std::mutex s_mutex;
      static std::weak_ptr<Pepper_Interface> s_instance;
      std::lock_guard<std::mutex> lock(s_mutex);
      if (auto existing = s_instance.lock()) return existing;
      auto fresh = std::make_shared<Pepper_Interface>();
      s_instance = fresh;
      return fresh;
    }

  private:
    std::unique_ptr<Fill_Worker> m_fill_worker {
      std::make_unique<Fill_Worker>()};
  };

  // Reads parton-level events from HDF5 databases that live entirely
  // in RAM and are filled by the Pepper library. A double-buffered
  // current/next pair lets consumption and refill overlap: events are
  // drawn from "current", while "next" is filled ahead of time so it is
  // ready to take over once "current" is depleted.
  class Pepper_Reader: public Event_Reader, public MPI_Object {
  private:

    std::shared_ptr<Pepper_Interface> p_pepper;
    std::unique_ptr<Pepper::Process> p_process;

    std::unique_ptr<File> m_current_file, m_next_file;
    std::unique_ptr<LHEFile> p_current, p_next;

    size_t m_ievt, m_iblock, m_trials, m_ncurrent, m_ncache;
    // Running total of events consumed in buffers that have already
    // been rotated out by Advance(). Combined with m_ievt this gives
    // the lifetime number of events served by the reader, which is
    // what PrintStatistics reports at end-of-run.
    size_t m_total_consumed {0};
    bool m_finished;
    bool m_warmed_up;
    bool m_finalized {false};

    // When true, refilling the next buffer runs on a worker thread so
    // Sherpa can consume events from the current buffer while Pepper
    // populates the next one (especially helpful with a GPU back-end).
    // Invariant: while m_pending_fill is valid, the main thread must not
    // touch m_next_file or p_next; only one worker is ever in flight, so
    // Pepper is never re-entered concurrently.
    bool m_async_fill;
    std::future<size_t> m_pending_fill;
    size_t m_last_fill_result {0};

    Vec4D_Vector m_ctmoms;

    std::unique_ptr<File> CreateInMemoryFile(const std::string &tag)
    {
      auto fapl = FileAccessProps::Empty();
      fapl.add(CoreFileAccess{1u << 20, /*backing_store=*/0});
      return std::make_unique<File>(
          tag, File::Truncate | File::Create | File::ReadWrite, fapl);
    }

    // HDF5 keys files by their string name in a global registry, even
    // for the in-memory core driver. Two simultaneously-live files with
    // the same name would collide both within one reader (current and
    // next buffers overlap during Advance) and across sibling readers
    // (one per jet multiplicity). A process-wide atomic counter gives
    // every buffer a unique name.
    std::string NextBufferName()
    {
      static std::atomic<std::uint64_t> s_counter {0};
      return "pepper_buffer_" + std::to_string(s_counter.fetch_add(1));
    }

    // Ask Pepper to populate the next in-memory database with up to
    // m_ncache events. Returning 0 signals Pepper exhaustion; the
    // reader stops once the current buffer is drained. This routine
    // is the unit of work scheduled by LaunchFill() and may run on a
    // worker thread when PEPPER_ASYNC_FILL is enabled.
    size_t PerformFill()
    {
      m_next_file = CreateInMemoryFile(NextBufferName());
      const size_t n_filled = p_process->fill_lheh5_buffer(*m_next_file, m_ncache);
      if (n_filled == 0) {
        m_next_file.reset();
        p_next.reset();
        return 0;
      }
      p_next = std::make_unique<LHEFile>();
      p_next->ReadHeader(*m_next_file);
      p_next->ReadEvents(*m_next_file, 0, n_filled);
      return n_filled;
    }

    // Kick off a fill of the next buffer. With async fill enabled, the
    // work is submitted to the shared Fill_Worker — a single FIFO thread
    // serving every Pepper_Reader instance — so that sibling readers do
    // not stack concurrent Pepper invocations onto the (often GPU)
    // backend. The result becomes visible via WaitForFill(); otherwise
    // it runs inline and m_last_fill_result is updated directly.
    // Precondition: no fill is currently in flight for this reader.
    void LaunchFill()
    {
      if (m_async_fill) {
        m_pending_fill = p_pepper->fill_worker().Submit(
            [this]() { return PerformFill(); });
      } else {
        m_last_fill_result = PerformFill();
      }
    }

    // Block until the most recently launched fill has completed and
    // return the number of events written. Safe to call when no fill
    // is in flight: it then just returns the cached previous result.
    size_t WaitForFill()
    {
      if (m_pending_fill.valid())
        m_last_fill_result = m_pending_fill.get();
      return m_last_fill_result;
    }

    void Advance()
    {
      // Make sure the next buffer is completely filled before swapping;
      // this is the invariant that keeps consumers from racing the
      // producer thread.
      const size_t filled = WaitForFill();
      // Accumulate consumption of the buffer we are about to discard
      // (zero on the very first Advance from WarmUp).
      m_total_consumed += m_ievt;
      m_current_file = std::move(m_next_file);
      p_current = std::move(p_next);
      m_ievt = 0;
      m_ncurrent = p_current ? p_current->NEvents() : 0;
      ++m_iblock;
      if (m_ncurrent > 0) {
        m_totalxs = p_current->TotalXS();
        m_unitwgt = p_current->UnitWeight();
        if (p_current->Version()[0]==2 &&
            p_current->Version()[1]==0 &&
            p_current->Version()[2]==0) m_unitwgt*=rpa->Picobarn();
      }
      if (filled == 0) m_finished = true;
      else LaunchFill();
    }

    bool FillAmplitude()
    {
      Event e(p_current->GetEvent(m_ievt));
      if (p_current->Version()[0]==2 &&
	  p_current->Version()[1]==0 &&
	  p_current->Version()[2]==0) e.pinfo.unitwgt*=rpa->Picobarn();
      msg_Debugging()<<e<<"\n";
      if (e.empty()) {
	m_trials+=e.trials;
	m_ievt++;
	return false;
      }
      if (e[0].pz<0 && e[1].pz>0) {
	msg_Debugging()<<"Flip initial states\n";
	std::swap<Particle>(e[0],e[1]);
	std::swap<double>(e.z1,e.z2);
	if (e.ctparts.size()) {
	  std::swap<Particle>(e.ctparts[0],e.ctparts[1]);
	  if (e.ijt<2) e.ijt=1-e.ijt;
	  if (e.kt<2) e.kt=1-e.kt;
	  if (e.i<2) e.i=1-e.i;
	  if (e.k<2) e.k=1-e.k;
	}
      }
      p_ampl = Cluster_Amplitude::New();
      for (size_t i(0);i<e.size();++i) {
	Flavour fl((long int)(e[i].id));
	Vec4D p(e[i].e,e[i].px,e[i].py,e[i].pz);
	ColorID cl(i<2?e[i].cl2:e[i].cl1,
		   i<2?e[i].cl1:e[i].cl2);
	p_ampl->CreateLeg(i<2?-p:p,i<2?fl.Bar():fl,cl);
      }
      p_ampl->SetNIn(2);
      p_ampl->SetMuR2(sqr(e.mur));
      p_ampl->SetMuF2(sqr(e.muf));
      p_ampl->SetMuQ2(sqr(e.muq));
      p_ampl->SetKT2(sqr(e.muq));
      p_ampl->SetLKF(e.wgts[0]);
      m_compute=1;
      m_trials+=e.trials;
      if (e.pinfo.npnlo>0) {
	p_ampl->SetLKF(e.psw);
	m_compute=2;
      }
      if (e.ctparts.empty()) p_sub->m_n=0;
      else {
	m_ctmoms.resize(p_sub->m_n=e.ctparts.size());
	for (size_t i(0);i<e.ctparts.size();++i)
	  m_ctmoms[i]=Vec4D(e.ctparts[i].e,e.ctparts[i].px,
			    e.ctparts[i].py,e.ctparts[i].pz);
	p_sub->p_mom=&m_ctmoms[0];
	p_sub->m_ijt=e.ijt;
	p_sub->m_kt=e.kt;
	p_sub->m_i=e.i;
	p_sub->m_j=e.j;
	p_sub->m_k=e.k;
	p_sub->m_x1=e.z1;
	p_sub->m_x2=e.z2;
	p_sub->m_result=e.bbpsw;
      }
      return true;
    }

  public:

    Pepper_Reader(const Event_Reader_Key &key):
      Event_Reader(key), p_pepper(Pepper_Interface::Acquire()),
      m_ievt(0), m_iblock(0), m_trials(0),
      m_ncurrent(0), m_ncache(0), m_finished(false), m_warmed_up(false)
    {
      Settings& s {Settings::GetMainSettings()};
      m_ncache = s["PEPPER_CACHE_SIZE"].SetDefault(10000).Get<int>();
      m_async_fill = s["PEPPER_ASYNC_FILL"].SetDefault(false).Get<bool>();

      p_sub = new NLO_subevt();
      s_objects.push_back(this);

      // The user may either pass an explicit Pepper process spec via
      // `Event_Source: Pepper[<spec>]`, or write a bare `Event_Source:
      // Pepper` and let us derive it from the surrounding Sherpa process
      // (forwarded as KF codes in `m_flavours` by Process_Base).
      const std::string spec {m_files.empty()
                                  ? BuildPepperProcessSpec(m_flavours, m_nin)
                                  : m_files.front()};
      p_process = std::make_unique<Pepper::Process>(spec);

      // The pipeline is primed lazily in WarmUp(), so that Pepper's
      // integration-grid optimisation and weight-maximum determination
      // run from the same hook that Sherpa uses for its own warm-up,
      // i.e. Matrix_Element_Handler::CalculateTotalXSecs.

      // Guard against a setting that will, in general, cause synchronization
      // to break in MPI runs, such that the entire run stalls.
      bool m_printmpixs {s["PRINT_MPI_XS"].Get<bool>()};
      if (m_printmpixs && mpi->MySize() > 1) {
        THROW(invalid_input,
              "`PRINT_MPI_XS: true` can not be used when reading parton-level "
              "events from Pepper. Please use `PRINT_MPI_XS: false`.");
      }
    }

    ~Pepper_Reader()
    {
      // Drain any in-flight worker before our members go out of scope, since
      // the worker captures `this` and touches p_process / the next-buffer
      // state. Normally Sherpa has already done this via Finalize(); this is a
      // fallback for potential shutdown paths that bypass it.
      Finalize();
      if (p_ampl) p_ampl->Delete();
      delete p_sub;
    }

    // Drain the async fill worker and release every Pepper-backed
    // resource while Pepper / Kokkos / MPI are still alive. Called from
    // Sherpa's shutdown path before any of those libraries finalise,
    // and also as a fallback from the destructor. Idempotent.
    void Finalize() override
    {
      if (m_finalized) return;
      m_finalized = true;
      if (m_pending_fill.valid()) m_pending_fill.wait();
      p_current.reset();
      p_next.reset();
      m_current_file.reset();
      m_next_file.reset();
      p_process.reset();
      // Releasing the last shared_ptr triggers ~Pepper_Interface, which
      // joins the fill worker and calls Pepper::finalize(); we want that
      // to happen here, not from an atexit handler after Kokkos is gone.
      p_pepper.reset();
    }

    // Drive Pepper's integration-grid optimisation and unit-weight
    // determination. Called once per process from Sherpa's cross-section
    // calculation step. Pepper performs its optimisation phases lazily
    // inside fill_lheh5_buffer on first use, so priming the double-buffered
    // pipeline here implicitly drives warm-up. Cached Pepper results are
    // reused on subsequent runs (Pepper handles its own cache).
    void WarmUp() override
    {
      if (m_warmed_up) return;
      m_warmed_up = true;
      // First fill is always the warm-up; do it synchronously so the
      // optimisation phase completes before Sherpa proceeds. Advance()
      // will then promote it to the current buffer and kick off the
      // first overlapped refill (async, if enabled).
      const bool saved_async = m_async_fill;
      m_async_fill = false;
      LaunchFill();
      m_async_fill = saved_async;
      if (m_last_fill_result == 0) m_finished = true;
      else Advance();
    }

    void PrintStatistics(std::ostream& o) override
    {
      const size_t total_consumed = m_total_consumed + m_ievt;
      o << "    Pepper: " << total_consumed << " events read"
        << " (current in-memory buffer #" << m_iblock << ": " << m_ievt
        << " of " << m_ncurrent;
      if (m_ncurrent > 0)
        o << ", " << m_ievt * 1000 / m_ncurrent / 10.0 << " %";
      o << ")";
      if (mpi->MySize() > 1) o << " on rank 0";
      o << '\n';
    }

    // Pepper feeds events independently on every rank, so there is
    // no shared file pointer to synchronize. The barrier is the right
    // moment to retry a refill if a previous attempt came up empty.
    // We must not race an in-flight async fill, so first drain any
    // pending future before inspecting p_next.
    void MPISync() override
    {
      WaitForFill();
      if (m_finished && !p_next) {
        // Retry inline (synchronously): we are at a barrier, so there
        // is no consumption to overlap with.
        const bool saved_async = m_async_fill;
        m_async_fill = false;
        LaunchFill();
        m_async_fill = saved_async;
        if (m_last_fill_result > 0) m_finished = false;
      }
    }

    void StepBackward() override
    {
      --m_ievt;
      if (!FillAmplitude())
	THROW(fatal_error, "Could not fill amplitude.");
      m_trials = 1;
    }

    Cluster_Amplitude *ReadEvent() override
    {
      DEBUG_FUNC("i="<<m_ievt<<"("<<m_ncurrent<<"),trial="<<m_trials);
      if (m_trials==0 && p_ampl!=NULL) {
	++m_trials;
	--m_ievt;
      }
      if (m_ievt >= m_ncurrent) {
        // p_next is owned by the worker thread until the fill completes;
        // drain before inspecting it from the consumer side.
        WaitForFill();
        if (!p_next || p_next->NEvents() == 0) {
          if (Communicate() < 0 || (!p_next && m_finished)) {
            msg_Info() << om::brown << om::bold << "WARNING" << om::reset
                       << ": Pepper has no more events to provide. "
                       << "Will stop event generation ...\n";
            rpa->gen.SetNumberOfEvents(rpa->gen.NumberOfGeneratedEvents());
            return NULL;
          }
        }
        Advance();
        if (m_ncurrent == 0) {
          msg_Info() << om::brown << om::bold << "WARNING" << om::reset
                     << ": Pepper has no more events to provide. "
                     << "Will stop event generation ...\n";
          rpa->gen.SetNumberOfEvents(rpa->gen.NumberOfGeneratedEvents());
          return NULL;
        }
      }
      if (p_ampl==NULL) {
	if (!FillAmplitude())
	  return NULL;
      }
      if (m_trials==1) {
	Cluster_Amplitude *ampl(p_ampl);
	p_ampl=NULL;
	--m_trials;
	m_ievt++;
	return ampl;
      }
      --m_trials;
      return NULL;
    }

  };// end of class Pepper_Reader

}// end of namespace LHEH5

using namespace LHEH5;

DECLARE_GETTER(Pepper_Reader,"Pepper",Event_Reader,Event_Reader_Key);

Event_Reader *ATOOLS::Getter<Event_Reader,Event_Reader_Key,Pepper_Reader>::
operator()(const Event_Reader_Key &args) const
{
  return new Pepper_Reader(args);
}

void ATOOLS::Getter<Event_Reader,Event_Reader_Key,Pepper_Reader>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Pepper reader (in-memory HDF5)";
}
