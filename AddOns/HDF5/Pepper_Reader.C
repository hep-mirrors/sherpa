#include "ATOOLS/Org/CXXFLAGS.H"

#include <mpi.h>
#include <iostream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
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

    ~Pepper_Interface() { Pepper::finalize(); }

    Pepper_Interface(const Pepper_Interface &) = delete;
    Pepper_Interface &operator=(const Pepper_Interface &) = delete;

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
    bool m_finished;

    Vec4D_Vector m_ctmoms;

    std::unique_ptr<File> CreateInMemoryFile(const std::string &tag)
    {
      auto fapl = FileAccessProps::Empty();
      fapl.add(CoreFileAccess{1u << 20, /*backing_store=*/0});
      return std::make_unique<File>(
          tag, File::Truncate | File::Create | File::ReadWrite, fapl);
    }

    // Ask Pepper to populate the next in-memory database with up to
    // m_ncache events. Returning 0 signals Pepper exhaustion; the
    // reader stops once the current buffer is drained.
    size_t FillNextBuffer()
    {
      m_next_file = CreateInMemoryFile("pepper_buffer_next");
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

    void Advance()
    {
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
      if (FillNextBuffer() == 0) m_finished = true;
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
      m_ncurrent(0), m_ncache(0), m_finished(false)
    {
      Settings& s {Settings::GetMainSettings()};
      m_ncache = s["PEPPER_CACHE_SIZE"].SetDefault(10000).Get<int>();

      p_sub = new NLO_subevt();
      s_objects.push_back(this);

      if (m_files.empty())
        THROW(invalid_input,
              "Pepper_Reader requires a process specification "
              "(e.g. \"p p -> j j\") as the first entry in EVENT_INPUT.");
      p_process = std::make_unique<Pepper::Process>(m_files.front());

      // Prime the pipeline: fill the first buffer, promote it to
      // current, then fill the look-ahead buffer for the next swap.
      if (FillNextBuffer() == 0) {
        m_finished = true;
      } else {
        Advance();
      }

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
      if (p_ampl) p_ampl->Delete();
      delete p_sub;
    }

    void PrintStatistics(std::ostream& o)
    {
      o << "    Pepper in-memory buffer #" << m_iblock << ": " << m_ievt
        << " of " << m_ncurrent << " events read";
      if (m_ncurrent > 0)
        o << " (" << m_ievt * 1000 / m_ncurrent / 10.0 << " %)";
      if (mpi->MySize() > 1) o << " on rank 0";
      o << '\n';
    }

    // Pepper feeds events independently on every rank, so there is
    // no shared file pointer to synchronize. The barrier is the right
    // moment to retry a refill if a previous attempt came up empty.
    void MPISync()
    {
      if (m_finished && !p_next) {
        if (FillNextBuffer() > 0) m_finished = false;
      }
    }

    void StepBackward()
    {
      --m_ievt;
      if (!FillAmplitude())
	THROW(fatal_error, "Could not fill amplitude.");
      m_trials = 1;
    }

    Cluster_Amplitude *ReadEvent()
    {
      DEBUG_FUNC("i="<<m_ievt<<"("<<m_ncurrent<<"),trial="<<m_trials);
      if (m_trials==0 && p_ampl!=NULL) {
	++m_trials;
	--m_ievt;
      }
      if (m_ievt >= m_ncurrent) {
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
