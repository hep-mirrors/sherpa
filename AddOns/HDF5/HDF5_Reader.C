#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/My_MPI.H"
#ifdef USING__HDF5
#ifdef USING__MPI

#include <mpi.h>
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>

#include <highfive/H5File.hpp>
#include <highfive/H5FileDriver.hpp>
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

  std::ostream &operator<<(std::ostream &s,const ProcInfo &p)
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

  std::ostream &operator<<(std::ostream &s,const Particle &p)
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

  std::ostream &operator<<(std::ostream &s,const Event &e)
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

    inline double TotalXS() const
    {
      double xs(0.);
      for (int i(0);i<pinfo.size();++i) xs+=pinfo[i][3];
      return xs;
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
      xfer_props.add(UseCollectiveIO{});

      file.getDataSet("version").read(version, xfer_props);
      file.getDataSet("procInfo").read(pinfo, xfer_props);
      DataSet init(file.getDataSet("init"));
      DataSet events(file.getDataSet("events"));
      auto attr_keys(events.listAttributeNames());
      Attribute a(events.getAttribute(attr_keys[0]));
      a.read(wgtnames);
      for (int i(0);i<9;++i) wgtnames.erase(wgtnames.begin());
    }
    void ReadEvents(File &file,size_t first_event,size_t n_events)
    {
      auto xfer_props = DataTransferProps{};
      xfer_props.add(UseCollectiveIO{});

      DataSet events(file.getDataSet("events"));
      std::vector<size_t> eoffsets{first_event,0};
      std::vector<size_t> ecounts{n_events,9+wgtnames.size()};
      evts.resize(n_events,std::vector<double>(9+wgtnames.size()));
      events.select(eoffsets,ecounts).read(evts, xfer_props);
      DataSet particles(file.getDataSet("particles"));
      std::vector<size_t> poffsets{(size_t)evts.front()[2],0};
      size_t nmax(0);
      for (size_t i(0);i<pinfo.size();++i)
	nmax=std::max((size_t)std::max(pinfo[i][1],pinfo[i][2]+1),nmax);
      std::vector<size_t> pcounts{n_events*nmax,13};
      parts.resize(n_events*nmax,std::vector<double>(13));
      particles.select(poffsets,pcounts).read(parts, xfer_props);
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

}// end of namespace LHEH5

#include "PHASIC++/Main/Event_Reader.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace LHEH5 {

  class HDF5_Reader: public Event_Reader, public MPI_Object {
  private:

    LHEFile *p_file;
    size_t   m_ievt, m_ifile, m_trials, m_nstart, m_ncache;

    Vec4D_Vector m_ctmoms;

  public:

    HDF5_Reader(const Event_Reader_Key &key):
      Event_Reader(key), m_ievt(0), m_ifile(0), m_trials(0),
      m_nstart(0), m_ncache(0)
    {
      Data_Reader read(" ",";","#","=");
      m_ncache=read.GetValue<int>("HDF5_CACHE_SIZE",10000);
      p_file = OpenFile(m_files[m_ifile]);
      p_sub = new NLO_subevt();
      s_objects.push_back(this);
    }

    ~HDF5_Reader()
    {
      if (p_ampl) p_ampl->Delete();
      delete p_sub;
    }

    void MPISync()
    {
      delete p_file;
      if (m_nstart==0) ++m_ifile;
      p_file = OpenFile(m_files[m_ifile]);
    }

    LHEFile *OpenFile(const std::string &fname)
    {
      m_ievt=0;
      int size(mpi->MySize()), rank(mpi->MyRank());
      FileAccessProps fapl;
      fapl.add(MPIOFileAccess{MPI_COMM_WORLD, MPI_INFO_NULL});
      fapl.add(MPIOCollectiveMetadata{});
      File file(fname,File::ReadOnly,fapl);
      LHEFile *e(new LHEFile());
      e->ReadHeader(file);
      m_totalxs=e->TotalXS();
      long int nevts(file.getDataSet("events").getSpace().
		     getDimensions().front());
      if (mpi->Rank()==0) {
	msg_Info()<<METHOD<<"(): File '"<<fname
		  <<"' contains "<<nevts
		  <<" events ( read "<<m_nstart<<" ).\n";
      }
      mpi->Bcast(&nevts,1,MPI_LONG_INT,0);
      size_t iStart(rank*nevts/size);
      size_t iStop((rank+1)*nevts/size-1);
      if (rank==size-1) iStop=nevts-1;
      size_t nread(std::min(m_ncache,iStop-iStart-m_nstart+1));
      e->ReadEvents(file,iStart+m_nstart,nread);
      if (iStart+(m_nstart+=nread)>iStop) m_nstart=0;
      return e;
    }

    Cluster_Amplitude *ReadEvent()
    {
      DEBUG_FUNC("i="<<m_ievt<<"("<<p_file->NEvents()<<"),trial="<<m_trials);
      if (m_trials==0 && p_ampl!=NULL) {
	++m_trials;
	--m_ievt;
      }
      if (m_ievt>=p_file->NEvents()) {
	if ((m_ifile+1>=m_files.size() && m_nstart == 0) || Communicate()<0) {
	  msg_Info()<<METHOD<<"(): No more events in "<<m_files<<".\n";
	  rpa->gen.SetNumberOfEvents(rpa->gen.NumberOfGeneratedEvents());
	  return NULL;
	}
	m_ievt=0;
      }
      if (p_ampl==NULL) {
	Event e(p_file->GetEvent(m_ievt));
	msg_Debugging()<<e<<"\n";
	if (e.empty()) {
	  m_trials+=e.trials;
	  m_ievt++;
	  return NULL;
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
    
  };// end of class HDF5_Reader
  
}// end of namespace LHEH5

using namespace LHEH5;

DECLARE_GETTER(HDF5_Reader,"HDF5",Event_Reader,Event_Reader_Key);

Event_Reader *ATOOLS::Getter<Event_Reader,Event_Reader_Key,HDF5_Reader>::
operator()(const Event_Reader_Key &args) const
{
  return new HDF5_Reader(args);
}

void ATOOLS::Getter<Event_Reader,Event_Reader_Key,HDF5_Reader>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HDF5 reader (version 2)";
}

#endif
#endif
