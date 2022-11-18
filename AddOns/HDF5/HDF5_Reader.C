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
    inline Event(const ProcInfo &_pinfo,
		 size_t _trials,std::vector<double> _wgts,
		 double _mur,double _muf,double _muq,
		 double _aqed,double _aqcd):
      pinfo(_pinfo), trials(_trials), wgts(_wgts),
      mur(_mur), muf(_muf), muq(_muq), aqed(_aqed), aqcd(_aqcd) {}
  };// end of struct Event

  std::ostream &operator<<(std::ostream &s,const Event &e)
  { s<<"Event "<<e.pinfo<<" {\n"
     <<"  trials="<<e.trials<<",weights=("<<e.wgts[0];
    for (size_t i(1);i<e.wgts.size();++i) s<<","<<e.wgts[i];
    s<<")\n  mur="<<e.mur<<", muf="<<e.muf<<", muq="<<e.muq
     <<",aqed="<<e.aqed<<",aqcd="<<e.aqcd<<"\n";
    for (size_t i(0);i<e.size();++i) s<<"  "<<e[i]<<"\n";
    return s<<"}"; }

  class LHEFile {
  private:
    
    std::vector<int> version;
    std::vector<std::vector<double> > evts, parts, pinfo;
    std::vector<std::string> wgtnames;

    inline Particle GetParticle(size_t i) const
    {
      return Particle(parts[i][0],parts[i][1],parts[i][2],parts[i][3],
		      parts[i][4],parts[i][5],parts[i][6],parts[i][7],
		      parts[i][8],parts[i][9],parts[i][10]);
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
      Event e(GetProcInfo(evts[i][0]-1),evts[i][3],wgts,
	      evts[i][6],evts[i][5],evts[i][4],
	      evts[i][7],evts[i][8]);
      for (int n(0);n<evts[i][1];++n)
	e.push_back(GetParticle(evts[i][2]-evts[0][2]+n));
      return e;
    }
 
    void ReadHeader(File &file)
    {
      file.getDataSet("version").read(version);
      file.getDataSet("procInfo").read(pinfo);
      DataSet init(file.getDataSet("init"));
      DataSet events(file.getDataSet("events"));
      auto attr_keys(events.listAttributeNames());
      Attribute a(events.getAttribute(attr_keys[0]));
      a.read(wgtnames);
      for (int i(0);i<9;++i) wgtnames.erase(wgtnames.begin());
    }
    void ReadEvents(File &file,size_t first_event,size_t n_events)
    {
      DataSet events(file.getDataSet("events"));
      std::vector<size_t> eoffsets{first_event,0};
      std::vector<size_t> ecounts{n_events,9+wgtnames.size()};
      events.select(eoffsets,ecounts).read(evts);
      DataSet particles(file.getDataSet("particles"));
      std::vector<size_t> poffsets{(size_t)evts.front()[2],0};
      size_t nps(evts.back()[2]-evts.front()[2]+evts.back()[1]);
      std::vector<size_t> pcounts{nps,13};
      particles.select(poffsets,pcounts).read(parts);
    }
    
  };// end of struct LHEFile

}// end of namespace LHEH5

#include "PHASIC++/Main/Event_Reader.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace LHEH5 {

  class HDF5_Reader: public Event_Reader {
  private:

    LHEFile *p_file;
    size_t   m_ievt, m_ifile, m_trials;

    Cluster_Amplitude *p_ampl;
    
  public:

    HDF5_Reader(const Event_Reader_Key &key):
      Event_Reader(key), m_ievt(0), m_ifile(0), p_ampl(NULL)
    {
      p_file = OpenFile(m_files[m_ifile]);
    }

    ~HDF5_Reader()
    {
      if (p_ampl) p_ampl->Delete();
    }

    LHEFile *OpenFile(const std::string &fname)
    {
      m_ievt=0;
      int size(mpi->Size()), rank(mpi->Rank());
      File file(fname,File::ReadOnly,
		MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));
      LHEFile *e(new LHEFile());
      e->ReadHeader(file);
      m_totalxs=e->TotalXS();
      long int nevts(file.getDataSet("events").getSpace().
		     getDimensions().front());
      if (rank==0) {
	msg_Info()<<METHOD<<"(): File '"<<fname
		  <<"' contains "<<nevts<<" events.\n";
      }
      mpi->Bcast(&nevts,1,MPI_LONG_INT,0);
      size_t iStart(rank*nevts/size);
      size_t iStop((rank+1)*nevts/size-1);
      if (rank==size-1) iStop=nevts-1;
      e->ReadEvents(file,iStart,iStop-iStart+1);
      return e;
    }

    Cluster_Amplitude *ReadEvent()
    {
      DEBUG_FUNC("");
      if (m_ievt>=p_file->NEvents()) {
	delete p_file;
	p_file=NULL;
	++m_ifile;
	if (m_ifile>=m_files.size()) {
	  msg_Info()<<METHOD<<"(): No more events in "<<m_files<<".\n";
	  rpa->gen.SetNumberOfEvents(rpa->gen.NumberOfGeneratedEvents());
	  return NULL;
	}
	p_file = OpenFile(m_files[m_ifile]);
	m_ievt=0;
      }
      if (p_ampl==NULL) {
	Event e(p_file->GetEvent(m_ievt++));
	if (e[0].pz<0 && e[1].pz>0) std::swap<Particle>(e[0],e[1]);
	msg_Debugging()<<e<<"\n";
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
	m_trials=e.trials;
      }
      if (m_trials==1) {
	Cluster_Amplitude *ampl(p_ampl);
	p_ampl=NULL;
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
