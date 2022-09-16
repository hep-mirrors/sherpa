#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
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
    
    std::vector<size_t> start, trials;
    std::vector<int> pid, np, nplo, npnlo, version;
    std::vector<int> id, st, mo1, mo2, cl1, cl2;
    std::vector<double> unitwgt, xsec, px, py, pz, e, m;
    std::vector<double> mur, muf, muq, aqed, aqcd;
    std::vector<std::vector<double> > wgts;
    hid_t dspace;
    std::vector<std::string> wgtnames;
    size_t hasmwts;

    template <typename Type> void Read
    (Group &group,Type &values,const std::string &name)
    {
      DataSet ds(group.getDataSet(name));
      ds.read(values);
    }

    template <typename Type> void ReadSlice
    (Group &group,Type &values,const std::string &name,
     const std::vector<size_t> &offset,const std::vector<size_t> &size)
    {
      DataSet ds(group.getDataSet(name));
      ds.select(offset,size).read(values);
    }

    inline Particle GetParticle(size_t i) const
    {
      return Particle(id[i],st[i],mo1[i],mo2[i],
		      cl1[i],cl2[i],px[i],py[i],pz[i],e[i],m[i]);
    }

  public:

    inline double TotalXS() const
    { double xs(0.); for (auto c: xsec) xs+=c; return xs; }

    inline const std::vector<std::string> &
    WeightNames() const { return wgtnames; }

    inline size_t NProcesses() const { return xsec.size(); }
    inline ProcInfo GetProcInfo(const size_t i) const
    {
      return ProcInfo(pid[i]-1,nplo[pid[i]-1],npnlo[pid[i]-1],
		      unitwgt[pid[i]-1],xsec[pid[i]-1]);
    }

    inline size_t NEvents() const { return start.size(); }
    inline Event GetEvent(size_t i) const
    {
      Event e(GetProcInfo(i),trials[i],wgts[i],
	      mur[i],muf[i],muq[i],aqed[i],aqcd[i]);
      for (int n(0);n<np[i];++n)
	e.push_back(GetParticle(start[i]-start[0]+n));
      return e;
    }
 
    void ReadHeader(File &file)
    {
      Group init(file.getGroup("init"));
      Group procinfo(file.getGroup("procInfo"));
      hid_t dspace(H5Dget_space(file.getDataSet("procInfo/xSection").getId()));
      long int n(H5Sget_simple_extent_npoints(dspace));
      Read(procinfo,unitwgt,"unitWeight");
      Read(procinfo,xsec,"xSection");
      version={0,1,0};
      if (file.exist("/init/version")) Read(init,version,"version");
      if (version[0]>=0) {
	Read(procinfo,nplo,"npLO");
	Read(procinfo,npnlo,"npNLO");
	dspace=H5Dget_space(file.getDataSet("event/start").getId());
      }
      else {
	nplo.resize(n,-1);
	npnlo.resize(n,-1);
	dspace=H5Dget_space(file.getDataSet("index/start").getId());
      }
      if (file.exist("event/weight")&&version[0]>=1) {
	Group event(file.getGroup("event"));
	DataSet weights(event.getDataSet("weight"));
	auto attr_keys(weights.listAttributeNames());
	Attribute a(weights.getAttribute(attr_keys[0]));
	a.read(wgtnames);
	hasmwts=true;
      }
      else {
	wgtnames={"NOMINAL"};
	hasmwts=false;
      }
    }
    void ReadEvents(File &file,size_t first_event,size_t n_events)
    {
      Group g_event(file.getGroup("event"));
      std::vector<size_t> offset_e{first_event}, size_e{n_events};
      DataSet _weight(g_event.getDataSet("weight"));
      if (hasmwts) {
	std::vector<size_t> offsets{first_event,0};
	std::vector<size_t> counts{n_events,wgtnames.size()};
	_weight.select(offsets,counts).read(wgts);
      }
      else {
	std::vector<double> weights;
	_weight.select(offset_e,size_e).read(weights);
	wgts.resize(n_events,std::vector<double>(1));
	for (size_t i(0);i<n_events;++i) wgts[i][0]=weights[i];
      }
      ReadSlice(g_event,start,"start",offset_e,size_e);
      ReadSlice(g_event,np,"nparticles",offset_e,size_e);
      ReadSlice(g_event,pid,"pid",offset_e,size_e);
      ReadSlice(g_event,trials,"trials",offset_e,size_e);
      ReadSlice(g_event,mur,"rscale",offset_e,size_e);
      ReadSlice(g_event,muf,"fscale",offset_e,size_e);
      ReadSlice(g_event,muq,"scale",offset_e,size_e);
      ReadSlice(g_event,aqed,"aqed",offset_e,size_e);
      ReadSlice(g_event,aqcd,"aqcd",offset_e,size_e);
      Group g_particle(file.getGroup("particle"));
      size_t nps(start.back()-start.front()+np.back());
      std::vector<size_t> size_p{nps}, offset_p{start.front()};
      ReadSlice(g_particle,id,"id",offset_p,size_p);
      ReadSlice(g_particle,st,"status",offset_p,size_p);
      ReadSlice(g_particle,mo1,"mother1",offset_p,size_p);
      ReadSlice(g_particle,mo2,"mother2",offset_p,size_p);
      ReadSlice(g_particle,cl1,"color1",offset_p,size_p);
      ReadSlice(g_particle,cl2,"color2",offset_p,size_p);
      ReadSlice(g_particle,px,"px",offset_p,size_p);
      ReadSlice(g_particle,py,"py",offset_p,size_p);
      ReadSlice(g_particle,pz,"pz",offset_p,size_p);
      ReadSlice(g_particle,e,"e",offset_p,size_p);
      ReadSlice(g_particle,m,"m",offset_p,size_p);
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
      long int nevts;
      hid_t dspace(H5Dget_space(file.getDataSet("event/start").getId()));
      if (rank==0) {
	nevts=H5Sget_simple_extent_npoints(dspace);
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
  str<<"HDF5 reader";
}

#endif
#endif
