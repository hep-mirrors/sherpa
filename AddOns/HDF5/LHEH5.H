#ifndef AddOns__HDF5__LHEH5_H
#define AddOns__HDF5__LHEH5_H

// In-memory representation of the LHEH5 parton-level event format, shared
// by every reader of that format (HDF5_Reader, which reads it from disk,
// and Pepper_Reader, which reads it from RAM-backed databases filled by
// the Pepper library).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>

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
 
    void ReadHeader(HighFive::File &file,
                    const HighFive::DataTransferProps &xfer_props
                    = HighFive::DataTransferProps{})
    {
      file.getDataSet("version").read(version,xfer_props);
      file.getDataSet("procInfo").read(pinfo,xfer_props);
      HighFive::DataSet events(file.getDataSet("events"));
      auto attr_keys(events.listAttributeNames());
      HighFive::Attribute a(events.getAttribute(attr_keys[0]));
      a.read(wgtnames);
      for (int i(0);i<9;++i) wgtnames.erase(wgtnames.begin());
    }
    void ReadEvents(HighFive::File &file,
                    size_t first_event,size_t n_events,
                    const HighFive::DataTransferProps &xfer_props
                    = HighFive::DataTransferProps{})
    {
      HighFive::DataSet events(file.getDataSet("events"));
      std::vector<size_t> eoffsets{first_event,0};
      std::vector<size_t> ecounts{n_events,9+wgtnames.size()};
      evts.resize(n_events,std::vector<double>(9+wgtnames.size()));
      events.select(eoffsets,ecounts).read(evts,xfer_props);
      HighFive::DataSet particles(file.getDataSet("particles"));
      std::vector<size_t> poffsets{(size_t)evts.front()[2],0};
      size_t nmax(0);
      for (size_t i(0);i<pinfo.size();++i)
	nmax=std::max((size_t)std::max(pinfo[i][1],pinfo[i][2]+1),nmax);
      std::vector<size_t> pcounts{n_events*nmax,13};
      parts.resize(n_events*nmax,std::vector<double>(13));
      particles.select(poffsets,pcounts).read(parts,xfer_props);
      if (file.exist("ctevents")) {
	HighFive::DataSet events(file.getDataSet("ctevents"));
	std::vector<size_t> eoffsets{first_event,0};
	std::vector<size_t> ecounts{n_events,9};
	ctevts.resize(n_events,std::vector<double>(9));
	events.select(eoffsets,ecounts).read(ctevts,xfer_props);
	HighFive::DataSet particles(file.getDataSet("ctparticles"));
	std::vector<size_t> poffsets{(size_t)evts.front()[2],0};
	std::vector<size_t> pcounts{n_events*nmax,4};
	ctparts.resize(n_events*nmax,std::vector<double>(4));
	particles.select(poffsets,pcounts).read(ctparts,xfer_props);
      }
    }
    
  };// end of struct LHEFile

}// end of namespace LHEH5

#endif
