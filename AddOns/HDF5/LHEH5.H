#ifndef AddOns__HDF5__LHEH5_H
#define AddOns__HDF5__LHEH5_H

// In-memory representation of the LHEH5 parton-level event format, shared
// by every reader of that format (HDF5_Reader, which reads it from disk,
// and Pepper_Reader, which reads it from RAM-backed databases filled by
// the Pepper library).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
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

  // The column names of an LHEH5 dataset. They are stored in an attribute
  // called "properties" since LHEH5 version 2.1.0, while earlier versions
  // named some of these attributes after their datasets instead.
  inline std::vector<std::string> ColumnNames(const HighFive::DataSet &ds)
  {
    const std::vector<std::string> keys(ds.listAttributeNames());
    if (keys.empty())
      throw std::runtime_error("LHEH5 dataset \""+ds.getPath()+
			       "\" does not name its columns.");
    std::string key(keys.front());
    for (const auto& k : keys) if (k=="properties") key=k;
    std::vector<std::string> names;
    ds.getAttribute(key).read(names);
    return names;
  }

  // The indices of the columns `names` among the columns `all` of the
  // dataset `ds`.
  inline std::vector<size_t>
  ColumnIndices(const std::vector<std::string> &all,
		const std::vector<std::string> &names,
		const std::string &ds)
  {
    std::vector<size_t> cols;
    for (const auto& name : names) {
      const auto it(std::find(all.begin(),all.end(),name));
      if (it==all.end())
	throw std::runtime_error("LHEH5 dataset \""+ds+
				 "\" has no column \""+name+"\".");
      cols.push_back(it-all.begin());
    }
    return cols;
  }

  // Read the rows [first, first + n) of the dataset `ds`, keeping only the
  // columns `cols`, in that order.
  inline std::vector<std::vector<double> >
  ReadColumns(const HighFive::DataSet &ds,size_t first,size_t n,
	      const std::vector<size_t> &cols,
	      const HighFive::DataTransferProps &xfer_props)
  {
    const size_t ncols(ds.getSpace().getDimensions()[1]);
    std::vector<std::vector<double> > rows(n,std::vector<double>(ncols));
    ds.select({first,0},{n,ncols}).read(rows,xfer_props);
    for (auto& row : rows) {
      std::vector<double> sel(cols.size());
      for (size_t i(0);i<cols.size();++i) sel[i]=row[cols[i]];
      row.swap(sel);
    }
    return rows;
  }

  // The columns are located by name, and stored in the order given here,
  // such that the reader does not depend on the order they are written in.
  inline const std::vector<std::string> &ProcInfoColumns()
  {
    static const std::vector<std::string> names
      {"procId","npLO","npNLO","xSection","error","unitWeight"};
    return names;
  }
  inline const std::vector<std::string> &EventColumns()
  {
    static const std::vector<std::string> names
      {"pid","nparticles","start","trials",
       "scale","fscale","rscale","aqed","aqcd"};
    return names;
  }
  inline const std::vector<std::string> &ParticleColumns()
  {
    static const std::vector<std::string> names
      {"id","status","mother1","mother2","color1","color2",
       "px","py","pz","e","m"};
    return names;
  }
  inline const std::vector<std::string> &CTEventColumns()
  {
    static const std::vector<std::string> names
      {"ijt","kt","i","j","k","z1","z2","bbpsw","tlpsw"};
    return names;
  }
  inline const std::vector<std::string> &CTParticleColumns()
  {
    static const std::vector<std::string> names {"px","py","pz","e"};
    return names;
  }

  class LHEFile {
  private:
    
    std::vector<int> version;
    std::vector<std::vector<double> > evts, parts, pinfo;
    std::vector<std::vector<double> > ctevts, ctparts;
    // The weight names, starting with the nominal one. In the event rows,
    // the weights follow the columns listed in EventColumns().
    std::vector<std::string> wgtnames;
    std::vector<size_t> ecols;
    std::string psgen;

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

    // Empty if the file does not specify its phase-space generator, which
    // files written before LHEH5 version 2.5.0 never do.
    inline const std::string &PhaseSpaceGenerator() const { return psgen; }

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
      HighFive::DataSet procinfo(file.getDataSet("procInfo"));
      pinfo=ReadColumns
	(procinfo,0,procinfo.getSpace().getDimensions()[0],
	 ColumnIndices(ColumnNames(procinfo),ProcInfoColumns(),"procInfo"),
	 xfer_props);
      HighFive::DataSet init(file.getDataSet("init"));
      if (init.hasAttribute("phasespaceGenerator"))
        init.getAttribute("phasespaceGenerator").read(psgen);
      HighFive::DataSet events(file.getDataSet("events"));
      const std::vector<std::string> enames(ColumnNames(events));
      ecols=ColumnIndices(enames,EventColumns(),"events");
      // All other columns are weights. The nominal one is called "NOMINAL",
      // and is assumed to be the first one if no column has that name.
      wgtnames.clear();
      for (size_t i(0);i<enames.size();++i)
	if (std::find(ecols.begin(),ecols.end(),i)==ecols.end())
	  wgtnames.push_back(enames[i]);
      if (wgtnames.empty())
	throw std::runtime_error("LHEH5 dataset \"events\" has no weights.");
      const auto nominal(std::find(wgtnames.begin(),wgtnames.end(),"NOMINAL"));
      if (nominal!=wgtnames.end())
	std::rotate(wgtnames.begin(),nominal,nominal+1);
      const std::vector<size_t> wcols
	(ColumnIndices(enames,wgtnames,"events"));
      ecols.insert(ecols.end(),wcols.begin(),wcols.end());
    }
    void ReadEvents(HighFive::File &file,
                    size_t first_event,size_t n_events,
                    const HighFive::DataTransferProps &xfer_props
                    = HighFive::DataTransferProps{})
    {
      HighFive::DataSet events(file.getDataSet("events"));
      evts=ReadColumns(events,first_event,n_events,ecols,xfer_props);
      HighFive::DataSet particles(file.getDataSet("particles"));
      const size_t pstart(evts.front()[2]);
      size_t nmax(0);
      for (size_t i(0);i<pinfo.size();++i)
	nmax=std::max((size_t)std::max(pinfo[i][1],pinfo[i][2]+1),nmax);
      parts=ReadColumns
	(particles,pstart,n_events*nmax,
	 ColumnIndices(ColumnNames(particles),ParticleColumns(),"particles"),
	 xfer_props);
      if (file.exist("ctevents")) {
	HighFive::DataSet ctevents(file.getDataSet("ctevents"));
	ctevts=ReadColumns
	  (ctevents,first_event,n_events,
	   ColumnIndices(ColumnNames(ctevents),CTEventColumns(),"ctevents"),
	   xfer_props);
	HighFive::DataSet ctparticles(file.getDataSet("ctparticles"));
	ctparts=ReadColumns
	  (ctparticles,pstart,n_events*nmax,
	   ColumnIndices(ColumnNames(ctparticles),CTParticleColumns(),
			 "ctparticles"),
	   xfer_props);
      }
    }
    
  };// end of struct LHEFile

}// end of namespace LHEH5

#endif
