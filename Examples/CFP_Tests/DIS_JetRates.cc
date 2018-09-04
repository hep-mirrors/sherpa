// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {

  class DIS_JetRates : public Analysis {
  private:

    FourMomentum _p_p;

    std::vector<Histo1DPtr> _h_log10_d;

  public:

    DIS_JetRates(): Analysis("DIS_JetRates") {}

    void init()
    {
      const DISKinematics &dk = addProjection(DISKinematics(), "Kinematics");
      addProjection(DISFinalState(dk, DISFinalState::BREIT), "FS");
      for (size_t i=0;i<4;++i) {
	string dname="log10_d_"+to_str(i?i+1:0)+to_str(i+2);
	_h_log10_d.push_back(bookHisto1D(dname,100,0.0,log10(0.5*sqrtS()/GeV)));
      }
    }

    void analyze(const Event& event)
    {
      const double weight=event.weight();
      const DISKinematics &dk = applyProjection<DISKinematics>(event, "Kinematics");
      const DISFinalState &fs = applyProjection<DISFinalState>(event, "FS");
      _p_p=dk.beamHadron().momentum();
      std::vector<FourMomentum> jets(fs.size());
      for (size_t i(0);i<fs.size();++i)
	jets[i]=fs.particles()[i].momentum();
      std::vector<double> dij2 = Cluster(jets);
      for (size_t i=0;i<min(_h_log10_d.size(),dij2.size()-1);++i)
	if (dij2[dij2.size()-2-i])
	  _h_log10_d[i]->fill(0.5*log10(dij2[dij2.size()-2-i]), weight);
    }

    void finalize()
    {
      const double xsec_unitw=crossSection()/picobarn/sumOfWeights();
      for (size_t i=0;i<_h_log10_d.size();++i) scale(_h_log10_d[i],xsec_unitw);
    }

    double R2(const FourMomentum &p, const FourMomentum &q) const
    {
      double pq(p.x()*q.x()+p.y()*q.y()+p.z()*q.z());
      double p2(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());
      double q2(q.x()*q.x()+q.y()*q.y()+q.z()*q.z());
      return 2.0*(1.0-min(max(pq/sqrt(p2*q2),-1.0),1.0));
    }

    double Q2i(const FourMomentum &p) const
    {
      return sqr(p.E())*R2(p,_p_p);
    }

    double Q2ij(const FourMomentum &p,const FourMomentum &q) const
    {
      return sqr(min(p.E(),q.E()))*R2(p,q);
    }

    std::vector<double> Cluster(std::vector<FourMomentum> &p)
    {
      std::vector<double> kt2;
      if (p.size()==1) kt2.push_back(Q2i(p[0]));
      if (p.size()<=1) return kt2;
      int ii=0, jj=0, n=p.size();
      std::vector<int> imap(p.size());
      for (size_t i(0);i<imap.size();++i) imap[i]=i;
      std::vector<std::vector<double> > kt2ij(n,std::vector<double>(n));
      double dmin=std::numeric_limits<double>::max();
      for (int i=0;i<n;++i) {
	double di=kt2ij[i][i]=Q2i(p[i]);
	if (di<dmin) { dmin=di; ii=i; jj=i; }
	for (int j=0;j<i;++j) {
	  double dij=kt2ij[i][j]=Q2ij(p[i],p[j]);
	  if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	}
      }
      while (n>0) {
	MSG_DEBUG("Q_{"<<n<<"->"<<n-1<<"} = "<<sqrt(dmin)
		  <<" <- "<<p[imap[jj]]<<" "<<p[imap[ii]]);
	if (ii!=jj) p[imap[jj]]+=p[imap[ii]];
	kt2.push_back(dmin); --n;
	for (int i=ii;i<n;++i) imap[i]=imap[i+1];
	int jjx=imap[jj];
	kt2ij[jjx][jjx]=Q2i(p[jjx]);
	for (int j=0;j<jj;++j) kt2ij[jjx][imap[j]]=Q2ij(p[jjx],p[imap[j]]);
	for (int i=jj+1;i<n;++i) kt2ij[imap[i]][jjx]=Q2ij(p[jjx],p[imap[i]]);
	ii=jj=0; dmin=kt2ij[imap[0]][imap[0]];
	for (int i=0;i<n;++i) {
	  int ix=imap[i]; double di=kt2ij[ix][ix];
	  if (di<dmin) { dmin=di; ii=jj=i; }
	  for (int j=0;j<i;++j) {
	    int jx=imap[j]; double dij=kt2ij[ix][jx];
	    if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	  }
	}
      }
      return kt2;
    }

  };

  DECLARE_RIVET_PLUGIN(DIS_JetRates);

}
