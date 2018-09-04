// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

  class LL_JetRates : public Analysis {
  private:

    std::vector<Histo1DPtr> _h_log10_d;

    Histo1DPtr _h_cosphi1, _h_cosphi2, _h_a34;

    double m_a34;

  public:

    LL_JetRates(): Analysis("LL_JetRates") {}

    void init()
    {
      addProjection(FinalState(), "FS");
      for (size_t i=0;i<4;++i) {
	string dname="log10_y_"+to_str(i+2)+to_str(i+3);
	_h_log10_d.push_back(bookHisto1D(dname,100,-4.3,-0.3));
      }
      _h_cosphi1 = bookHisto1D("cos_phi1",100,-1.0,1.0);
      _h_cosphi2 = bookHisto1D("cos_phi2",100,-1.0,1.0);
      _h_a34 = bookHisto1D("alpha_34",100,-1.0,1.0);
    }

    void analyze(const Event& event)
    {
      const double weight=event.weight();
      const FinalState &fs = applyProjection<FinalState>(event, "FS");
      std::vector<FourMomentum> jets(fs.size());
      for (size_t i(0);i<fs.size();++i) {
	MSG_DEBUG(i<<": "<<fs.particles()[i].pdgId()<<" "<<fs.particles()[i].momentum());
	jets[i]=fs.particles()[i].momentum();
      }
      if (jets.size()==4) {
	double si1=jets[0]*jets[2], sj1=jets[1]*jets[2], s12=jets[2]*jets[3];
	double si2=jets[0]*jets[3], sj2=jets[1]*jets[3], sij=jets[0]*jets[1];
	double cp1=(si1*sj2+si2*sj1-s12*sij)/sqrt(4.*si1*sj1*si2*sj2);
	double cp2=(si1*sj2-si2*sj1)/sqrt(s12*sij*(si1+si2)*(sj1+sj2));
	double z1z2=(si1+sj1+s12)*(si2+sj2+s12)/sqr(si1+si2+sj1+sj2+sij+s12);
	_h_cosphi1->fill(cp1, weight);
	_h_cosphi2->fill(cp2/sqrt(4.*z1z2), weight);
      }
      std::vector<double> dij2 = Cluster(jets, weight);
      for (size_t i=0;i<min(_h_log10_d.size(),dij2.size());++i)
	if (dij2[dij2.size()-1-i])
	  _h_log10_d[i]->fill(log10(dij2[dij2.size()-1-i]), weight);
    }

    void finalize()
    {
      const double xsec_unitw=crossSection()/picobarn/sumOfWeights();
      scale(_h_cosphi1,xsec_unitw);
      scale(_h_cosphi2,xsec_unitw);
      for (size_t i=0;i<_h_log10_d.size();++i) scale(_h_log10_d[i],xsec_unitw);
      scale(_h_a34,xsec_unitw);
    }

    double Yij(const FourMomentum &p, const FourMomentum &q) const
    {
      double pq(p.x()*q.x()+p.y()*q.y()+p.z()*q.z());
      double p2(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());
      double q2(q.x()*q.x()+q.y()*q.y()+q.z()*q.z());
      return 2.0*sqr(min(p.E(),q.E()))*(1.0-min(max(pq/sqrt(p2*q2),-1.0),1.0))/sqr(sqrtS());
    }

    std::vector<double> Cluster(std::vector<FourMomentum> &p, const double &weight)
    {
      std::vector<double> kt2;
      if (p.size()<=1) return kt2;
      int ii=0, jj=0, n=p.size();
      std::vector<int> imap(p.size());
      for (size_t i(0);i<imap.size();++i) imap[i]=i;
      std::vector<std::vector<double> > kt2ij(n,std::vector<double>(n));
      double dmin=std::numeric_limits<double>::max();
      for (int i=0;i<n;++i)
	for (int j=0;j<i;++j) {
	  double dij=kt2ij[i][j]=Yij(p[i],p[j]);
	  if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	}
      for (--n;n>1;--n) {
	MSG_DEBUG("y_{"<<n+1<<"->"<<n<<"} = "<<dmin
		  <<" <- "<<p[imap[jj]]<<" "<<p[imap[ii]]);
	if (n==3) {
	  FourMomentum pp(p[imap[ii]]), qq(p[imap[jj]]);
	  double pq(pp.x()*qq.x()+pp.y()*qq.y()+pp.z()*qq.z());
	  double p2(pp.x()*pp.x()+pp.y()*pp.y()+pp.z()*pp.z());
	  double q2(qq.x()*qq.x()+qq.y()*qq.y()+qq.z()*qq.z());
	  m_a34=pq/sqrt(p2*q2);
	  MSG_DEBUG("\\alpha_{34} = "<<m_a34);
	  _h_a34->fill(m_a34, weight);
	}
	kt2.push_back(dmin);
	int jjx=imap[jj]; p[jjx]+=p[imap[ii]];
	for (int i=ii;i<n;++i) imap[i]=imap[i+1];
	for (int j=0;j<jj;++j) kt2ij[jjx][imap[j]]=Yij(p[jjx],p[imap[j]]);
	for (int i=jj+1;i<n;++i) kt2ij[imap[i]][jjx]=Yij(p[jjx],p[imap[i]]);
	ii=jj=0; dmin=std::numeric_limits<double>::max();
	for (int i=0;i<n;++i)
	  for (int j=0;j<i;++j) {
	    double dij=kt2ij[imap[i]][imap[j]];
	    if (dij<dmin) { dmin=dij; ii=i; jj=j; }
	  }
      }
      return kt2;
    }

  };

  DECLARE_RIVET_PLUGIN(LL_JetRates);

}
