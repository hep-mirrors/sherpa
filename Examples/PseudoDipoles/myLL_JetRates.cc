// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

  class LL_JetRates : public Analysis {
  private:

    std::vector<Histo1DPtr> _h_log10_d;
    std::vector<Histo1DPtr> _h_log10_d_b_or_g;
    std::vector<Histo1DPtr> _h_log10_d_bb_or_g;

  public:

    LL_JetRates(): Analysis("myLL_JetRates") {}

    void init()
    {
      addProjection(FinalState(), "FS");
      for (size_t i=0;i<1;++i) {
	string dname="log10_y_"+to_str(i+2)+to_str(i+3);
	string dname_b="log10_y_"+to_str(i+2)+to_str(i+3)+"_b_or_g";
	string dname_bb="log10_y_"+to_str(i+2)+to_str(i+3)+"_bb_or_g";
	_h_log10_d.push_back(bookHisto1D(dname,100,-8.0,-0.3));
	_h_log10_d_b_or_g.push_back(bookHisto1D(dname_b,100,-8.0,-0.3));
	_h_log10_d_bb_or_g.push_back(bookHisto1D(dname_bb,100,-8.0,-0.3));
      }
    }

    void analyze(const Event& event)
    {
      const double weight=event.weight();
      const FinalState &fs = applyProjection<FinalState>(event, "FS");
      Particles partons;
      for (size_t i(0);i<fs.size();++i) {
        if(fs.particles()[i].pdgId()!=24 && fs.particles()[i].pdgId()!=-24)
          partons.push_back(fs.particles()[i]);
      }
      Particles partons_b_or_g;
      for (size_t i(0);i<fs.size();++i) {
        if(fs.particles()[i].pdgId()!=24 && fs.particles()[i].pdgId()!=-24
        && fs.particles()[i].pdgId()!=-5)
          partons_b_or_g.push_back(fs.particles()[i]);
      }
      Particles partons_bb_or_g;
      for (size_t i(0);i<fs.size();++i) {
        if(fs.particles()[i].pdgId()!=24 && fs.particles()[i].pdgId()!=-24
        && fs.particles()[i].pdgId()!=5)
          partons_bb_or_g.push_back(fs.particles()[i]);
      }

//if(partons.size()>3)
//  std::cout << "partons = " << partons << "\n";
//  std::cout << "partons_b_or_g = " << partons_b_or_g << "\n";
//  std::cout << "partons_bb_or_g = " << partons_bb_or_g << "\n";

      std::vector<FourMomentum> jets(partons.size());
      for (size_t i(0);i<partons.size();++i) {
        MSG_DEBUG(i<<": "<<partons[i].pdgId()<<" "<<partons[i].momentum());
        jets[i]=partons[i].momentum();
      }
      std::vector<FourMomentum> jets_b_or_g(partons_b_or_g.size());
      for (size_t i(0);i<partons_b_or_g.size();++i) {
        MSG_DEBUG(i<<": "<<partons_b_or_g[i].pdgId()<<" "<<partons_b_or_g[i].momentum());
        jets_b_or_g[i]=partons_b_or_g[i].momentum();
      }
      std::vector<FourMomentum> jets_bb_or_g(partons_bb_or_g.size());
      for (size_t i(0);i<partons_bb_or_g.size();++i) {
        MSG_DEBUG(i<<": "<<partons_bb_or_g[i].pdgId()<<" "<<partons_bb_or_g[i].momentum());
        jets_bb_or_g[i]=partons_bb_or_g[i].momentum();
      }
//std::cout << "jets = " << jets << "\n";
//std::cout << "jets_b_or_g = " << jets_b_or_g << "\n";
//std::cout << "jets_bb_or_g = " << jets_bb_or_g << "\n";

      std::vector<double> dij2 = Cluster(jets, weight);
//std::cout << "dij2 = " << dij2 << "\n";
      if (dij2.size()>0)
        _h_log10_d[0]->fill(log10(dij2[0]), weight);

      std::vector<double> dij2_b_or_g = littleCluster(jets_b_or_g, weight);
      std::vector<double> dij2_bb_or_g = littleCluster(jets_bb_or_g, weight);
//std::cout << "dij2_b_or_g = " << dij2_b_or_g << "\n";
//std::cout << "dij2_bb_or_g = " << dij2_bb_or_g << "\n";
      if(dij2_b_or_g.size()>0 && dij2_bb_or_g.size()>0) {
        if (dij2_b_or_g[0] < dij2_bb_or_g[0])
          _h_log10_d_b_or_g[0]->fill(log10(dij2_b_or_g[0]), weight);
        else if (dij2_bb_or_g[0] < dij2_b_or_g[0])
          _h_log10_d_bb_or_g[0]->fill(log10(dij2_bb_or_g[0]), weight);
      }
    }

    void finalize()
    {
      const double xsec_unitw=crossSection()/picobarn/sumOfWeights();
      for (size_t i=0;i<_h_log10_d.size();++i) scale(_h_log10_d[i],xsec_unitw);
      for (size_t i=0;i<_h_log10_d_b_or_g.size();++i) scale(_h_log10_d_b_or_g[i],xsec_unitw);
      for (size_t i=0;i<_h_log10_d_bb_or_g.size();++i) scale(_h_log10_d_bb_or_g[i],xsec_unitw);
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
//          std::cout << "dij = " << dij << ", dmin = " << dmin << "\n";
        }

      for (--n;n>1;--n) {
//        std::cout << "y_{"<<n+1<<"->"<<n<<"} = "<<dmin
//                  <<" <- "<<p[imap[jj]]<<" "<<p[imap[ii]] << "\n";
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
    std::vector<double> littleCluster(std::vector<FourMomentum> &p, const double &weight)
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
//          std::cout << "dij = " << dij << ", dmin = " << dmin << "\n";
        }

      for (--n;n>0;--n) {
//        std::cout << "y_{"<<n+1<<"->"<<n<<"} = "<<dmin
//                  <<" <- "<<p[imap[jj]]<<" "<<p[imap[ii]] << "\n";
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
