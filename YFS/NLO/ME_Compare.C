#include "YFS/NLO/ME_Compare.H"

#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"

#include <cmath>
#include <iomanip>

using namespace YFS;
using namespace ATOOLS;

ME_Compare::ME_Compare(bool check, const std::string &histdir,
                        const std::string &label, double tol)
  : m_check(check), m_tol(tol), m_histdir(histdir), m_label(label),
    m_npoints(0), m_nbad(0), m_nbad_is(0), m_nbad_fs(0), m_maxdev(0.),
    m_scatter_ok(0), m_scatter_max(200000)
{
  if(m_check){
    if (!ATOOLS::DirectoryExists("./"+m_histdir)) ATOOLS::MakeDir("./"+m_histdir);
    // per-photon scatter: one row per photon per checked point, so energy can
    // be plotted against the opening angle to the nearest charged particle.
    // Per-rank filename so MPI ranks do not clobber a shared file.
    int rank = 0;
#ifdef USING__MPI
    if (mpi->Size() > 1) rank = mpi->Rank();
#endif
    m_scatter.open((std::string("./")+m_histdir+"/photon_angle_energy_rank"
                    +std::to_string(rank)+".dat").c_str());
    m_scatter<<"# E_gamma  pT_gamma  theta_charged[rad]  nearest_kf  isIS  "
             <<"ratio  dev  bad\n";
    // log10 of photon energy/pT [GeV], of the opening angle [rad] to the
    // nearest charged leg, and of the angle [rad] to the nearest beam
    // direction, filled separately for points that agreed vs. disagreed
    // with the reference ME.
    m_histograms1d["photon_logE_ok"]       = new Histogram(0,-10.,3.,130);
    m_histograms1d["photon_logE_bad"]      = new Histogram(0,-10.,3.,130);
    m_histograms1d["photon_logPT_ok"]      = new Histogram(0,-10.,3.,130);
    m_histograms1d["photon_logPT_bad"]     = new Histogram(0,-10.,3.,130);
    m_histograms1d["photon_logTheta_ok"]   = new Histogram(0,-14.,1.,150);
    m_histograms1d["photon_logTheta_bad"]  = new Histogram(0,-14.,1.,150);
    m_histograms1d["photon_logThetaZ_ok"]  = new Histogram(0,-14.,1.,150);
    m_histograms1d["photon_logThetaZ_bad"] = new Histogram(0,-14.,1.,150);
    // log10 of the minimum photon-photon opening angle [rad], when there
    // are 2+ photons (e.g. RealReal) - the photon-photon collinear/soft
    // region has no analogue in a single-photon (Real) process.
    m_histograms1d["photon_logThetaGG_ok"]  = new Histogram(0,-14.,1.,150);
    m_histograms1d["photon_logThetaGG_bad"] = new Histogram(0,-14.,1.,150);
  }
}

ME_Compare::~ME_Compare()
{
  if(m_check && m_npoints>0){
    msg_Out()<<ATOOLS::om::bold<<ATOOLS::om::blue
             <<"###############################################"<<std::endl
             <<m_label<<" ME comparison summary:"<<ATOOLS::om::reset<<std::endl
             <<"  points checked    = "<<m_npoints<<std::endl
             <<"  points mismatched = "<<m_nbad<<" ("
             <<std::setprecision(3)<<100.*m_nbad/m_npoints<<"%)"<<std::endl
             <<"  max |1-ratio|     = "<<std::setprecision(6)<<m_maxdev<<std::endl;
    if(m_nbad_is+m_nbad_fs>0)
      msg_Out()<<"  mismatched, nearest charged leg: IS = "<<m_nbad_is
               <<", FS = "<<m_nbad_fs<<std::endl;
    if(m_nbad>0)
      msg_Out()<<ATOOLS::om::bold<<ATOOLS::om::red
               <<"WARNING: "<<m_nbad<<" / "<<m_npoints
               <<" points disagreed beyond tolerance!"<<ATOOLS::om::reset<<std::endl;
    else
      msg_Out()<<ATOOLS::om::bold<<ATOOLS::om::green
               <<"All points agreed within tolerance."<<ATOOLS::om::reset<<std::endl;
    msg_Out()<<ATOOLS::om::bold<<ATOOLS::om::blue
             <<"###############################################"
             <<ATOOLS::om::reset<<std::endl;
    for(auto hit: m_histograms1d){
      hit.second->MPISync();
      hit.second->Finalize();
      hit.second->Output(std::string("./")+m_histdir+"/"+hit.first+".dat");
      delete hit.second;
    }
  }
  if(m_scatter.is_open()) m_scatter.close();
}

void ME_Compare::CheckAgreement(const ATOOLS::Vec4D_Vector &p, double test,
                                 double reference,
                                 const ATOOLS::Flavour_Vector &flavs,
                                 size_t nin)
{
  double ratio = test/reference;
  double dev = std::abs(1.-ratio);
  m_npoints++;
  if(dev > m_maxdev) m_maxdev = dev;
  bool bad = dev > m_tol;
  if(bad) m_nbad++;

  double softe(-1.), softpt(-1.);
  int softi(-1);
  for(size_t i(0);i<flavs.size() && i<p.size();++i){
    if(flavs[i].IsPhoton() && (softe<0. || p[i][0]<softe)){
      softe = p[i][0];
      softpt = p[i].PPerp();
      softi = i;
    }
  }
  double mintheta(-1.);
  int mini(-1);
  if(softi>=0){
    for(size_t j(0);j<flavs.size() && j<p.size();++j){
      if((int)j==softi || flavs[j].IntCharge()==0) continue;
      double theta = p[softi].Theta(p[j]);
      if(mintheta<0. || theta<mintheta){
        mintheta = theta;
        mini = j;
      }
    }
    if(bad && mini>=0){
      if((size_t)mini<nin) m_nbad_is++;
      else m_nbad_fs++;
    }
  }
  double thetaz(-1.);
  if(softi>=0){
    double t = p[softi].Theta();
    thetaz = std::min(t, M_PI-t);
  }
  double ggtheta(-1.);
  {
    std::vector<size_t> photons;
    for(size_t i(0);i<flavs.size() && i<p.size();++i)
      if(flavs[i].IsPhoton()) photons.push_back(i);
    for(size_t a(0);a<photons.size();++a)
      for(size_t b(a+1);b<photons.size();++b){
        double theta = p[photons[a]].Theta(p[photons[b]]);
        if(ggtheta<0. || theta<ggtheta) ggtheta = theta;
      }
  }
  if(softi>=0){
    const std::string suffix = bad?"_bad":"_ok";
    m_histograms1d["photon_logE"+suffix]->Insert(log10(Max(softe,1e-30)));
    m_histograms1d["photon_logPT"+suffix]->Insert(log10(Max(softpt,1e-30)));
    if(mini>=0)
      m_histograms1d["photon_logTheta"+suffix]->Insert(log10(Max(mintheta,1e-30)));
    m_histograms1d["photon_logThetaZ"+suffix]->Insert(log10(Max(thetaz,1e-30)));
    if(ggtheta>=0.)
      m_histograms1d["photon_logThetaGG"+suffix]->Insert(log10(Max(ggtheta,1e-30)));
  }

  // per-photon scatter file: for EVERY photon (not just the softest) record its
  // energy and the angle to its nearest charged particle, so energetic photons
  // that are collinear-enhanced show up as high-E / small-angle points. Always
  // write points that disagree (bad); cap the agreeing ones.
  if (m_scatter.is_open() && (bad || m_scatter_ok < m_scatter_max)) {
    if (!bad) ++m_scatter_ok;
    for(size_t i(0);i<flavs.size() && i<p.size();++i){
      if(!flavs[i].IsPhoton()) continue;
      double th(-1.); int ni(-1);
      for(size_t j(0);j<flavs.size() && j<p.size();++j){
        if(j==i || flavs[j].IntCharge()==0) continue;
        double t = p[i].Theta(p[j]);
        if(th<0. || t<th){ th = t; ni = j; }
      }
      m_scatter<<std::setprecision(10)
               <<p[i][0]<<" "<<p[i].PPerp()<<" "<<th<<" "
               <<(ni>=0?(flavs[ni].IsAnti()?-(long)flavs[ni].Kfcode():(long)flavs[ni].Kfcode()):0)<<" "
               <<(ni>=0 && (size_t)ni<nin?1:0)<<" "
               <<ratio<<" "<<dev<<" "<<(bad?1:0)<<"\n";
    }
  }

  if(bad){
    msg_Out()<<ATOOLS::om::bold<<ATOOLS::om::red<<"WARNING: "<<ATOOLS::om::reset
             <<ATOOLS::om::red<<m_label<<" ME mismatch: ratio = "
             <<std::setprecision(10)<<ratio<<", |1-ratio| = "<<dev;
    if(softe>=0.) msg_Out()<<", softest photon: E = "<<softe<<", pT = "<<softpt;
    if(mini>=0) msg_Out()<<", min angle to charged leg ["<<flavs[mini].IDName()
                          <<((size_t)mini<nin?" (IS)":" (FS)")
                          <<"] = "<<mintheta<<" rad";
    if(softi>=0) msg_Out()<<", angle to beam axis = "<<thetaz<<" rad";
    if(ggtheta>=0.) msg_Out()<<", min photon-photon angle = "<<ggtheta<<" rad";
    msg_Out()<<ATOOLS::om::reset<<std::endl;
  }
  else {
    msg_Debugging()<<ATOOLS::om::green<<m_label<<" ME OK: ratio = "
                    <<std::setprecision(10)<<ratio<<ATOOLS::om::reset<<std::endl;
  }
}
