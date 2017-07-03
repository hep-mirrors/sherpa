#include "ATOOLS/Phys/Particle_Dresser.H"

#include <limits>
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ATOOLS;

Particle_Dresser::Particle_Dresser(const Flavour * fl,
                                 const size_t& nin, const size_t& nout,
                                 std::string algo, double dR, double exp) :
  m_on(true), p_fl(fl), m_nin(nin), m_nout(nout), m_n(m_nin+m_nout),
  m_algo(0), m_exp(exp)
{
  DEBUG_FUNC(m_nin<<" -> "<<m_nout<<", algo="<<algo<<", dR="<<dR);
  if      (algo=="Cone")          m_algo=0;
  else if (algo=="Recombination") m_algo=1;
  else                            m_algo=100;

  for (size_t i(m_nin);i<m_n;++i) {
    if (p_fl[i].IsPhoton()) m_photons.push_back(i);
    if (p_fl[i].Charge())   m_charges.push_back(i);
  }
  m_dR2.resize(m_charges.size(),dR*dR);
  m_di.resize(m_charges.size(),0.);
  m_dj.resize(m_photons.size(),0.);
  m_dij.resize(m_charges.size());
  for (size_t i(0);i<m_dij.size();++i) m_dij[i].resize(m_photons.size(),0.);

  if (!m_charges.size() || !m_photons.size()) m_on=false;
  if (!m_on) msg_Debugging()<<"switched off"<<std::endl;
}

Particle_Dresser::~Particle_Dresser()
{
}

void Particle_Dresser::CompleteConeLists()
{
  if (m_kfdR2s.empty()) return;
  DEBUG_FUNC("");
  for (size_t i(0);i<m_dR2.size();++i) {
    kf_code kf(p_fl[m_charges[i]].Kfcode());
    if (m_kfdR2s.find(kf)!=m_kfdR2s.end()) m_dR2[i]=m_kfdR2s[kf];
    msg_Debugging()<<i<<": "<<kf<<" -> dR="<<sqrt(m_dR2[i])<<std::endl;
  }
}

Vec4D_Vector Particle_Dresser::Dress(const Vec4D_Vector& p)
{
  DEBUG_FUNC("N_P="<<m_photons.size()<<", N_C="<<m_charges.size());
  if (!m_on) return p;
  switch (m_algo) {
  case 0:
    return ConeDress(p);
  case 1:
    return RecombinationDress(p);
  default:
    THROW(fatal_error,"Unknown dressing algorithm.");
    return p;
  }
  return p;
}

Vec4D_Vector Particle_Dresser::ConeDress(const Vec4D_Vector& p)
{
  Vec4D_Vector pp(p);
  size_t n(m_photons.size());
  std::vector<bool> valid(n,true);
  double maxd(std::numeric_limits<double>::max()),dmin(maxd);
  size_t ii(0),jj(0),max(std::numeric_limits<size_t>::max());
  // calculate initial dijs=dR(i,j)^2/dR_i^2
  for (size_t i(0);i<m_charges.size();++i) {
    for (size_t j(0);j<m_photons.size();++j) {
      double dij(m_dij[i][j]=DeltaR2(pp[m_charges[i]],
                                     pp[m_photons[j]])/m_dR2[i]);
      if (dij<dmin) { dmin=dij; ii=i; jj=j; }
    }
  }
  while (dmin<1.) {
    if (msg_LevelIsDebugging()) {
      msg_Out()<<"ktij: ";
      for (size_t i(0);i<m_dij.size();++i) {
        msg_Out()<<m_dij[i]<<"\n      ";
      }
      msg_Out()<<"-> i: "<<ii<<" , j: "<<jj<<" , dmin="<<dmin<<std::endl;
    }
    // mark photon that is recombined
    valid[jj]=false;
    // recombine, do not recompute always with respect to bare axis
    pp[m_charges[ii]]+=pp[m_photons[jj]];
    pp[m_photons[jj]]=Vec4D(0.,0.,0.,0.);
    for (size_t i(0);i<m_charges.size();++i) m_dij[i][jj]=maxd;
    // find new dmin
    dmin=maxd;
    for (size_t i(0);i<m_charges.size();++i) {
      for (size_t j(0);j<m_photons.size();++j) if (valid[j]) {
        double dij(m_dij[i][j]);
        if (dij<dmin) { dmin=dij; ii=i; jj=j; }
      }
    }
  }
  return pp;
}

Vec4D_Vector Particle_Dresser::RecombinationDress(const Vec4D_Vector& p)
{
  Vec4D_Vector pp(p);
  size_t n(m_photons.size());
  std::vector<bool> valid(n,true);
  double maxd(std::numeric_limits<double>::max()),dmin(maxd);
  size_t ii(0),jj(0),max(std::numeric_limits<size_t>::max());
  // calculate initial di, dj, dijs
  for (size_t i(0);i<m_charges.size();++i) {
    double di(Pow(pp[m_charges[i]].PPerp2(),m_exp));
    m_di[i]=di;
    for (size_t j(0);j<m_photons.size();++j) if (valid[j]) {
      double dj(Pow(pp[m_photons[j]].PPerp2(),m_exp));
      double dij(Min(di,dj)*DeltaR2(pp[m_charges[i]],
                                    pp[m_photons[j]])/m_dR2[i]);
      m_dj[j]=dj;
      m_dij[i][j]=dij;
      if (dj<dmin)  { dmin=dj;  ii=max; jj=j; }
      if (dij<dmin) { dmin=dij; ii=i;   jj=j; }
    }
  }
  while (true) {
    if (msg_LevelIsDebugging()) {
      msg_Out()<<"ktj:  "<<m_dj<<std::endl;
      msg_Out()<<"ktij: ";
      for (size_t i(0);i<m_dij.size();++i) {
        msg_Out()<<m_dij[i]<<"\n      ";
      }
      msg_Out()<<"-> i: "<<ii<<" , j: "<<jj<<" , dmin="<<dmin<<std::endl;
    }
    // mark photon that is either recombined or removed
    valid[jj]=false;
    // if dmin is dij, then recombine, recompute
    if (ii<max) {
      pp[m_charges[ii]]+=pp[m_photons[jj]];
      pp[m_photons[jj]]=Vec4D(0.,0.,0.,0.);
      m_di[ii]=Pow(pp[m_charges[ii]].PPerp2(),m_exp);
      m_dj[jj]=maxd;
      for (size_t i(0);i<m_charges.size();++i) m_dij[i][jj]=maxd;
      for (size_t j(0);j<m_photons.size();++j) if (valid[j])
        m_dij[ii][j]=Min(m_di[ii],m_dj[j])
                      *DeltaR2(pp[m_charges[ii]],
                               pp[m_photons[j]])/m_dR2[ii];
    }
    // if dmin is dj, then remove
    else {
      m_dj[jj]=maxd;
      for (size_t i(0);i<m_charges.size();++i) m_dij[i][jj]=maxd;
    }
    n--;
    // if no photon left, nothing more to do
    if (n==0) break;
    // else find new dmin
    dmin=maxd;
    for (size_t i(0);i<m_charges.size();++i) {
      for (size_t j(0);j<m_photons.size();++j) if (valid[j]) {
        double dj(m_dj[j]),dij(m_dij[i][j]);
        if (dj<dmin)  { dmin=dj;  ii=max; jj=j; }
        if (dij<dmin) { dmin=dij; ii=i;   jj=j; }
      }
    }
  }
  return pp;
}

double Particle_Dresser::Pow(const double& x, const double& exp)
{
  if      (exp== 0.) return 1.;
  else if (exp== 1.) return x;
  else if (exp==-1.) return 1./x;
  else               return std::pow(x,exp);
}

double Particle_Dresser::DeltaPhi(const Vec4D& p1, const Vec4D& p2)
{
  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  return acos((p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2));
}

double Particle_Dresser::DeltaR2(const Vec4D& p1, const Vec4D& p2)
{
  return sqr(p1.Y()-p2.Y()) + sqr(DeltaPhi(p1,p2));
}

