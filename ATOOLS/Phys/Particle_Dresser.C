#include "ATOOLS/Phys/Particle_Dresser.H"

#include <limits>
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ATOOLS;

Particle_Dresser::Particle_Dresser(const Flavour * fl,
                                   const size_t& nin, const size_t& nout,
                                   std::string algo, double dR) :
  m_on(true), p_sl(new Selector_List()),
  m_algo(0), m_exp(1.), m_dR2global(dR*dR)
{
  // This constructor sets up the dressing algorithm parameters
  // and the flavour config assumed constant
  DEBUG_FUNC(nin<<" -> "<<nout<<", algo="<<algo<<", dR="<<dR);
  SetAlgo(algo);

  for (size_t i(nin);i<nin+nout;++i) {
    p_sl->AddFlavour(fl[i]);
    if (p_sl->back().Flavour().IsPhoton()) m_photons.push_back(i);
    if (p_sl->back().Flavour().Charge())   m_charges.push_back(i);
  }
  p_sl->SetNIn(nin);
  m_dR2.resize(m_charges.size(),m_dR2global);
  m_di.resize(m_charges.size(),0.);
  m_dj.resize(m_photons.size(),0.);
  m_dij.resize(m_charges.size());
  for (size_t i(0);i<m_dij.size();++i) m_dij[i].resize(m_photons.size(),0.);
  if (!m_charges.size() || !m_photons.size()) m_on=false;
  if (!m_on) msg_Debugging()<<"switched off"<<std::endl;
}

Particle_Dresser::Particle_Dresser(std::string algo, double dR) :
  m_on(true), p_sl(NULL),
  m_algo(0), m_exp(1.), m_dR2global(dR*dR)
{
  // This constructor sets up the dressing algorithm parameters only
  DEBUG_FUNC("algo="<<algo<<", dR="<<dR);
  SetAlgo(algo);
}

Particle_Dresser::~Particle_Dresser()
{
  if (p_sl) delete p_sl;
}

void Particle_Dresser::SetAlgo(std::string algo)
{
  if      (algo=="Cone")          { m_algo=0; }
  else if (algo=="kt")            { m_algo=1; m_exp=1.; }
  else if (algo=="CA")            { m_algo=1; m_exp=0.; }
  else if (algo=="antikt")        { m_algo=1; m_exp=-1.; }
  else                            { m_algo=100; }
}

void Particle_Dresser::CompleteConeLists()
{
  if (m_kfdR2s.empty()) return;
  DEBUG_FUNC("");
  for (size_t i(0);i<m_dR2.size();++i) {
    kf_code kf((*p_sl)[m_charges[i]].Flavour().Kfcode());
    if (m_kfdR2s.find(kf)!=m_kfdR2s.end()) m_dR2[i]=m_kfdR2s[kf];
    msg_Debugging()<<i<<": "<<kf<<" -> dR="<<sqrt(m_dR2[i])<<std::endl;
  }
}

Vec4D_Vector Particle_Dresser::Dress(const Vec4D_Vector& p)
{
  // This method assumes there is a constant flavour config set in p_sl,
  // only update momenta
  DEBUG_FUNC("N_P="<<m_photons.size()<<", N_C="<<m_charges.size());
  if (!m_on) return p;
  p_sl->SetMomenta(p);
  switch (m_algo) {
  case 0:
    ConeDress(*p_sl); break;
  case 1:
    RecombinationDress(*p_sl); break;
  default:
    THROW(fatal_error,"Unknown dressing algorithm."); break;
  }
  return p_sl->ExtractMomenta();
}

void Particle_Dresser::Dress(Selector_List &sl)
{
  // This method assumes no information about the process and
  // operates on the flavours and momenta given
  m_photons.clear();
  m_charges.clear();
  for (size_t i(sl.NIn());i<sl.size();++i) {
    if (sl[i].Flavour().IsPhoton()) m_photons.push_back(i);
    if (sl[i].Flavour().Charge())   m_charges.push_back(i);
  }
  DEBUG_FUNC("N_P="<<m_photons.size()<<", N_C="<<m_charges.size());
  m_dR2.resize(m_charges.size(),m_dR2global);
  m_di.resize(m_charges.size(),0.);
  m_dj.resize(m_photons.size(),0.);
  m_dij.resize(m_charges.size());
  for (size_t i(0);i<m_dij.size();++i) m_dij[i].resize(m_photons.size(),0.);
  if (!m_charges.size() || !m_photons.size()) m_on=false;
  if (!m_on) msg_Debugging()<<"switched off"<<std::endl;
  if (!m_on) return;
  msg_Debugging()<<sl<<std::endl;

  // dress
  switch (m_algo) {
  case 0:
    ConeDress(sl); break;
  case 1:
    RecombinationDress(sl); break;
  default:
    THROW(fatal_error,"Unknown dressing algorithm."); break;
  }

  // remove clustered photons
  for (Selector_List::iterator it(sl.begin());it<sl.end();) {
    if (it->Momentum()==Vec4D(0.,0.,0.,0.)) it=sl.erase(it);
    else ++it;
  }
  msg_Debugging()<<sl<<std::endl;
}

void Particle_Dresser::ConeDress(Selector_List& sl)
{
  // implements cone dressing with flavour dependent cone sizes dR
  // momenta are added, cone axis remains constant
  size_t n(m_photons.size());
  std::vector<bool> valid(n,true);
  double maxd(std::numeric_limits<double>::max()),dmin(maxd);
  size_t ii(0),jj(0),max(std::numeric_limits<size_t>::max());
  // calculate initial dijs=dR(i,j)^2/dR_i^2
  for (size_t i(0);i<m_charges.size();++i) {
    for (size_t j(0);j<m_photons.size();++j) {
      double dij(m_dij[i][j]=DeltaR2(sl[m_charges[i]].Momentum(),
                                     sl[m_photons[j]].Momentum())/m_dR2[i]);
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
    sl[m_charges[ii]].SetMomentum(sl[m_charges[ii]].Momentum()
                                  +sl[m_photons[jj]].Momentum());
    sl[m_photons[jj]].SetMomentum(Vec4D(0.,0.,0.,0.));
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
}

void Particle_Dresser::RecombinationDress(Selector_List &sl)
{
  // implements sequential recombination dressing with flavour
  // dependent angular parameter dR
  // momenta are added, dij are updated, operates until either
  // all photons removed or clustered
  size_t n(m_photons.size());
  std::vector<bool> valid(n,true);
  double maxd(std::numeric_limits<double>::max()),dmin(maxd);
  size_t ii(0),jj(0),max(std::numeric_limits<size_t>::max());
  // calculate initial di, dj, dijs
  for (size_t i(0);i<m_charges.size();++i) {
    double di(Pow(sl[m_charges[i]].Momentum().PPerp2(),m_exp));
    m_di[i]=di;
    for (size_t j(0);j<m_photons.size();++j) if (valid[j]) {
      double dj(Pow(sl[m_photons[j]].Momentum().PPerp2(),m_exp));
      double dij(Min(di,dj)*DeltaR2(sl[m_charges[i]].Momentum(),
                                    sl[m_photons[j]].Momentum())/m_dR2[i]);
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
      sl[m_charges[ii]].SetMomentum(sl[m_charges[ii]].Momentum()
                                    +sl[m_photons[jj]].Momentum());
      sl[m_photons[jj]].SetMomentum(Vec4D(0.,0.,0.,0.));
      m_di[ii]=Pow(sl[m_charges[ii]].Momentum().PPerp2(),m_exp);
      m_dj[jj]=maxd;
      for (size_t i(0);i<m_charges.size();++i) m_dij[i][jj]=maxd;
      for (size_t j(0);j<m_photons.size();++j) if (valid[j])
        m_dij[ii][j]=Min(m_di[ii],m_dj[j])
                      *DeltaR2(sl[m_charges[ii]].Momentum(),
                               sl[m_photons[j]].Momentum())/m_dR2[ii];
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

