#include "PHOTONS++/PhotonSplitter/Splitting_Functions.H"

#include "PHOTONS++/Main/Photons.H" // DO NOT MOVE THIS TO THE HEADER. CAUSES CIRCULAR DEPENDENCIES
#include "ATOOLS/Math/Random.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Spectator::Spectator(const ATOOLS::Vec4D &p, const int &id, const ATOOLS::Flavour &flav) :
      m_id(id), m_p(p), m_flav(flav)
{
  m_charge = m_flav.Charge();
}

// begin Splitting_Function class
Splitting_Function::Splitting_Function(ATOOLS::Particle * splitter, int fla, 
      int flb, int flc, int intspin, const size_t &id) 
      : m_id(id), m_intspin(intspin), m_on(true)
{
  p_splitter = splitter;
  m_flavs[0] = Flavour(fla);
  m_flavs[1] = Flavour(flb);
  m_flavs[2] = Flavour(flc);
  m_mij2 = sqr(m_flavs[0].Mass(1));
  m_mi2  = sqr(m_flavs[1].Mass(1));
  m_mj2  = sqr(m_flavs[2].Mass(1));
  msg_Debugging()<<m_flavs[0]<<"("<<m_mij2<<") -> "
                 <<m_flavs[1]<<"("<<m_mi2<<") "
                 <<m_flavs[2]<<"("<<m_mj2<<")"<<std::endl;
  
  // possiblity to enhance splitting here by multiplying alpha 
  m_alpha = PHOTONS::Photons::s_alpha;
};

double Splitting_Function::Lambda(const double &a, const double &b, const double &c) const
{
  return a*a+b*b+c*c - 2.*(a*b+a*c+b*c);
}

double Splitting_Function::JFF(const double Q2, const double y) const
{
  return (1.-y)*sqr(1.0-m_mi2/Q2-m_mj2/Q2-m_mk2/Q2)/sqrt(Lambda(1.0,m_mij2/Q2,m_mk2/Q2));
}

void Splitting_Function::SetSpec(Spectator *spec)
{
  p_spec = spec;
  m_flspec = spec->GetFlavour();
  m_mk2 = sqr(m_flspec.Mass(1));
}

void Splitting_Function::AddSpec(Spectator *spec)
{
  m_specs.push_back(spec);
}

double Splitting_Function::operator()
  (const double t, const double z, const double y, const double Q2)
{
  if (t < 4*m_mi2) return 0.; // mass cutoff
  double fac, viji, vijk, frac, zm, zp, coupling;
  // only photons massless, so no point including simplified form here. 
  // coupling 
  coupling = 1./m_specs.size() * m_alpha;
  // determinants
  fac = 1.-m_mi2/Q2-m_mj2/Q2-m_mk2/Q2;
  viji = sqr(fac*y)-4.*m_mi2*m_mj2/sqr(Q2); // mistake in Krauss, Schumann (2007), this is correct taken from Catani-Seymour (2002)
  vijk = sqr(2.*m_mk2/Q2+fac*(1.-y))-4.*m_mk2/Q2;
  // check kinematically accessible 
  if (viji<0.0 || vijk<=0.0) return 0.0; // stops y being too close to 1 
  viji = sqrt(viji)/(fac*y+2.*m_mi2/Q2);
  vijk = sqrt(vijk)/(fac*(1.-y));
  frac = (2.*m_mi2/Q2+fac*y)/(2.*(m_mi2/Q2+m_mj2/Q2+fac*y));
  zm = frac*(1.- viji*vijk);
  zp = frac*(1.+ viji*vijk);

  if (zp*zm > z*(1.-z)) { msg_Debugging() << "Rejected due to z limits\n"; return 0; } // z is not in allowed range

  if (m_intspin == 1) {
    return 2.0 * coupling * JFF(Q2,y) * (1.- 2.*(z*(1.-z) - zp*zm));
  }
  else if (m_intspin == 0) {
    return 2.0 * coupling * JFF(Q2,y) * (1.- 2.*(z*(1.-z) - zp*zm));
  }
  else if (m_intspin == 2) {
    // vector 
    return -2.; // not implemented yet 
  }
  else {
    return -2.;
  }
}

double Splitting_Function::OverIntegrated(const double zmin,const double zmax)
{
  SetLimits(zmin,zmax);
  return 2.0/m_specs.size() * (m_zmax-m_zmin) * m_alpha;
}
    
double Splitting_Function::OverEstimated()
{
  return 2.0/m_specs.size() * m_alpha;
}
    
    
double Splitting_Function::Z()
{
  return m_zmin + (m_zmax-m_zmin)*ATOOLS::ran->Get();
}
