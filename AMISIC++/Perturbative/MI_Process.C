#include "AMISIC++/Perturbative/MI_Process.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

// This is the base class for the "naked" matrix elements, which do not need to
// know their flavours, as this is part of their explicit implementation.
// We've reimplemented them to make sure we can call them depending on
// Mandelstam variables, store the results, and they are also equipped with a
// somewhat improved colour handling.

XS_Base::XS_Base() :
  m_name(""), m_lastxs(0.), m_Ms(0.), m_Mt(0.), m_Mu(0.) {
  m_masses.resize(4);
  m_masses2.resize(m_masses.size());
  m_colours.resize(m_masses.size());
  for (size_t i=0;i<m_masses.size();i++) {
    m_colours[i].resize(2);
    m_colours[i][0] = m_colours[i][1] = 0;
    m_masses2[i]    = m_masses[i] = 0;
  }
}
  

XS_Base::XS_Base(const vector<double> & masses) :
  m_name(""), m_lastxs(0.), m_Ms(0.), m_Mt(0.), m_Mu(0.), m_masses(masses) {
  m_masses2.resize(m_masses.size());
  m_colours.resize(m_masses.size());
  for (size_t i=0;i<m_masses.size();i++) {
    m_colours[i].resize(2);
    m_colours[i][0] = m_colours[i][1] = 0;
    m_masses2[i]    = m_masses[i] = 0;
  }
}

XS_Base::~XS_Base() {}

///////////////////////////////////////////////////////////////////////////////
// MI_Process
///////////////////////////////////////////////////////////////////////////////

// This is the base class for the parton-level processes with the matrix elements
// only externally set.  This makes sure we calculate the cross section for every
// generic ME only once and do not iterate over flavour permutations unless necessary.
// The MI_Process class knows the flavours and is mainly used as container and to
// construct scattering kinematics, i.e. the four-vector of the outgoing particles.

MI_Process::MI_Process(const vector<Flavour> & flavs) :
  m_name(flavs[0].IDName()+" "+flavs[1].IDName()+" --> "+
	 flavs[2].IDName()+" "+flavs[3].IDName()),
  p_me2(NULL), m_emin(0.),
  m_massless(flavs[0].HadMass()<1.0 && flavs[1].HadMass()<1.0 &&
	     flavs[2].HadMass()<1.0 && flavs[3].HadMass()<1.0),
  m_masslessIS((flavs[0].Kfcode()<4 || flavs[0].Kfcode()==21) &&
	       (flavs[1].Kfcode()<4 || flavs[1].Kfcode()==21))
{
  if (flavs.size()!=4) {
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   Tried to initialize MPI process with wrong number of "
	       <<"flavours ("<<m_flavs.size()<<" --> "<<m_name<<".\n"
	       <<"   Will exit the run.\n";
    exit(1);
  }
  m_flavs.resize(flavs.size());
  m_momenta.resize(m_flavs.size());
  m_masses.resize(m_flavs.size());
  m_PSmasses.resize(m_flavs.size());
  m_masses2.resize(m_flavs.size());
  for (size_t i=0;i<m_flavs.size();i++) {
    m_flavs[i]    = flavs[i];
    m_masses[i]   = flavs[i].Mass();
    m_PSmasses[i] = flavs[i].HadMass();
    m_masses2[i]  = sqr(m_masses[i]);
  }
}

MI_Process::~MI_Process() {}

bool MI_Process::MakeKinematics(const double & pt2,
				const double & y3,const double & y4,
				const double & Ehat) {
  if (!AllowedKinematics(Ehat)) return false;
  double phi = 2.*M_PI*ran->Get();
  // Until now we only have massless initial state partons.
  if (m_masslessIS) MasslessKinematics(pt2,phi,y3,y4);
  else return false;
  // If the final state is massive, we use the momenta stretcher to push
  // particles onto their mass shells.  The logic is to go to the c.m. system
  // of the scatter, rescale momenta there, and boost back.
  if (m_flavs[2].Kfcode()==5 || m_flavs[2].Kfcode()==4) {
    Vec4D cms = m_momenta[0]+m_momenta[1];
    Poincare scattercms(cms);
    for (size_t i=2;i<m_momenta.size();i++) scattercms.Boost(m_momenta[i]);
    if (!m_stretcher.ZeroThem(2,m_momenta) ||
	!m_stretcher.MassThem(2,m_momenta,m_masses)) {
      return false;
    }
    for (size_t i=2;i<m_momenta.size();i++) {
      scattercms.BoostBack(m_momenta[i]);
    }
  }
  return (m_momenta[0][0] > m_emin && m_momenta[1][0] > m_emin);
}

bool MI_Process::AllowedKinematics(const double & Ehat) {
  // making sure that the c.m. energy of the scatter is larger than the
  // IS or FS sum of masses.
  if (m_massless) return true;
  return (m_PSmasses[0]+m_PSmasses[1]<Ehat &&
	  m_PSmasses[2]+m_PSmasses[3]<Ehat); 
}

void MI_Process::MasslessKinematics(const double & pt2,const double & phi,
				    const double & y3,const double & y4) {
  // reconstruct kinematics from transverse momentum of outgoing particles and
  // their individual rapidities
  double mt2 = sqrt(pt2+m_masses2[2]);
  double mt3 = sqrt(pt2+m_masses2[3]);
  double pt  = sqrt(pt2);
  double cosphi = cos(phi), sinphi = sin(phi);
  m_momenta[2] = Vec4D(mt2*cosh(y3), pt*cosphi, pt*sinphi, mt2*sinh(y3));
  m_momenta[3] = Vec4D(mt3*cosh(y4),-pt*cosphi,-pt*sinphi, mt3*sinh(y4));
  double E = m_momenta[2][0]+m_momenta[3][0];
  double p = m_momenta[2][3]+m_momenta[3][3];
  m_momenta[0] = (E+p)/2. * Vec4D(1, 0, 0, 1);
  m_momenta[1] = (E-p)/2. * Vec4D(1, 0, 0,-1);
}

Particle * MI_Process::GetParticle(const size_t & i) {
  Particle * part = new Particle(-1,m_flavs[i],m_momenta[i],(i<2?'I':'F'));
  part->SetNumber();
  for(size_t j=0;j<2;j++) part->SetFlow(j+1,p_me2->Colour(i,j));
  return part;
}

