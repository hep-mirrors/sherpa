#include "AMISIC++/Perturbative/MI_Process.H"
#include "AMISIC++/Perturbative/MI_Integrator.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE AMISIC::XS_Base
#define PARAMETER_TYPE ATOOLS::Flavour_Vector
#include "ATOOLS/Org/Getter_Function.C"

using namespace AMISIC;
using namespace ATOOLS;

// This is the base class for the "naked" matrix elements, which do not need to
// know their flavours, as this is part of their explicit implementation.
// We've reimplemented them to make sure we can call them depending on
// Mandelstam variables, store the results, and they are also equipped with a
// somewhat improved colour handling.

XS_Base::XS_Base() :
  m_name(""), m_Ms(0.), m_Mt(0.), m_Mu(0.), m_lastxs(0.) {
  m_masses.resize(4);
  m_masses2.resize(m_masses.size());
  m_colours.resize(m_masses.size());
  for (size_t i=0;i<m_masses.size();i++) {
    m_colours[i].resize(2);
    m_colours[i][0] = m_colours[i][1] = 0;
    m_masses2[i]    = m_masses[i] = 0;
  }
}


XS_Base::XS_Base(const std::vector<double> & masses) :
  m_name(""), m_Ms(0.), m_Mt(0.), m_Mu(0.), m_lastxs(0.), m_masses(masses) {
  m_masses2.resize(m_masses.size());
  m_colours.resize(m_masses.size());
  for (size_t i=0;i<m_masses.size();i++) {
    m_colours[i].resize(2);
    m_colours[i][0] = m_colours[i][1] = 0;
    m_masses2[i]    = m_masses[i] = 0;
  }
}

XS_Base::~XS_Base() = default;

///////////////////////////////////////////////////////////////////////////////
// MI_Process
///////////////////////////////////////////////////////////////////////////////

// The base class for the parton-level processes with the matrix elements
// only externally set.  This makes sure we calculate the cross section
// for every generic ME only once and do not iterate over flavour
// permutations unless necessary.
// The MI_Process class knows the flavours and is mainly used as container
// and to construct scattering kinematics, i.e. the four-vectors of the
// outgoing particles.

MI_Process::MI_Process(const std::vector<Flavour>& flavs)
    : m_name(flavs[0].IDName() + " " + flavs[1].IDName() + " --> " +
             flavs[2].IDName() + " " + flavs[3].IDName()),
      m_stretcher(Momenta_Stretcher(std::string("AMISIC: ") + m_name)),
      p_me2(NULL), m_Emin((*mipars)("E_min")),
      m_masslessIS((flavs[0].Kfcode() < 4 || flavs[0].Kfcode() == 21) &&
                   (flavs[1].Kfcode() < 4 || flavs[1].Kfcode() == 21))
{
  if (flavs.size() != 4) {
    msg_Error() << "Error in " << METHOD << ":\n"
                << "   Tried to initialize MPI process with wrong number of "
                << "flavours = " << m_flavs.size() << " --> " << m_name
                << ".\n";
    THROW(fatal_error,
          "Tried to initialize MPI process with wrong number of flavours.");
  }
  m_flavs.resize(flavs.size());
  m_momenta.resize(m_flavs.size());
  m_masses.resize(m_flavs.size());
  m_masses2.resize(m_flavs.size());
  for (size_t i = 0; i < m_flavs.size(); i++) {
    m_flavs[i]   = flavs[i];
    m_masses[i]  = flavs[i].Mass();
    m_masses2[i] = sqr(m_masses[i]);
    if (i < 2) m_sumInMasses += m_masses[i];
    else
      m_sumOutMasses += m_masses[i];
  }
}

MI_Process::~MI_Process() = default;

bool MI_Process::MakeKinematics(
        MI_Integrator*                                         integrator,
        std::array<std::shared_ptr<REMNANTS::Remnant_Base>, 2> remnants)
{
  if (!AllowedKinematics(sqrt(integrator->SHat()))) return false;
  ///////////////////////////////////////////////////////////////////////////
  // Until now only have massless initial state partons.
  ///////////////////////////////////////////////////////////////////////////
  if (!m_masslessIS) return false;
  MasslessISKinematics(integrator);
  ///////////////////////////////////////////////////////////////////////////
  // If the final state is massive, we use the momenta stretcher to push
  // particles onto their mass shells.  The logic is to go to the c.m. system
  // of the scatter, rescale momenta there, and boost back.
  ///////////////////////////////////////////////////////////////////////////
  if (m_flavs[2].Kfcode() == 5 || m_flavs[2].Kfcode() == 4) {
    Vec4D    cms = m_momenta[0] + m_momenta[1];
    Poincare scattercms(cms);
    for (size_t i = 2; i < m_momenta.size(); i++)
      scattercms.Boost(m_momenta[i]);
    if (!m_stretcher.ZeroThem(2, m_momenta) ||
        !m_stretcher.MassThem(2, m_momenta, m_masses))
      return false;
    for (size_t i = 2; i < m_momenta.size(); i++)
      scattercms.BoostBack(m_momenta[i]);
  }
  return (m_momenta[0][0] > Max(m_flavs[0].HadMass(), m_Emin) &&
          m_momenta[1][0] > Max(m_flavs[1].HadMass(), m_Emin) &&
          AllowedRemnants(remnants));
}

const bool MI_Process::AllowedRemnants(
        std::array<std::shared_ptr<REMNANTS::Remnant_Base>, 2> remnants) const
{
  ///////////////////////////////////////////////////////////////////////////
  // Make sure there is enough energy left in the remnants
  ///////////////////////////////////////////////////////////////////////////
  // msg_Out()<<"   * "<<METHOD<<" for "
  //	   <<m_flavs[0]<<" ("<<m_momenta[0][0]<<" vs. "
  //	   <<(*remnants)[0]->ResidualE()<<") and "
  //	   <<m_flavs[1]<<" ("<<m_momenta[1][0]<<" vs. "
  //	   <<(*remnants)[1]->ResidualE()<<").\n";
  return (remnants[0]->TestExtract(m_flavs[0], m_momenta[0]) &&
          remnants[1]->TestExtract(m_flavs[1], m_momenta[1]));
}

void MI_Process::MasslessISKinematics(MI_Integrator* integrator)
{
  ///////////////////////////////////////////////////////////////////////////
  // Reconstruct kinematics from transverse momentum of outgoing particles
  // and their individual rapidities
  ///////////////////////////////////////////////////////////////////////////
  double pt2 = integrator->PT2(), pt = sqrt(pt2);
  double y3 = integrator->Y(0), y4 = integrator->Y(1);
  double mt2 = sqrt(pt2 + m_masses2[2]), mt3 = sqrt(pt2 + m_masses2[3]);
  double phi    = 2. * M_PI * ran->Get();
  double cosphi = cos(phi), sinphi = sin(phi);
  m_momenta[2] =
          Vec4D(mt2 * cosh(y3), pt * cosphi, pt * sinphi, mt2 * sinh(y3));
  m_momenta[3] =
          Vec4D(mt3 * cosh(y4), -pt * cosphi, -pt * sinphi, mt3 * sinh(y4));
  double E = m_momenta[2][0] + m_momenta[3][0];
  double p = m_momenta[2][3] + m_momenta[3][3];
  ///////////////////////////////////////////////////////////////////////////
  // correspond to momenta
  // p0 = (mt2*exp(y3)+mt3*exp(y4))/2 * (1, 0, 0, 1) and
  // p1 = (mt2*exp(-y3)+mt3*exp(-y4))/2 * (1, 0, 0, -1)
  ///////////////////////////////////////////////////////////////////////////
  m_momenta[0] = (E + p) / 2. * Vec4D(1, 0, 0, 1);
  m_momenta[1] = (E - p) / 2. * Vec4D(1, 0, 0, -1);
  // msg_Out()<<"   * "<<METHOD<<"(y1,2 = "<<y3<<", "<<y4<<", pt^2 = "<<pt2<<"):
  // "
  //	   <<"x1 = "<<(2.*m_momenta[0][0]/rpa->gen.Ecms())<<", "
  //	   <<"x2 = "<<(2.*m_momenta[1][0]/rpa->gen.Ecms())<<".\n";
}

Particle* MI_Process::GetParticle(const size_t& i)
{
  Particle* part =
          new Particle(-1, m_flavs[i], m_momenta[i], (i < 2 ? 'I' : 'F'));
  part->SetNumber();
  for (size_t j = 0; j < 2; j++) part->SetFlow(j + 1, p_me2->Colour(i, j));
  return part;
}
