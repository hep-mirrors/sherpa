#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "MODEL/Main/Model_Base.H"
#include "PHOTONS++/Main/Photons.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHOTONS::PHOTONS_ME_Base
#define PARAMETER_TYPE PHOTONS::Particle_Vector_Vector
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

PHOTONS_ME_Base::PHOTONS_ME_Base(const Particle_Vector_Vector& pvv) :
  m_alpha(Photons::s_alpha),
  m_e(sqrt(4.*M_PI*m_alpha)),
  m_sW(sqrt(std::abs(MODEL::s_model->ComplexConstant("csin2_thetaW")))),
  m_cW(sqrt(std::abs(MODEL::s_model->ComplexConstant("ccos2_thetaW")))),
  m_GF(1.16639e-5),
  m_sqrt2(1.41421356237),
  m_i(Complex(0.,1.)),
  p_boost(NULL), p_rot(NULL),
  m_pvv_zero(pvv) {}

PHOTONS_ME_Base::~PHOTONS_ME_Base() {
  if (p_boost) delete p_boost;
  if (p_rot) delete p_rot;
}

PHOTONS_ME_Base * PHOTONS_ME_Base::GetIRsubtractedME
(const Particle_Vector_Vector& pvv)
{
  PHOTONS_ME_Getter::Getter_List glist(PHOTONS_ME_Getter::GetGetters());
  for (PHOTONS_ME_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    PHOTONS_ME_Base * pme = (*git)->GetObject(pvv);
    if (pme) return pme;
  }
  return NULL;
}

PHOTONS_ME_Base * PHOTONS_ME_Base::GetIRsubtractedME
(const std::string& tag, const Particle_Vector_Vector& pvv)
{
  PHOTONS_ME_Base * pme = PHOTONS_ME_Getter::GetObject(tag, pvv);
  if (!pme) THROW(fatal_error, "Did not find IR subtracted ME "+tag);
  return pme;
}

