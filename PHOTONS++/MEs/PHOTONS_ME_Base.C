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
  m_GF(1.16639e-5),
  m_sqrt2(1.41421356237),
  m_i(Complex(0.,1.)),
  p_boost(NULL), p_rot(NULL),
  m_pvv_zero(pvv),
  m_ew(Photons::s_ew_corr),
  m_ew_scheme(Photons::s_ew_scheme),
  m_nnlo_qed(Photons::s_nnlo_qed),
  m_dev(Photons::s_dev), m_xi(Photons::s_xi)
{
  double  MW  = Flavour(kf_Wplus).Mass();
  double  MZ  = Flavour(kf_Z).Mass();
  double  MH  = Flavour(kf_h0).Mass();
  double  GH  = Flavour(kf_h0).Width();
  double  GW  = Flavour(kf_Wplus).Width();
  double  GZ  = Flavour(kf_Z).Width();
  Complex I   = Complex(0.,1.);
  Complex sw2 = 1.-sqr(MW/MZ);
  Complex cw2 = 1.-sw2;
  if (MODEL::s_model->ScalarNumber("WidthScheme")) {
    Complex muW2(MW*(MW-I*GW)), muZ2(MZ*(MZ-I*GZ)), muH2(MH*(MH-I*GH));
    sw2=muW2/muZ2;
    cw2=1.-cw2;
  }
  m_sW = sqrt(std::abs(sw2));
  m_cW = sqrt(std::abs(cw2));
  #ifndef USING__YFS_NNLO
  if (m_nnlo_qed == 1) {
    msg_Error() << METHOD << " YFS NNLO corrections not compiled properly. Please reconfigure "
    << "with option '--enable-yfsnnlo'. Will calculate using NLO corrections for the "
    << "being.\n";
  }
#endif
}

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
    if (pme && pme->Name() != "Collinear_Approximation_FF" && pme->Name() != "Collinear_Approximation_FI") return pme;
  }
  return NULL;
}

PHOTONS_ME_Base * PHOTONS_ME_Base::GetIRsubtractedME
(const std::string& tag, const Particle_Vector_Vector& pvv)
{
  PHOTONS_ME_Base * pme = PHOTONS_ME_Getter::GetObject(tag, pvv);
  if (!pme && (tag == "Collinear_Approximation_FF" || tag == "Collinear_Approximation_FI")) return NULL;
  if (!pme) THROW(fatal_error, "Did not find IR subtracted ME "+tag);
  return pme;
}

