#include "AddOns/OpenLoops/OpenLoops_Born.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

namespace OpenLoops {

OpenLoops_Interface* OpenLoops_Born::s_interface=NULL;

OpenLoops_Born::OpenLoops_Born(const Process_Info& pi,
                               const Flavour_Vector& flavs,
                               Amp2Func amp2,
                               PermutationFunc permutationfunc,
                               std::vector<int> permutation,
                               std::string functag) :
  ME2_Base(pi, flavs),
  m_amp2(amp2), m_permutationfunc(permutationfunc),
  m_permutation(permutation)
{
  m_oew=pi.m_oew;
  m_oqcd=pi.m_oqcd;
}
  
OpenLoops_Born::~OpenLoops_Born()
{
}

//double OpenLoops_Born::Calc(const Vec4D_Vector& momenta)
double OpenLoops_Born::operator()(const Vec4D_Vector& momenta)
{
  Vec4D_Vector m_moms(momenta);

  double alpha_QED=AlphaQED();
  double alpha_S=AlphaQCD();
  double mur2(10000.0);
  s_interface->OpenLoopsInit(mur2, alpha_QED, alpha_S);

  double B(0.0), V_finite(0.0), V_eps(0.0), V_eps2(0.0), I_finite(0.0), I_eps(0.0), I_eps2(0.0);

  m_permutationfunc(&m_permutation[0]);
  m_amp2(&m_moms[0][0], &B, &V_finite, &V_eps, &V_eps2, &I_finite, &I_eps, &I_eps2);

  if (IsZero(V_eps) && IsZero(V_eps2) && IsZero(I_eps) && IsZero(I_eps2)) {
    if (IsZero(B)) return I_finite;
    if (IsZero(I_finite)) return B;
    PRINT_INFO("B!=0 and I_finite!=0. Returning 0.");
    PRINT_VAR(B);
    PRINT_VAR(I_finite);
    return 0.0;
  }
  else {
    PRINT_INFO("Poles non-zero. Returning 0.");
    PRINT_VAR(B);
    PRINT_VAR(V_finite);
    PRINT_VAR(V_eps);
    PRINT_VAR(V_eps2);
    PRINT_VAR(I_finite);
    PRINT_VAR(I_eps);
    PRINT_VAR(I_eps2);
    return 0.0;
  }
}

}
