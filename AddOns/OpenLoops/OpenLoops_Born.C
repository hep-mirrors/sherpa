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
  Tree_ME2_Base(pi, flavs),
  m_amp2(amp2), m_permutationfunc(permutationfunc),
  m_permutation(permutation)
{
  m_symfac=pi.m_fi.FSSymmetryFactor();
  m_symfac*=pi.m_ii.ISSymmetryFactor();
  m_oew=pi.m_oew;
  m_oqcd=pi.m_oqcd;
}

OpenLoops_Born::~OpenLoops_Born()
{
}

double OpenLoops_Born::Calc(const Vec4D_Vector& momenta)
{
  Vec4D_Vector m_moms(momenta);

  double alpha_QED=AlphaQED();
  double alpha_S=AlphaQCD();
  double mur2(10000.0);
  s_interface->OpenLoopsInit(mur2, alpha_QED, alpha_S);

  double M2L0;
  std::vector<double> M2L1(3), M2L2(5), IRL1(3), IRL2(5);

  m_permutationfunc(&m_permutation[0]);
  m_amp2(&m_moms[0][0], &M2L0, &M2L1[0], &IRL1[0], &M2L2[0], &IRL2[0]);

  if (IsZero(M2L1[1]) && IsZero(M2L1[2]) &&
      IsZero(IRL1[1]) && IsZero(IRL1[2]) &&
      IsZero(M2L2[1]) && IsZero(M2L2[2]) && IsZero(M2L2[3]) && IsZero(M2L2[4]) &&
      IsZero(IRL2[1]) && IsZero(IRL2[2]) && IsZero(IRL2[3]) && IsZero(IRL2[4])) {
    // OL returns ME2 including 1/symfac, but Calc is supposed to return it
    // without 1/symfac, thus multiplying with symfac here
    if (IsZero(M2L0)) return m_symfac*M2L2[0];
    else return m_symfac*M2L0;
  }
  else {
    PRINT_INFO("Poles non-zero. Returning 0.");
    return 0.0;
  }
}

}
