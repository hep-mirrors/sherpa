#include "AddOns/OpenLoops/OpenLoops_Born.H"

#include "AddOns/OpenLoops/OpenLoops_Interface.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Library_Loader.H"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

namespace OpenLoops {

OpenLoops_Born::OpenLoops_Born(const Process_Info& pi,
                               const Flavour_Vector& flavs,
                               int ol_id, int amptype) :
  Tree_ME2_Base(pi, flavs), m_ol_id(ol_id), m_amptype(amptype)
{
  m_symfac=pi.m_fi.FSSymmetryFactor();
  m_symfac*=pi.m_ii.ISSymmetryFactor();
}

double OpenLoops_Born::Calc(const Vec4D_Vector& momenta)
{
  OpenLoops_Interface::SetParameter("alpha", AlphaQED());
  OpenLoops_Interface::SetParameter("alphas", AlphaQCD());
  
  double result;
  if (m_amptype==1) OpenLoops_Interface::EvaluateTree(m_ol_id, momenta, result);
  if (m_amptype==12) OpenLoops_Interface::EvaluateLoop2(m_ol_id, momenta, result);

  // OL returns ME2 including 1/symfac, but Calc is supposed to return it
  // without 1/symfac, thus multiplying with symfac here
  return m_symfac*result;
}

}

using namespace OpenLoops;

DECLARE_TREEME2_GETTER(OpenLoops_Born,"OpenLoops_Born")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,OpenLoops_Born>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="OpenLoops") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo && pi.m_fi.m_nloewtype!=nlo_type::real) return NULL;
  if (pi.m_fi.m_nloqcdtype!=nlo_type::lo && pi.m_fi.m_nloqcdtype!=nlo_type::real) return NULL;

  OpenLoops_Interface::SetParameter("coupling_qcd_0", (int) pi.m_maxcpl[0]);
  OpenLoops_Interface::SetParameter("coupling_qcd_1", 0);
  OpenLoops_Interface::SetParameter("coupling_ew_0", (int) pi.m_maxcpl[1]);
  OpenLoops_Interface::SetParameter("coupling_ew_1", 0);

  int born_types[2] = {12, 1};
  for (size_t i=0; i<2; ++i) {
    int id = OpenLoops_Interface::RegisterProcess(pi.m_ii, pi.m_fi, born_types[i]);
    if (id>0) {
      Flavour_Vector flavs = pi.ExtractFlavours();
      return new OpenLoops_Born(pi, flavs, id, born_types[i]);
    }
  }

  return NULL;
}


