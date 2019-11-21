#include "AddOns/OpenLoops/OpenLoops_Born.H"
#include "AddOns/OpenLoops/OpenLoops_Interface.H"

#include "PHASIC++/Process/External_ME_Args.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Library_Loader.H"

using namespace OpenLoops;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

OpenLoops_Born::OpenLoops_Born(const External_ME_Args& args,
                               int ol_id,
                               AmplitudeType type) :
  Tree_ME2_Base(args), m_ol_id(ol_id), m_amplitudetype(type)
{
  m_symfac =Flavour::FSSymmetryFactor(args.m_outflavs);
  m_symfac*=Flavour::ISSymmetryFactor(args.m_inflavs);
  
  m_order_qcd = args.m_orders[0];
  m_order_ew  = args.m_orders[1];
}

double OpenLoops_Born::Calc(const Vec4D_Vector& momenta)
{
  OpenLoops_Interface::SetParameter("alpha", AlphaQED());
  OpenLoops_Interface::SetParameter("alphas", AlphaQCD());

  double result(0.0);
  switch (m_amplitudetype) {
    case Tree:
      OpenLoops_Interface::EvaluateTree(m_ol_id, momenta, result);
      break;
    case Loop2:
      OpenLoops_Interface::EvaluateLoop2(m_ol_id, momenta, result);
      break;
  }

  // OL returns ME2 including 1/symfac, convention
  // is not to include it in Tree_ME2_Base, however
  return m_symfac*result;
}

DECLARE_TREEME2_GETTER(OpenLoops_Born,
		       "OpenLoops_Born")

Tree_ME2_Base* ATOOLS::Getter<Tree_ME2_Base,
			      External_ME_Args,
			      OpenLoops_Born>::
operator()(const External_ME_Args &args) const
{
  if (!args.m_source.empty() && args.m_source != "OpenLoops")
    return NULL;

  OpenLoops_Interface::SetParameter("coupling_qcd_0", (int) args.m_orders[0]);
  OpenLoops_Interface::SetParameter("coupling_qcd_1", 0);
  OpenLoops_Interface::SetParameter("coupling_ew_0" , (int) args.m_orders[1]);
  OpenLoops_Interface::SetParameter("coupling_ew_1" , 0);

  AmplitudeType types[2] = {Loop2, Tree};
  for (size_t i=0; i<2; ++i) 
    {
      int id = OpenLoops_Interface::RegisterProcess(args.m_inflavs, 
						    args.m_outflavs, 
						    (int)(types[i]));
      if (id>0) return new OpenLoops_Born(args, id, types[i]);
    }

  return NULL;
}

int OpenLoops_Born::OrderQCD(const int &id) const
{ return m_order_qcd; }

int OpenLoops_Born::OrderEW(const int &id) const
{ return m_order_ew; }
