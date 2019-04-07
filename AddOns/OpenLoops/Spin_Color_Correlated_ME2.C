#include "AddOns/OpenLoops/OpenLoops_Interface.H"
#include "AddOns/OpenLoops/Spin_Color_Correlated_ME2.H"
#include "AddOns/OpenLoops/OpenLoops_Born.H"

#include <assert.h>

using namespace OpenLoops;


Spin_Color_Correlated_ME2::
Spin_Color_Correlated_ME2(const PHASIC::External_ME_Args& args,
			  int ol_id, const AmplitudeType& type) :
  PHASIC::Spin_Color_Correlated_ME2(args),
  m_ol_id(ol_id), m_amptype(type)
{
  m_norm =ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
  m_norm*=ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
}


double Spin_Color_Correlated_ME2::
CalcColorCorrelator(const ATOOLS::Vec4D_Vector& born_moms,
		    const size_t& born_ij,
		    const size_t& born_k) const
{
  OpenLoops_Interface::SetParameter("alphas", AlphaQCD());
  OpenLoops_Interface::SetParameter("alpha" , AlphaQED());
  
  /* <1,...,m;a,b| T_ij T_k |b,a;m,...m1> * T_ij^{-2} */
  return OpenLoops_Interface::EvaluateColorCorrelator(m_ol_id,
						      born_moms,
						      born_ij, born_k,
						      AmpType()) * m_norm;
}


double Spin_Color_Correlated_ME2::
CalcSpinCorrelator(const ATOOLS::Vec4D_Vector& born_moms,
		   const ATOOLS::Vec4D& p_tilde,
		   const size_t& born_ij,
		   const size_t& born_k) const
{
  OpenLoops_Interface::SetParameter("alphas", AlphaQCD());
  OpenLoops_Interface::SetParameter("alpha" , AlphaQED());
  
  /* <1,...,m;a,b| ptilde^\mu T_ij T_k ptilde^\nu |b,a;m,...m1> * T_ij^{-2} * ptilde^{-2} */
  return OpenLoops_Interface::EvaluateSpinCorrelator(m_ol_id,
						     born_moms,
						     p_tilde,
						     born_ij, born_k,
						     AmpType()) * m_norm;
}


DECLARE_SPINCOLORCORRELATEDME2_GETTER(Spin_Color_Correlated_ME2,
				      "OpenLoops::Spin_Color_Correlated_ME2")


PHASIC::Spin_Color_Correlated_ME2 *ATOOLS::Getter
<PHASIC::Spin_Color_Correlated_ME2,
 PHASIC::External_ME_Args,
 Spin_Color_Correlated_ME2>::
operator()(const PHASIC::External_ME_Args& args) const
{
  DEBUG_FUNC(this);

  if (!args.m_source.empty() && args.m_source != "OpenLoops")
    return NULL;

  OpenLoops_Interface::SetParameter("coupling_qcd_0", (int) args.m_orders[0]);
  OpenLoops_Interface::SetParameter("coupling_qcd_1", 0);
  OpenLoops_Interface::SetParameter("coupling_ew_0" , (int) args.m_orders[1]);
  OpenLoops_Interface::SetParameter("coupling_ew_1" , 0);
  
  AmplitudeType types[2] = {Loop2, Tree};
  for (size_t i=0; i<2; ++i) {
    int id = OpenLoops_Interface::RegisterProcess(args.m_inflavs, 
						  args.m_outflavs, 
						  (int)(types[i]));
    if (id>0) return new Spin_Color_Correlated_ME2(args, id, types[i]);
  }
  return NULL;
}
