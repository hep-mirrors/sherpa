#include "AddOns/OpenLoops/OpenLoops_Interface.H"
#include "AddOns/OpenLoops/Color_Correlated_ME2.H"
#include "AddOns/OpenLoops/OpenLoops_Born.H"

using namespace OpenLoops;

Color_Correlated_ME2::
Color_Correlated_ME2(const PHASIC::External_ME_Args& args,
		     int ol_id, const AmplitudeType& type) :
  PHASIC::Color_Correlated_ME2(args), m_amptype(type), m_ol_id(ol_id),
  m_dim(args.m_inflavs.size()+args.m_outflavs.size())
{
  m_symfac =ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
  m_symfac*=ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
  
  m_ccmatrix = new double[m_dim*m_dim];
}


Color_Correlated_ME2::~Color_Correlated_ME2()
{
  delete [] m_ccmatrix;
}


void Color_Correlated_ME2::Calc(const ATOOLS::Vec4D_Vector &p)
{
  OpenLoops_Interface::SetParameter("alpha",  AlphaQED());
  OpenLoops_Interface::SetParameter("alphas", AlphaQCD());
  OpenLoops_Interface::PopulateColorCorrelatorMatrix(m_ol_id,p, m_born2, m_ccmatrix, AmpType());
}


double Color_Correlated_ME2::GetValue(const size_t& i, const size_t& j) const
{
  /* Fortran representation of 2D matrix is via a column-major.
     Stick to the convention of Tree_ME2_Base etc and normalize in
     such a way that the symmetry factor is NOT taken into account */
  return m_ccmatrix[ j*m_dim + i ] * m_symfac;
}


double Color_Correlated_ME2::GetBorn2() const
{
  /* Symmetry factor not to be taken into account here */
  return m_born2 * m_symfac;
}

    
DECLARE_COLORCORRELATEDME2_GETTER(Color_Correlated_ME2,
				  "OpenLoops::Color_Correlated_ME2")


PHASIC::Color_Correlated_ME2 *ATOOLS::Getter
<PHASIC::Color_Correlated_ME2,
 PHASIC::External_ME_Args,
 Color_Correlated_ME2>::
operator()(const PHASIC::External_ME_Args &args) const
{
  OpenLoops_Interface::SetParameter("coupling_qcd_0", (int) args.m_orders[0]);
  OpenLoops_Interface::SetParameter("coupling_qcd_1", 0);
  OpenLoops_Interface::SetParameter("coupling_ew_0",  (int) args.m_orders[1]);
  OpenLoops_Interface::SetParameter("coupling_ew_1",  0);

  if (args.m_source != "" && args.m_source != "OpenLoops")
    return NULL;

  AmplitudeType types[2] = {Loop2, Tree};
  for (size_t i=0; i<2; ++i) {
    int id = OpenLoops_Interface::RegisterProcess(args.m_inflavs, 
						  args.m_outflavs, 
						  (int)(types[i]));
    if (id>0) return new Color_Correlated_ME2(args, id, types[i]);
  }
  return NULL;
}
