#include "METOOLS/HadronCurrents/V_0_Isoscalar3Pi.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"
#include "ATOOLS/Org/Exception.H"
#include <map>

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

V_0_Isoscalar3Pi::V_0_Isoscalar3Pi(const ATOOLS::Flavour_Vector& flavs,
				   const std::vector<int>& indices,
				   const std::string& name) :
  Current_Base(flavs, indices, name), p_ff(NULL)
{
  if (p_i.size()!=3)
    THROW(fatal_error,"V_0_Isoscalar3Pi needs exactly three pions.");
}

V_0_Isoscalar3Pi::~V_0_Isoscalar3Pi() {
  if (p_ff) { delete p_ff; p_ff = NULL; }
}

void V_0_Isoscalar3Pi::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti) {
  // All three pions are spinless, so there is a single helicity entry and the
  // charge conjugate is the same current -- anti is unused.
  Insert(p_ff->Current(moms), size_t(0));
}

void V_0_Isoscalar3Pi::SetModelParameters(struct GeneralModel model) {
  m_model = model;
  // The line-shape registry is normally created by Hadron_Decay_Handler.  A
  // run with no hadronisation never constructs one, so it is NULL here and
  // every Get() would segfault; make it if it is missing.
  if (LineShapes==NULL) { LineShapes = new Line_Shapes(); LineShapes->Init(); }
  // This branch's FF_Parameters carries a parameter map; nothing here uses
  // it, so an empty one is passed.
  std::map<std::string,double> pmap;
  FF_Parameters params(ff_model::none,m_flavs,p_i,pmap,"FS_0_EE3Pi",&m_model);
  p_ff = new FF_0_Isoscalar3Pi(params);
}

DEFINE_CURRENT_GETTER(METOOLS::V_0_Isoscalar3Pi,"V_0_Isoscalar3Pi")

void ATOOLS::Getter<METOOLS::Current_Base,
		    METOOLS::ME_Parameters,METOOLS::V_0_Isoscalar3Pi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ \\gamma^* \\rightarrow \\pi^+ \\pi^- \\pi^0 $ \n\n"
    <<"Order: 0 = $\\pi^+$, 1 = $\\pi^-$, 2 = $\\pi^0$ \n\n"
    <<"\\[ \\epsilon^{\\mu\\alpha\\beta\\gamma} p_{+\\alpha} p_{-\\beta} "
    <<"p_{0\\gamma}\\, a(q^2) \\sum_{ij} BW_\\rho(s_{ij}) \\] \n\n"
    <<"Reference: Hoferichter, Hoid and Kubis, arXiv:1907.01556 \n"
    <<std::endl;
}
