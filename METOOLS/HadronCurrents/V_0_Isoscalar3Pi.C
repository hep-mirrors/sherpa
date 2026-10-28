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
  // Exchange_ME calls this a second time with only the default v and a.
  // Those are not used here, and must not rebuild the form factor with
  // default settings.
  if (p_ff) return;
  // Merged into what the current already has, not replacing it. A later call
  // carrying only some keys must not wipe the settings an earlier one made:
  // Exchange_ME used to call every current a second time with only v and a,
  // and this current then rebuilt itself from those alone and silently fell
  // back to its default form factor.
  for (const auto & kv : model) m_model[kv.first] = kv.second;

  // Only the electromagnetic current is implemented: gamma* -> pi+ pi- pi0
  // through the isoscalar part of j_em. A Z couples to the same final state
  // through about -2 sin^2(theta_W) times this PLUS a strange-quark piece
  // (the phi), which is a form factor of its own; a W does not produce this
  // final state through a vector current at all. kf_none is let through: it
  // is what a contact term carries, and also a current built outside
  // Exchange_ME, which cannot be told apart from it here.
  const Flavour & bos = ExchangedBoson();
  if (bos.Kfcode()!=kf_none && bos.Kfcode()!=kf_photon)
    THROW(not_implemented,"V_0_Isoscalar3Pi implements the electromagnetic "
          "current only (gamma* -> pi+ pi- pi0), not exchange of "
          +ToString(bos)+". A Z would also need the strange-quark (phi) "
          "part of the current as its own form factor.");

  // The line-shape registry is normally created by Hadron_Decay_Handler.  A
  // run with no hadronisation never constructs one, so it is NULL here and
  // every Get() would segfault; make it if it is missing.
  if (LineShapes==NULL) { LineShapes = new Line_Shapes(); LineShapes->Init(); }
  // This branch's FF_Parameters carries a parameter map; nothing here uses
  // it, so an empty one is passed.
  std::map<std::string,double> pmap;
  // The parameterisation arrives as a FORM_FACTOR number, local to the
  // FF_0_Isoscalar3Pi family: 10, the default, is hep-ph/0512180 Eq. (10);
  // 11 is the Khuri-Treiman form. FORM_FACTOR is the key every other current
  // reads; EE3Pi_ff_model, the name this current used to read, is still
  // accepted so older cards keep working. Both, disagreeing, is an error
  // rather than a silent choice.
  const bool hasff(m_model.find("FORM_FACTOR")!=m_model.end());
  const bool hasee(m_model.find("EE3Pi_ff_model")!=m_model.end());
  const int ff(int(m_model("FORM_FACTOR",10.)+0.5));
  const int ee(int(m_model("EE3Pi_ff_model",10.)+0.5));
  if (hasff && hasee && ff!=ee)
    THROW(fatal_error,"V_0_Isoscalar3Pi: FORM_FACTOR = "+ToString(ff)+
          " and EE3Pi_ff_model = "+ToString(ee)+" disagree. Give one; "
          "FORM_FACTOR is the preferred name.");
  const int ffm(hasff ? ff : ee);
  FF_Parameters params(ffm,m_flavs,p_i,pmap,"FS_0_EE3Pi",&m_model);
  // Rebuilt from the merged settings; the previous one, if any, is freed
  // rather than leaked.
  if (p_ff) { delete p_ff; p_ff = NULL; }
  p_ff = new FF_0_Isoscalar3Pi(params);
}

DEFINE_CURRENT_GETTER(METOOLS::V_0_Isoscalar3Pi,"V_0_Isoscalar3Pi")

void ATOOLS::Getter<METOOLS::Current_Base,
		    METOOLS::ME_Parameters,METOOLS::V_0_Isoscalar3Pi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ \\gamma^* \\rightarrow \\pi^+ \\pi^- \\pi^0 $ \n\n"
    <<"Order: 0 = $\\pi^+$, 1 = $\\pi^-$, 2 = $\\pi^0$ \n\n"
    <<"\\[ \\epsilon^{\\mu\\alpha\\beta\\gamma} p_{+\\alpha} p_{-\\beta} "
    <<"p_{0\\gamma}\\, F(q^2, s_{ij}) \\] \n\n"
    <<"FORM\\_FACTOR = 10 (default): a sum of products of isoscalar Breit-Wigners "
    <<"in $q^2$ with rho line shapes in the two-pion sub-masses, "
    <<"hep-ph/0512180 Eq. (10). \n\n"
    <<"FORM\\_FACTOR = 11: the Khuri-Treiman form, a(q^2) times a tabulated "
    <<"dispersive basis, for $0.454 < \\sqrt{q^2} < 1.060$ GeV only; "
    <<"arXiv:2212.11767. \n\n"
    <<"EE3Pi\\_ff\\_model is accepted as an alias for FORM\\_FACTOR. \n"
    <<"Photon exchange only: a Z or W throws. \n"

    <<std::endl;
}
