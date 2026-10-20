
#include "METOOLS/HadronCurrents/V_Dirac_Dirac.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


V_Dirac_Dirac::V_Dirac_Dirac(
    const ATOOLS::Flavour_Vector& flavs,
    const std::vector<int>& indices,
    const std::string& name) :
  Current_Base(flavs, indices, name),
  m_cR(1.,0.), m_cL(1.,0.),
  p_f1(NULL)
{}


V_Dirac_Dirac::~V_Dirac_Dirac(){
  if (p_f1) {delete p_f1;p_f1 = NULL;}
}


void V_Dirac_Dirac::CheckTopology() const {
  if (p_i.size()!=2)
    THROW(fatal_error,"Current "+m_name+" needs exactly two legs, got "+ToString(p_i.size())+".");
}


void V_Dirac_Dirac::SetModelParameters(struct GeneralModel model) {
  m_cR = Complex(1.,0.);
  m_cL = Complex(1.,0.);

  map<string,double> pmap;
  FF_Parameters params(int(model[string("FORM_FACTOR")]),
            m_flavs,p_i,pmap,"F1_Dirac_Dirac",&model);
  p_f1 = FF_Getter::GetObject("FF_Dirac_Dirac_1",params);
}


void V_Dirac_Dirac::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{

  const Complex f1 = (*p_f1)(moms);

  XYZFunc F(moms, m_flavs, anti, p_i);
  for (int h0 = 0; h0 < 2; h0++) {
    for (int h1 = 0; h1 < 2; h1++) {
      // Equal left/right couplings give a pure vector current.
      const Vec4C amp = f1 * F.L(0, h0, 1, h1, m_cR, m_cL);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,h0));
      spins.push_back(make_pair(1,h1));
      Insert(amp, spins);
    }
  }
}


DEFINE_CURRENT_GETTER(METOOLS::V_Dirac_Dirac,"V_Dirac_Dirac")

void ATOOLS::Getter<
    METOOLS::Current_Base,
    METOOLS::ME_Parameters,
    METOOLS::V_Dirac_Dirac>::
PrintInfo(std::ostream &st, const size_t width) const
{
  st << "Vector fermion current with the Dirac form factor F1(Q^2).\n\n"
     << "Order: 0 = barred spinor, 1 = unbarred spinor.\n\n"
     << "Requires two proton legs or two neutron legs.\n\n"
     << "Form factor supplied by FF_Dirac_Dirac_1.\n"
     << std::endl;
}

