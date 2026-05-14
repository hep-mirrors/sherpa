#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Phys/FormFactor.H"
#include "AddOns/Analysis/Tools/Particle_Qualifier.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Model_Base.H"

using namespace ATOOLS;

namespace METOOLS {

class FF_Pion : public Form_Factor {
private:
  double m_lambda, m_mp, m_muP;
  int m_mode;
  std::unique_ptr<ATOOLS::FormFactor> p_pionformfactor;

public:
  FF_Pion(const Vertex_Key &key)
      : Form_Factor("FF_Pion", key), m_mode(-1),
        p_pionformfactor(std::unique_ptr<FormFactor>(new FormFactor())) {
    p_pionformfactor->Init();
  }
  Complex FF() const override {
    if (p_pionformfactor->Type() == ATOOLS::ff::factored)
      return 1.;
    Current *j(m_mode < 0 ? p_v->JC() : p_v->J(m_mode));
    Vec4D Q;
    Complex fact = 1.0;
    int npions = 0;
    // msg_Out()<<"===================================="<<std::endl;
    for (size_t i(0); i < p_v->J().size(); ++i) {
      // PRINT_VAR(p_v->J(i)->Flav());
      if (p_v->J(i)->Flav().IsChPion()) {
        npions++;
        Q += p_v->J(i)->P();
        // PRINT_VAR(p_v->J(i)->P());
        // PRINT_VAR(p_v->J(i)->P().Mass());
      }
    }
    if (npions != 2)
      return 1.;
    double Q2(Q.Abs2());
    // if(Q.Mass()<0.3) return 1.;
    // if(p_v->J().size()==3){
    //   PRINT_VAR(j->P().Abs2());
    //   PRINT_VAR(Q2);
    // }

    // #ifdef DEBUG__BG
    // msg_Out() << METHOD << "(" << j->Id() << "," << j->Flav() << "): {\n"
    //                 << "  p = " << j->P() << " -> Q^2 = " << Q2 << "\n"
    //                 << "\n";
    // #endif
    if (p_v->J().size() == 3) {
      p_pionformfactor->Calculate(Q2);
      return p_pionformfactor->Form() * 0.5;
    }
    p_pionformfactor->Calculate(Q2);
    return p_pionformfactor->Form();
  }
}; // end of class FF_Pion

} // end of namespace METOOLS

using namespace METOOLS;

DECLARE_GETTER(FF_Pion, "FF_Pion", Form_Factor, Vertex_Key);
Form_Factor *Getter<Form_Factor, Vertex_Key, FF_Pion>::operator()(
    const Vertex_Key &args) const {
  return new FF_Pion(args);
}
void ATOOLS::Getter<Form_Factor, Vertex_Key, FF_Pion>::PrintInfo(
    std::ostream &str, const size_t width) const {
  str << "FF_Pion";
}

// DECLARE_GETTER(F2p,"F2p",Form_Factor,Vertex_Key);
// Form_Factor *Getter<Form_Factor,Vertex_Key,F2p>::
// operator()(const Vertex_Key &args) const
// { return new F2p(args); }
// void ATOOLS::Getter<Form_Factor,Vertex_Key,F2p>::
// PrintInfo(std::ostream &str,const size_t width) const
// { str<<"F2p"; }
