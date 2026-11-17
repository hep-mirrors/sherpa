#ifndef METOOLS_Explicit_Form_Factor_H
#define METOOLS_Explicit_Form_Factor_H

#include "METOOLS/Explicit/Vertex_Key.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Math/MyComplex.H"

namespace METOOLS {
  
  class Vertex;

  class Form_Factor {
  protected:
    Vertex *p_v;
    std::string m_id;
  public:
    Form_Factor(const std::string &id,const Vertex_Key &key);
    virtual ~Form_Factor();

    virtual Complex FF() = 0;
    /*!
      The tag under which the model stores a hadron-photon vertex's form
      factor parameters, e.g. "_2212_2212_22" for Lambda2_2212_2212_22: the
      kf codes of the model vertex's legs with the hadrons first (in their
      vertex order) and the other legs after them. Comix builds the same
      vertex with its legs in other orders when the photon is not the
      current it computes (e+e- -> p pbar gamma), so the raw leg order does
      not name the parameters.
    */
    static std::string HadronFirstTag(const Vertex_Key &key);

    inline const std::string &ID() const { return m_id; }
  };// end of class FormFactor

  /*!
    A hadronic current linear in the momenta of its vertex's hadron legs,
      J^mu = sum_l c_l p_l^mu,
    with the coefficients c_l computed per event by the form factor and the
    Lorentz structure applied by the calculator (VPPPP_LC.C), so FF() is 1.
    Indexed by the declared leg of the model vertex (photon legs carry 0).
  */
  class Linear_Current_FF: public Form_Factor {
  public:
    Linear_Current_FF(const std::string &id,const Vertex_Key &key):
      Form_Factor(id,key) {}
    virtual const std::vector<Complex> &Coefficients() = 0;
    Complex FF() { return Complex(1.,0.); }
  };

  typedef ATOOLS::Getter_Function<Form_Factor,Vertex_Key,
				  std::less<std::string> > FF_Getter;

  std::ostream &operator<<(std::ostream &str,const Form_Factor &c);

  typedef std::vector<Form_Factor*> FF_Vector;
}// end of namespace METOOLS

#endif
