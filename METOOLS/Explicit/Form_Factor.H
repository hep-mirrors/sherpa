#ifndef METOOLS_Explicit_Form_Factor_H
#define METOOLS_Explicit_Form_Factor_H

#include "METOOLS/Explicit/Vertex_Key.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace METOOLS {

  class Vertex;

  class Form_Factor {
  protected:

    Vertex *p_v;

    std::string m_id;

  public:

    Form_Factor(const std::string &id,const Vertex_Key &key);

    virtual ~Form_Factor();

    virtual double FF(double) const = 0;

    inline const std::string &ID() const { return m_id; }    

  };// end of class FormFactor
  
  typedef ATOOLS::Getter_Function<Form_Factor,Vertex_Key,
				  std::less<std::string> > FF_Getter;

  std::ostream &operator<<(std::ostream &str,const Form_Factor &c);

  typedef std::vector<Form_Factor*> FF_Vector;

}// end of namespace METOOLS

#endif
