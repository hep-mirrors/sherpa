#ifndef MODEL_Main_Coupling_Data_H
#define MODEL_Main_Coupling_Data_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/STL_Tools.H"

#include <map>

namespace ATOOLS { class Function_Base; }

namespace MODEL {

  class Coupling_Data {
  protected:

    ATOOLS::Function_Base *p_cpl;

    std::string m_id;

    double m_def, m_fac, *p_scl;

  public:

    inline Coupling_Data(ATOOLS::Function_Base *const cpl,
			 const std::string &id,const double &def):
      p_cpl(cpl), m_id(id), m_def(def), m_fac(1.0), p_scl(NULL) {}

    void Calculate();

    inline void SetScale(double *const scl) { p_scl=scl; }

    inline const std::string &ID() const { return m_id; }

    inline double Default() const { return m_def; }
    inline double Factor() const  { return m_fac; }

  };// end of class Coupling_Data

  std::ostream &operator<<(std::ostream &str,const Coupling_Data &cd);

  struct Coupling_Map: public std::map
  <std::string,Coupling_Data*,ATOOLS::String_Sort> {
    void Calculate() const;
    Coupling_Data *Get(const std::string &tag) const;
  };

}// end of namespace MODEL

#endif
