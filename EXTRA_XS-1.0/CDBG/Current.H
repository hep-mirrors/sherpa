#ifndef Current_H
#define Current_H

#include "Vertex.H"

namespace EXTRAXS {

  template <typename CType> 
  class Current: public Current_Base {
  protected:

    std::vector<CType> m_j;

  public:

    inline Current(const ATOOLS::Flavour &fl):
      Current_Base(fl) {}

    // inline functions
    void ResetJ();
    void AddJ(const CType &j);

    const std::vector<CType> &J() const;
    
    void Evaluate();

    void Print() const;

  };// end of class Current

}// end of namespace EXTRAXS

#endif
