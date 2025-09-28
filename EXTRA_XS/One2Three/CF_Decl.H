#ifndef EXTRA_XS_CF_DECL_H
#define EXTRA_XS_CF_DECL_H

#include "METOOLS/Explicit/Current.H"
#include "ATOOLS/Math/Vector.H"
#include <complex>
#include <vector>
#include <iostream>

namespace METOOLS {
  template <typename> class CSpinor;
  template <typename SType>
  class CF : public Current, public Current_Contractor<SType> {
  public:
    typedef std::complex<SType>   SComplex;
    typedef std::vector<SComplex> SComplex_Vector;
    typedef CSpinor<SType>        CSpinorType;

    void ConstructJ(const ATOOLS::Vec4D &p, const int ch,
                    const int cr, const int ca, const int mode);
    void SetGauge(const ATOOLS::Vec4D &k);

    inline void Print() const {
      const Current* c = static_cast<const Current*>(this);
      const_cast<Current*>(c)->Print();
    }

    inline void AddPropagator() {
      // fallback no-op to avoid undefined symbol if real implementation is not linked
    }

    void SContract(const Current &c, const Int_Vector &pols,
                   SComplex_Vector &ress, const size_t &offset) const;
    std::string Format(const CObject *c) const;
    char Type() const;
  };

  // removed: extern template class CF<double>;
}

#endif // EXTRA_XS_CF_DECL_H