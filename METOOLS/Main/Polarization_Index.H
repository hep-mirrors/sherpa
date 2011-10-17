#ifndef METOOLS__Main__Polarization_Index_H
#define METOOLS__Main__Polarization_Index_H

#include <vector>
#include <cstddef>

namespace METOOLS {

  class Polarization_Index {
  protected:

    std::vector<int> m_spins;

    size_t m_n;

  public:

    inline Polarization_Index(): m_n(0) {}
    inline Polarization_Index(const std::vector<int> &spins) { Init(spins); }

    void Init(const std::vector<int> &spins);

    size_t operator()(const std::vector<int> &spins) const;

    std::vector<int> operator()(size_t number) const;
    
    inline const std::vector<int> &Spins() const { return m_spins; }

    inline size_t N() const { return m_n; }

  };// end of class Polarization_Index

}// end of namespace METOOLS

#endif
