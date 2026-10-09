#ifndef YFS_ACRAIC_Acraic_Amplitude_H
#define YFS_ACRAIC_Acraic_Amplitude_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MyComplex.H"
#include "YFS/Main/Leg_Set.H"

#include <vector>

namespace YFS {

  /*!
    A helicity amplitude over the external legs: one complex number per
    helicity configuration of the legs with two states (Leg_Set.H encoding).
  */
  class Amplitude{
  public:
    static constexpr int s_maxlegs  = 10;  //!< 2 -> 8; SetLegs refuses more
    static constexpr int s_warnlegs = 7;   //!< validated up to 2 -> 4 plus a Born photon
    static constexpr int s_maxhel  = NumLegSubsets<int>(s_maxlegs);
    static int s_nlegs;            //!< legs carrying two helicity states
    static inline int NHel() { return NumLegSubsets<int>(s_nlegs); }
    static void SetLegs(int n);
    //! leg i's bit in the packed index, or -1 for a leg with one spin state
    static std::vector<int> s_bit;
    //! packs only the legs with more than one spin state
    static void SetLegs(const ATOOLS::Flavour_Vector &flavs);
    static inline int Bit(size_t i)
    { return i < s_bit.size() ? s_bit[i] : -1; }
    //! four legs in states h0..h3; leg i's state is binary digit i
    static inline int Idx(int h0, int h1, int h2, int h3)
    { return h0 + 2*h1 + 4*h2 + 8*h3; }
    Amplitude();
    //! copy only the active entries (the capacity is not touched)
    Amplitude(const Amplitude &o);
    Amplitude &operator=(const Amplitude &o);
    ~Amplitude(){}
    Complex m_A[s_maxhel];
    Complex m_U[2][2];
    Complex m_V[2][2];
  };

}

#endif
