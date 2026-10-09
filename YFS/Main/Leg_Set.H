#ifndef YFS_Main_Leg_Set_H
#define YFS_Main_Leg_Set_H

#include <cstddef>

/*!
  \file Leg_Set.H

  Named operations on sets of legs (or photons) stored as bits: element i is
  in the set when bit i is set. The same encoding is used for

   - a helicity configuration: bit i set <=> leg i has its second helicity
     (for ACRAIC's Amplitude: helicity -1 of fermion leg i);
   - a flip mask: the legs whose helicity label is reversed between two
     conventions (e.g. Comix's and ACRAIC's, Acraic_Base::m_comixflip);
   - a subset of the event's photons (the NLO reals sum over them);
   - a set of external legs handed to Comix's propagator shifts.

*/

namespace YFS {

  //! The set {i}.
  template <typename T = std::size_t>
  constexpr T LegBit(std::size_t i) { return T(1) << i; }

  //! The number of subsets of n legs (= helicity configurations of n
  //! two-state legs): 2^n.
  template <typename T = std::size_t>
  constexpr T NumLegSubsets(std::size_t n) { return T(1) << n; }

  //! The set {0, ..., n-1}.
  template <typename T = std::size_t>
  constexpr T AllLegs(std::size_t n) { return NumLegSubsets<T>(n) - 1; }

  //! Whether leg i is in the set.
  template <typename T>
  constexpr bool HasLeg(T set, std::size_t i) { return (set >> i) & 1; }

  //! Leg i's bit of the set as 0 or 1 (for flipping one leg by another
  //! set's state: FlipLegs(cfg, LegBitOf(mask, i))).
  template <typename T>
  constexpr T LegBitOf(T set, std::size_t i) { return (set >> i) & 1; }

  //! The set with leg i added.
  template <typename T>
  constexpr T AddLeg(T set, std::size_t i) { return set | LegBit<T>(i); }

  //! The union of two sets.
  template <typename T>
  constexpr T LegUnion(T a, T b) { return a | b; }

  //! The configuration with the legs in flips reversed (symmetric
  //! difference); for flips inside set, the set without them.
  template <typename T>
  constexpr T FlipLegs(T set, T flips) { return set ^ flips; }

  //! The part of the set among legs 0 .. n-1.
  template <typename T>
  constexpr T LowLegs(T set, std::size_t n) { return set & AllLegs<T>(n); }

  //! The part of the set from leg n on, renumbered from 0.
  template <typename T>
  constexpr T HighLegs(T set, std::size_t n) { return set >> n; }

  //! low on legs 0 .. n-1 and high from leg n on: the inverse of LowLegs /
  //! HighLegs. Comix stores a helicity configuration this way, the n
  //! fermion legs first and the photons above them.
  template <typename T>
  constexpr T JoinLegs(T low, T high, std::size_t n) { return low | (high << n); }

  //! The number of legs in the set.
  template <typename T>
  constexpr std::size_t LegCount(T set)
  { std::size_t n(0); for (; set; set &= set - 1) ++n; return n; }

  //! The set without its smallest leg.
  template <typename T>
  constexpr T WithoutLowestLeg(T set) { return set & (set - 1); }

  //! The set without the legs of sub.
  template <typename T>
  constexpr T RemoveLegs(T set, T sub) { return set & ~sub; }

  //! The subset of `of` that follows sub when every non-empty subset is
  //! visited in decreasing order, 0 after the last:
  //!   for (T s(of); s; s = NextSubset(s, of)) { ... }
  template <typename T>
  constexpr T NextSubset(T sub, T of) { return (sub - 1) & of; }

  //! The smallest leg in a non-empty set.
  template <typename T>
  constexpr std::size_t LowestLeg(T set)
  { std::size_t i(0); for (; !(set & 1); set >>= 1) ++i; return i; }

  //! A helicity +-1 as a configuration bit: 0 for +1, 1 for -1.
  template <typename T = std::size_t>
  constexpr T HelicityBit(int helicity) { return helicity > 0 ? 0 : 1; }

  //! The helicity +-1 of leg i in configuration set.
  template <typename T>
  constexpr int LegHelicity(T set, std::size_t i) { return HasLeg(set, i) ? -1 : +1; }

}

#endif
