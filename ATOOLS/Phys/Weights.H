#ifndef ATOOLS_Phys_Event_Weights_H
#define ATOOLS_Phys_Event_Weights_H

#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Variations.H"

#include <algorithm>
#include <cassert>
#include <functional>

namespace ATOOLS {

  class Weights_Map;

  class Weights {

  public:

    Weights(double w = 1.0) : type {Variations_Type::custom}, weights {w} {};
    Weights(Variations_Type t, double w=1.0);

    explicit operator double() const { return Nominal(); }

    // general information getters
    Variations_Type Type() const { return type; }
    size_t Size() const { return weights.size(); }
    bool HasVariations() const { return weights.size() > 1; }
    bool IsZero() const;

    // content getters/setters to access individual weights, by index or name
    std::string Name(size_t, Variations_Source source = Variations_Source::all,
                     Variations_Name_Type name_type =
                         Variations_Name_Type::weight_name_convention) const;
    double& Nominal();
    double Nominal() const;
    double& Variation(size_t i);
    const double& Variation(size_t i) const;
    double& operator[](size_t i) { return weights[i]; }
    double operator[](size_t i) const { return weights[i]; }
    double& operator[](const std::string& name);

    /// Assign the same value to all weights.
    Weights& operator=(double);

    /// Multiply or divide all weights by the same value.
    Weights& operator*=(double);
    friend Weights operator*(Weights, double);
    Weights& operator/=(double rhs) { return operator*=(1.0 / rhs); }
    friend Weights operator/(Weights, double);

    /// Multiply by another set of weights weight-by-weight. This is only
    /// allowed if one set of weights has only a single entry, or both sets of
    /// weights have the same number of entries and the same variations type.
    Weights& operator*=(const Weights&);
    friend Weights operator*(Weights, Weights);

    /// Multiply by a vector of values element-by-element. This is only
    /// allowed if the weights or the vector have only a single entry, or
    /// they both have the same number of elements. Note that for managed
    /// variation types like QCD or Qcut variations, the number of vector
    /// elements must be either 1 or match the number of variations plus one
    /// (for the nominal entry).
    template <class T> Weights& operator*=(const std::vector<T>&);

    /// Add or subtract another set of weights weight-by-weight. This is only
    /// allowed if both sets of weights have the same number of entries and the
    /// same variations type.
    Weights& operator+=(const Weights&);
    friend Weights operator+(Weights, Weights);
    Weights& operator-=(const Weights&);
    friend Weights operator-(Weights, Weights);

    /// Create a copy with the sign of all weights swapped.
    Weights operator-() const;

    // Reweight functions for variation types managed by the Variations class,
    // that make it easy to iterate over all variations. The variants with the
    // suffix "All" also include the nominal entry in the iteration; as such
    // they pass the variation parameters as a pointer that is NULL for the
    // nominal entry.

    // QCD variations
    friend void Reweight(
        Weights&,
        std::function<double(double, QCD_Variation_Params&)>);
    friend void Reweight(
        Weights&,
        std::function<double(double, size_t varindex, QCD_Variation_Params&)>);
    friend void ReweightAll(
        Weights&,
        std::function<double(double, size_t varindex, QCD_Variation_Params*)>);

    // Qcut variations
    friend void Reweight(
        Weights&,
        std::function<double(double, Qcut_Variation_Params&)>);
    friend void Reweight(
        Weights&,
        std::function<double(double, size_t varindex, Qcut_Variation_Params&)>);
    friend void ReweightAll(
        Weights&,
        std::function<double(double, size_t varindex, Qcut_Variation_Params*)>);

    friend std::ostream& operator<<(std::ostream&, const Weights&);

    friend Weights_Map;

  private:

    /// Check if only a single value is present, and that no type has been set.
    bool IsUnnamedScalar() const;

    Variations_Type type;
    std::vector<double> weights;
    std::vector<std::string> names;
  };

  class Weights_Map : public std::map<std::string, Weights> {

  public:

    Weights_Map(double w=1.0) : base_weight{w} {};

    explicit operator double() const { return Nominal(); }

    /// Clears the content and sets the base weight to 1.0. The object is then
    /// equal to a newly default-constructed instance. Use this instead of the
    /// parent class member function `clear`.
    void Clear();

    // general information getters
    double Get(const std::string&, size_t) const;
    bool HasVariations() const;
    bool IsZero() const;
    double BaseWeight() const { return base_weight; }

    // access elements
    const Weights& at(const char* key) const {
      return this->std::map<std::string, Weights>::at(key);
    }
    const Weights& at(const std::string& key) const {
      return this->std::map<std::string, Weights>::at(key);
    }
    const Weights& at(Variations_Type type,
                      Variations_Source source = Variations_Source::all) const {
      if (type == Variations_Type::qcd) {
        return this->at(ToString(source));
      } else {
        return (*this).at("QCUT");
      }
    }

    // find elements
    const_iterator find(const char *key) const {
      return this->std::map<std::string, Weights>::find(key);
    }
    iterator find(const char *key) {
      return this->std::map<std::string, Weights>::find(key);
    }
    const_iterator find(const std::string& key) const {
      return this->std::map<std::string, Weights>::find(key);
    }
    iterator find(const std::string& key) {
      return this->std::map<std::string, Weights>::find(key);
    }
    const_iterator
    find(Variations_Type type,
         Variations_Source source = Variations_Source::all) const {
      if (type == Variations_Type::qcd) {
        return this->std::map<std::string, Weights>::find(ToString(source));
      } else {
        return this->std::map<std::string, Weights>::find("QCUT");
      }
    }

    // Convenience method to fill e.g. an ordinary map with variations; passing
    // Variations_Source::main only fills ME QCD-like variations (discarding
    // ones that originate from the shower), while passing
    // Variations_Source::all fill all variation (including custom one).
    template <typename T>
    void FillVariations(
        T &, Variations_Source source = ATOOLS::Variations_Source::all,
        Variations_Name_Type name_type =
            ATOOLS::Variations_Name_Type::weight_name_convention) const;

    // calculate nominal values, either including all entries, for a single
    // entry `k`, or excluding those that have a certain variation type
    double Nominal() const;
    double Nominal(const std::string& k) const;
    double NominalIgnoringVariationType(Variations_Type) const;

    // return relative Weights entries for a given name
    Weights& RelativeValues(const std::string& k) { return (*this)[k]; }
    Weights RelativeValues(const std::string& k) const;

    /// Set all values to exactly 0.0 that are close to zero. The tolerance
    /// is the maximum absolute distance from 0.0 defining this "closeness".
    /// This might be useful e.g. before calculating the square root.
    void SetZeroIfCloseToZero(double tolerance);

    /// Assign a value to the base weight
    Weights_Map& operator=(double w) { base_weight = w; return *this; }

    /// Multiply or divide the base weight
    Weights_Map& operator*=(double w) { base_weight *= w; return *this; }
    friend Weights_Map operator*(Weights_Map, double);
    Weights_Map& operator/=(double w) { base_weight /= w; return *this; }
    friend Weights_Map operator/(Weights_Map, double);

    /// Multiply by another weights map. If both maps share a key, the
    /// associated set of weights are multiplied with each other. If a key
    /// exists only in one of the maps, the associated set of weights will be
    /// present unchanged in the result. This is in accordance with the
    /// implicit assumption that an entry that does not exist explicitly has a
    /// weight of unit.
    Weights_Map& operator*=(const Weights_Map&);
    friend Weights_Map operator*(Weights_Map, const Weights_Map&);

    /// Add or subtract another weights map. This also makes the same implicit
    /// assumption that non-existing values are implicitly set to unity.
    Weights_Map& operator+=(const Weights_Map&);
    Weights_Map& operator-=(const Weights_Map&);
    friend Weights_Map operator+(Weights_Map, const Weights_Map&);
    friend Weights_Map operator-(Weights_Map, const Weights_Map&);

    friend Weights_Map sqrt(const Weights_Map&);

    friend std::ostream& operator<<(std::ostream&, const Weights_Map&);

    // generate variants where entries are absolute, or relative (the default);
    // this can be used to simplify operations dealing with individual entries
    // one-by-one, in particular when doing additions/subtractions
    void MakeRelative();
    void MakeAbsolute();

#ifdef USING__MPI
    // Using MPI Allreduce to combine weight maps from different MPI ranks;
    // this should not done too often, as it requires to synchronise the entire
    // structure of the weights maps, which in turn requires a nontrivial
    // amount of communication between the MPI ranks.
    void MPI_Allreduce();
#endif

  private:

    // make sure class users use our Clear() instead map's clear()
    using std::map<std::string, Weights>::clear;

    // same as Nominal(), but excluding the nominals_prefactor
    double NominalIgnoringPrefactor() const;

    /// While the Weights object in the map store relative factors due to
    /// specific variations or additional factors (like branching ratios), the
    /// base weight gives the overall scale of all weights. This is always
    /// taken as a prefactor when retrieving weights through Weights_Map member
    /// functions. A weights map that is implicitly constructed from a double
    /// (e.g. by returning 0.0 from a function that has Weights_Map as its
    /// return type) will have its base_weight set to this double, while it is
    /// otherwise empty.
    double base_weight;

    /// This is similar to the base weight, but is only applied to the nominal
    /// values of each Weights entry. This should usually be 1.0 or 0.0. While
    /// this seems redundant at first, it is not. Since all nominal values will
    /// always be applied as prefactors when calculating variations, it can be
    /// useful to store a zero in this member variable and set the actual
    /// nominal to a finite value like 1.0, even though it is formally zero. By
    /// this trick, we can have non-zero variations while all nominals are
    /// zero.
    double nominals_prefactor {1.0};

    bool is_absolute {false};
  };

  template <class T> Weights& Weights::operator*=(const std::vector<T>& rhs)
  {
    const auto size = weights.size();
    // either we have n*1, or 1*n, or n*n
    assert(rhs.size() == 1 || weights.size() == 1 ||
           weights.size() == rhs.size());
    if (size == 1) {
      // 1*n
      assert(type != Variations_Type::custom);
      assert(type == Variations_Type::qcd
                 ? rhs.size() == s_variations->Size(Variations_Type::qcd) + 1
                 : true);
      assert(type == Variations_Type::qcut
                 ? rhs.size() == s_variations->Size(Variations_Type::qcut) + 1
                 : true);
      const auto weight = Nominal();
      weights.clear();
      weights.reserve(rhs.size());
      std::copy(rhs.begin(), rhs.end(), std::back_inserter(weights));
      (*this) *= weight;
    } else if (rhs.size() > 1) {
      // n*n
      for (int i {0}; i < weights.size(); ++i) {
        weights[i] *= rhs[i];
      }
    } else {
      // n*1
      (*this) *= rhs[0];
    }
    return *this;
  }

  Weights operator*(Weights, double);
  Weights operator/(Weights, double);
  Weights operator*(Weights, Weights);
  Weights operator+(Weights, Weights);
  Weights operator-(Weights, Weights);

  std::ostream& operator<<(std::ostream&, const Weights&);

  void Reweight(
      Weights&,
      std::function<double(double, QCD_Variation_Params&)>);
  void Reweight(
      Weights&,
      std::function<double(double, size_t varindex, QCD_Variation_Params&)>);
  void ReweightAll(
      Weights&,
      std::function<double(double, size_t varindex, QCD_Variation_Params*)>);
  void Reweight(
      Weights&,
      std::function<double(double, Qcut_Variation_Params&)>);
  void Reweight(
      Weights&,
      std::function<double(double, size_t varindex, Qcut_Variation_Params&)>);
  void ReweightAll(
      Weights&,
      std::function<double(double, size_t varindex, Qcut_Variation_Params*)>);

  template <typename T>
  void Weights_Map::FillVariations(T &map, Variations_Source source,
                                   Variations_Name_Type name_type) const {
    DEBUG_FUNC(source);

    // Fill "managed" QCD variations
    for (const auto type : s_variations->ManagedVariationTypes()) {
      // calculate contributions
      const_iterator it {this->find(type, source)};
      if (it != this->end()) {
        const Weights& weights {it->second};
        const double relfac {NominalIgnoringVariationType(type)};
        // do remaining combination and store resulting weights
        size_t num_vars = weights.Size() - 1;
        for (size_t i(0); i < num_vars; ++i) {
          const std::string varname {weights.Name(i + 1, source, name_type)};
          map[varname] = weights.Variation(i) * relfac;
          msg_Debugging() << varname << ": " << weights.Variation(i) * relfac
                          << '\n';
        }
      }
    }

    // Fill other "custom" variations
    if (source == Variations_Source::all) {
      for (const auto& kv : *this) {
        if (kv.second.Type() != Variations_Type::custom)
          continue;
        const auto size = kv.second.Size();
        // Iterate weights, skipping the nominal entry (at position 0).
        for (size_t i {1}; i < size; ++i) {
          map[kv.first + "." + kv.second.Name(i)] =
              Nominal() / kv.second.Nominal() * kv.second[i];
        }
      }
    }
  }

  std::ostream& operator<<(std::ostream& out, const Weights_Map& w);
  Weights_Map operator+(Weights_Map lhs, const Weights_Map& rhs);
  Weights_Map operator-(Weights_Map lhs, const Weights_Map& rhs);
  Weights_Map operator*(Weights_Map lhs, const Weights_Map& rhs);
  Weights_Map operator*(Weights_Map lhs, double rhs);
  Weights_Map operator/(Weights_Map lhs, const Weights_Map& rhs);
  Weights_Map operator/(Weights_Map lhs, double rhs);

} // namespace ATOOLS

#endif
