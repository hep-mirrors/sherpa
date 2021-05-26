#ifndef ATOOLS_Phys_Event_Weights_H
#define ATOOLS_Phys_Event_Weights_H

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Variations.H"

#include <algorithm>
#include <cassert>
#include <valarray>

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
    std::string Name(size_t) const;
    double& Nominal();
    double Nominal() const;
    double& Variation(size_t i);
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

    friend std::ostream& operator<<(std::ostream& out, const Weights& w);

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
    bool HasVariations() const;
    bool IsZero() const;

    // calculate nominal values, either including all entries or excluding
    // those that have a certain variation type
    double Nominal() const;
    double NominalIgnoringVariationType(Variations_Type) const;

    /// Calculate the product of all entries that share the same type, e.g. to
    /// get the combined QCD or Qcut variations from different sources. Note
    /// that this returns a Weights object, that contains the nominal product
    /// and the entry-by-entry product of all variations.
    Weights Combine(Variations_Type type) const;

    /// Multiply by another weights map. If both maps share a key, the
    /// associated set of weights are multiplied with each other. If a key
    /// exists only in one of the maps, the associated set of weights will be
    /// present unchanged in the result. This is in accordance with the
    /// implicit assumption that an entry that does not exist explicitly has a
    /// weight of unit.
    Weights_Map& operator*=(const Weights_Map&);

    /// Add or subtract another weights map. This makes the same implicit
    /// assumption as stated for multiplication above. Also note that only the
    /// entry with the key "ME" are considered to contain absolute values,
    /// whereas all non-"ME" entries are treated as relative values (i.e. they
    /// are added together weighted by the nominal value in "ME"). To make this
    /// behaviour as explicit as possible and to prevent suprises, an "ME"
    /// entry must be present in both maps.
    Weights_Map& operator+=(const Weights_Map&);
    Weights_Map& operator-=(const Weights_Map&);

    friend std::ostream& operator<<(std::ostream&, const Weights_Map&);

  private:
    // make sure class users use our Clear() instead map's clear()
    using std::map<std::string, Weights>::clear;

    /// This member variable is mainly used to allow functions that return
    /// Weights_Map instances to just return 0.0, which is then used to
    /// implicitly construct an instance with `base_weight = 0.0` via the
    /// default constructor. Otherwise it should usually be 1.0.
    double base_weight;
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
} // namespace ATOOLS

#endif
