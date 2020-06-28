#ifndef ATOOLS_Phys_Event_Weights_H
#define ATOOLS_Phys_Event_Weights_H

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Variations.H"

#include <algorithm>
#include <cassert>

namespace ATOOLS {

  class Weights_Map;

  class Weights {
  public:
    Weights(double w=1.0) : weights{w} {};
    Variations_Type Type() const { return type; }
    size_t Size() const { return weights.size(); }
    std::string Name(size_t i) const
    {
      if (type == Variations_Type::custom) {
        if (names.empty()) {
          return "Nominal";
        }
        return names[i];
      } else {
        if (i == 0)
          return "Nominal";
        return s_variations->GetVariationNameAt(i - 1, type);
      }
    }
    bool HasVariations() const { return weights.size() > 1; }
    Weights& operator=(double w)
    {
      if (weights.empty()) {
        weights.push_back(w);
        return *this;
      } else {
        for (auto& weight : weights)
          weight = w;
        return *this;
      }
    }
    Weights& operator*=(const Weights&);
    friend Weights operator*(Weights lhs, Weights rhs)
    {
      lhs *= rhs;
      return lhs;
    }
    Weights operator-()
    {
      Weights ret = *this;
      const auto size = weights.size();
      for (size_t i {0}; i < size; ++i) {
        ret[i] = -weights[i];
      }
      return ret;
    }
    Weights& operator+=(const Weights& rhs)
    {
      assert(type == rhs.type);
      const auto size = weights.size();
      for (size_t i {0}; i < size; ++i) {
        weights[i] += rhs.weights[i];
      }
      return *this;
    }
    friend Weights operator+(Weights lhs, Weights rhs)
    {
      lhs += rhs;
      return lhs;
    }
    Weights& operator-=(const Weights& rhs)
    {
      assert(type == rhs.type);
      const auto size = weights.size();
      for (size_t i {0}; i < size; ++i) {
        weights[i] -= rhs.weights[i];
      }
      return *this;
    }
    friend Weights operator-(Weights lhs, Weights rhs)
    {
      lhs -= rhs;
      return lhs;
    }
    template <class T> Weights& operator*=(const std::vector<T>& rhs)
    {
      const auto size = weights.size();
      // either we have n*1, or 1*n, or n*n
      assert(rhs.size() == 1 || weights.size() == 1 ||
             weights.size() == rhs.size());
      if (size == 1) {
        // 1*n
        assert(type != Variations_Type::custom);
        assert(type == Variations_Type::qcd
                   ? rhs.size() == s_variations->Size(Variations_Type::qcd)
                   : true);
        assert(type == Variations_Type::qcut
                   ? rhs.size() == s_variations->Size(Variations_Type::qcut)
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
    Weights& operator*=(double rhs)
    {
      for (auto& w : weights) {
        w *= rhs;
      }
      return *this;
    }
    friend Weights operator*(Weights lhs, double rhs)
    {
      lhs *= rhs;
      return lhs;
    }
    Weights& operator/=(double rhs) { return operator*=(1.0 / rhs); }
    friend Weights operator/(Weights lhs, double rhs)
    {
      lhs /= rhs;
      return lhs;
    }
    double& Nominal()
    {
      if (weights.empty()) {
        *this = 1.0;
      }
      return weights[0];
    }
    double Nominal() const
    {
      if (weights.empty()) {
        return 1.0;
      }
      if (names.empty()) {
        return weights[0];
      }
      if (names.front() != "Nominal") {
        return 1.0;
      }
      return weights[0];
    }
    bool IsZero() const
    {
      if (weights.empty()) {
        return false;
      }
      for (const auto& w : weights)
        if (w != 0.0)
          return false;
      return true;
    }
    double& operator[](size_t i) { return weights[i]; }
    double operator[](size_t i) const { return weights[i]; }
    double& operator[](const std::string& key)
    {
      auto it = std::find(names.begin(), names.end(), key);
      if (it == names.end()) {
        names.push_back(key);
        weights.push_back(weights.empty() ? 1.0 : weights.front());
        return weights.back();
      } else {
        return weights[it - names.begin()];
      }
    }
    double& var(size_t i) {
      if (i >= weights.size())
        THROW(fatal_error, "variation index " + ToString<int>(i) + " does not exist.");
      return weights[i + 1];
    }
    friend void Reweight(Weights&,
                         std::function<double(double, QCD_Variation_Params&)> f);
    friend void Reweight(Weights&,
                         std::function<double(double, size_t varindex, QCD_Variation_Params&)> f);
    friend void ReweightAll(Weights&,
                         std::function<double(double, size_t varindex, QCD_Variation_Params*)> f);
    friend void Reweight(Weights&,
                         std::function<double(double, Qcut_Variation_Params&)> f);
    friend void ReweightAll(Weights&,
                            std::function<double(double, size_t varindex, Qcut_Variation_Params*)> f);
    friend Weights MakeWeights(Variations_Type t);

    friend std::ostream& operator<<(std::ostream& out, const Weights& w)
    {
      for (size_t i {0}; i < w.weights.size(); ++i) {
        out << w.Name(i) << '=' << w.weights[i] << '\n';
      }
      return out;
    }

    friend Weights_Map;

  private:
    Variations_Type type {Variations_Type::custom};
    std::vector<double> weights;
    std::vector<std::string> names;
    bool IsUnnamedScalar() const {
      return (weights.size() == 1 && type == Variations_Type::custom &&
              (names.empty() || names[0] == "Nominal"));
    }
  };

  void Reweight(Weights& w,
                std::function<double(double, QCD_Variation_Params&)> f);
  void Reweight(Weights& w,
                std::function<double(double, size_t varindex, QCD_Variation_Params&)> f);
  void ReweightAll(Weights& w,
                   std::function<double(double, size_t varindex, QCD_Variation_Params*)> f);
  void Reweight(Weights& w,
                std::function<double(double, Qcut_Variation_Params&)> f);
  void ReweightAll(Weights&,
                   std::function<double(double, size_t varindex, Qcut_Variation_Params*)> f);
  Weights MakeWeights(Variations_Type t);

  class Weights_Map : public std::map<std::string, Weights> {
  public:
    Weights_Map(double w=1.0) : base_weight{w} {};
    explicit operator double() const {
      assert(empty() || (size() == 1 && !begin()->second.HasVariations()));
      return Nominal();
    }
    void Clear() {
      clear();
      base_weight = 1.0;
    }
    bool HasVariations() const {
      for (const auto& kv : *this)
        if (kv.second.HasVariations())
          return true;
      return false;
    }
    double Nominal() const
    {
      double w {base_weight};
      for (const auto& kv : *this) {
        w *= kv.second.Nominal();
      }
      return w;
    }
    double NominalIgnoringVariationType(Variations_Type type) const
    {
      double w {base_weight};
      for (const auto& kv : *this) {
        if (kv.second.type != type)
          w *= kv.second.Nominal();
      }
      return w;
    }
    Weights_Map& operator*=(const Weights_Map& rhs)
    {
      base_weight *= rhs.base_weight;
      for (const auto& kv : rhs) {
        auto it = find(kv.first);
        if (it != end()) {
          // both lhs and rhs have this key, hence we use the product
          it->second *= kv.second;
        } else {
          // lhs does not yet have this key, so we can just copy rhs' values
          insert(kv);
        }
      }
      return *this;
    }
    Weights_Map& operator+=(const Weights_Map&);
    Weights_Map& operator-=(const Weights_Map&);
    bool IsZero() const {
      if (base_weight == 0.0)
        return true;
      if (empty())
        return false;  // empty is interpreted as a trivial weight, i.e. 1.0
      for (const auto& kv : *this) {
        if (kv.second.IsZero())
          return true;
      }
      return false;
    }
    Weights Combine(Variations_Type type) const {
      auto w = MakeWeights(type);
      for (const auto& kv : *this) {
        if (kv.second.type == type)
          w *= kv.second;
      }
      return w;
    }
    friend std::ostream& operator<<(std::ostream& out, const Weights_Map& w)
    {
      out << w.base_weight << ":\n";
      for (const auto& e : w)
        out << e.first << "\n" << e.second << '\n';
      return out;
    }

  private:
    // make sure class users use our Clear() instead map's clear()
    using std::map<std::string, Weights>::clear;
    double base_weight {1.0};
  };

} // namespace ATOOLS

#endif
