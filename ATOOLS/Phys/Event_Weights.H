#ifndef ATOOLS_Phys_Event_Weights_H
#define ATOOLS_Phys_Event_Weights_H

#include "ATOOLS/Phys/Variations.H"

#include <algorithm>
#include <cassert>

namespace ATOOLS {

  enum class Weight_Type {
    all = -1,
    nominal = 0
  };

  class Event_Weights {

  public:

    Event_Weights(double weight=1.0);
    Event_Weights(size_t num_variations, double weight=1.0);

    double Nominal() const;
    double& Nominal();
    double Variation(size_t i) const { return weights[1 + i]; }
    double& Variation(size_t i) { return weights[1 + i]; }
    double operator[](size_t i) const { return weights[i]; }
    double& operator[](size_t i) { return weights[i]; }
    bool ContainsVariations() const { return containsvariations; }
    bool Size() const { return weights.size(); }

    // explicit conversion, throws if more than one weight is stored
    explicit operator double() const;

    Event_Weights& operator=(double);
    Event_Weights& operator+=(const Event_Weights&);
    Event_Weights& operator-=(const Event_Weights&);
    Event_Weights& operator*=(const Event_Weights&);
    Event_Weights& operator/=(const Event_Weights&);
    Event_Weights& operator*=(double);
    template <class T> Event_Weights& operator*=(const std::vector<T>&);
    friend Event_Weights operator*(Event_Weights lhs, double rhs)
    {
      lhs *= rhs;
      return lhs;
    }
    friend Event_Weights operator*(Event_Weights lhs, Event_Weights rhs)
    {
      lhs *= rhs;
      return lhs;
    }
    Event_Weights& operator/=(double);

    friend bool operator==(const Event_Weights& lhs, double rhs)
    {
      if (std::all_of(lhs.weights.cbegin(),
                      lhs.weights.cend(),
                      [&rhs](double i) { return i == rhs; })) {
        return true;
      }
      return false;
    }
    friend bool operator!=(const Event_Weights& lhs, double rhs)
    {
      return !(lhs == rhs);
    }

    template <typename Handler>
    void Apply(Handler h);

    template <typename Handler>
    void ApplyAll(Handler h);

    friend std::ostream& operator<<(std::ostream&, const Event_Weights&);

  private:

    /// one nominal weight and one additional weight for each variation
    std::vector<double> weights;

    /// whether there is more than just a nominal weight
    bool containsvariations;

  };

  template <class T>
  Event_Weights& Event_Weights::operator*=(const std::vector<T>& rhs)
  {
    assert(rhs.size() == 1 || !ContainsVariations() ||
           weights.size() == rhs.size());
    if (!ContainsVariations()) {
      const auto weight = Nominal();
      weights.clear();
      weights.reserve(rhs.size());
      std::copy(rhs.begin(), rhs.end(), std::back_inserter(weights));
      containsvariations = rhs.size() > 1;
      (*this) *= weight;
    } else if (rhs.size() > 1) {
      for (int i {0}; i < weights.size(); ++i) {
        weights[i] *= rhs[i];
      }
    } else {
      (*this) *= rhs[0];
    }
    return *this;
  }

  template <typename Handler> void Event_Weights::Apply(Handler h)
  {
    for (int i {1}; i < weights.size(); ++i) {
      weights[i] = h(weights[i], i - 1, s_variations->Parameters(i - 1));
    }
  }

  template <typename Handler> void Event_Weights::ApplyAll(Handler h)
  {
    for (int i {0}; i < weights.size(); ++i) {
      weights[i] =
          h(weights[i], i, i == 0 ? nullptr : &s_variations->Parameters(i - 1));
    }
  }

} // namespace ATOOLS

#endif
