#include "ATOOLS/Phys/Event_Weights.H"
#include "ATOOLS/Phys/Blob_Data.C"

using namespace ATOOLS;

Event_Weights::Event_Weights(double weight)
    : Event_Weights(s_variations ? s_variations->Size() : 0, weight)
{
}

Event_Weights::Event_Weights(size_t num_variations, double weight)
    : weights(num_variations + 1, weight),
      containsvariations {num_variations != 0}
{
}

double Event_Weights::Nominal() const
{
  return weights[0];
}

double& Event_Weights::Nominal()
{
  return weights[0];
}

Event_Weights::operator double() const
{
  assert(weights.size() == 1);
  return weights[0];
}

Event_Weights& Event_Weights::operator=(double c)
{
  std::fill(weights.begin(), weights.end(), c);
  return *this;
}

Event_Weights& Event_Weights::operator+=(const Event_Weights& other)
{
  if (!ContainsVariations() && (*this)[0] == 0.0) {
    *this = other;
  } else if (other.ContainsVariations() || other[0] != 0.0) {
    assert(weights.size() == other.weights.size());
    for (int i {0}; i < weights.size(); ++i)
      weights[i] += other.weights[i];
  }
  return *this;
}

Event_Weights& Event_Weights::operator-=(const Event_Weights& other)
{
  assert(weights.size() == other.weights.size());
  for (int i {0}; i < weights.size(); ++i)
    weights[i] -= other.weights[i];
  return *this;
}

Event_Weights& Event_Weights::operator*=(const Event_Weights& other)
{
  (*this) *= other.weights;
  return *this;
}

Event_Weights& Event_Weights::operator/=(const Event_Weights& other)
{
  assert(!other.ContainsVariations() || weights.size() == other.weights.size());
  if (other.ContainsVariations()) {
    for (int i {0}; i < weights.size(); ++i) {
      weights[i] /= other.weights[i];
    }
  } else {
    (*this) /= other.Nominal();
  }
  return *this;
}

Event_Weights& Event_Weights::operator*=(double c)
{
  for (auto& w : weights)
    w *= c;
  return *this;
}

Event_Weights& Event_Weights::operator/=(double c)
{
  *this *= 1.0/c;
  return *this;
}

namespace ATOOLS {

  std::ostream& operator<<(std::ostream& out, const Event_Weights& w) {
    out << "nominal: " << w.weights[0] << "\n";
    out << "variations:";
    if (w.ContainsVariations()) {
      for (int i {1}; i < w.weights.size(); ++i) {
        out << "\n- " << w.weights[i];
      }
    } else {
      out << " none";
    }
    return out;
  }

  template class Blob_Data<Event_Weights>;
  template Event_Weights& Blob_Data_Base::Get<Event_Weights>();

}
