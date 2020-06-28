#include "ATOOLS/Phys/Weights.H"
#include "ATOOLS/Phys/Blob.H"

using namespace ATOOLS;

Weights::Weights(Variations_Type t, double w):
  type {t}
{
  size_t required_size = 1;
  if (type != Variations_Type::custom)
    required_size += s_variations->Size(t);
  weights.resize(required_size, w);
}

bool Weights::IsZero() const
{
  for (const auto& w : weights)
    if (w != 0.0)
      return false;
  return true;
}

std::string Weights::Name(size_t i) const
{
  if (i == 0) {
    return "Nominal";
  } else if (type == Variations_Type::custom) {
    return names[i];
  } else {
    return s_variations->GetVariationNameAt(i - 1, type);
  }
}

double& Weights::Nominal()
{
  return weights[0];
}

double Weights::Nominal() const
{
  return weights[0];
}

double& Weights::Variation(size_t i)
{
  assert(i + 1 < weights.size());
  return weights[i + 1];
}

double& Weights::operator[](const std::string& name)
{
  const auto it = std::find(names.begin(), names.end(), name);
  if (it == names.end()) {
    if (names.empty())
      names.push_back("Nominal");
    names.push_back(name);
    weights.push_back(weights.front());
    return weights.back();
  } else {
    return weights[it - names.begin()];
  }
}

Weights& Weights::operator=(double w)
{
  for (auto& weight : weights)
    weight = w;
  return *this;
}

Weights& Weights::operator*=(double rhs)
{
  for (auto& w : weights)
    w *= rhs;
  return *this;
}

Weights ATOOLS::operator*(Weights lhs, double rhs)
{
  lhs *= rhs;
  return lhs;
}

Weights ATOOLS::operator/(Weights lhs, double rhs)
{
  lhs /= rhs;
  return lhs;
}

Weights& Weights::operator*=(const Weights& rhs)
{
  if (this->IsUnnamedScalar()) {
    const auto w = weights[0];
    *this = rhs;
    *this *= w;
    return *this;
  }
  if (rhs.IsUnnamedScalar()) {
    *this *= rhs.weights[0];
    return *this;
  }
  assert(type == rhs.type);
  const auto size = weights.size();
  for (size_t i {0}; i < size; ++i) {
    weights[i] *= rhs.weights[i];
  }
  return *this;
}

Weights ATOOLS::operator*(Weights lhs, Weights rhs)
{
  lhs *= rhs;
  return lhs;
}

Weights& Weights::operator+=(const Weights& rhs)
{
  assert(type == rhs.type);
  const auto size = weights.size();
  for (size_t i {0}; i < size; ++i) {
    weights[i] += rhs.weights[i];
  }
  return *this;
}

Weights ATOOLS::operator+(Weights lhs, Weights rhs)
{
  lhs += rhs;
  return lhs;
}

Weights& Weights::operator-=(const Weights& rhs)
{
  assert(type == rhs.type);
  const auto size = weights.size();
  for (size_t i {0}; i < size; ++i) {
    weights[i] -= rhs.weights[i];
  }
  return *this;
}

Weights ATOOLS::operator-(Weights lhs, Weights rhs)
{
  lhs -= rhs;
  return lhs;
}

Weights Weights::operator-() const
{
  Weights ret = *this;
  const auto size = weights.size();
  for (size_t i {0}; i < size; ++i) {
    ret[i] = -weights[i];
  }
  return ret;
}

void ATOOLS::Reweight(Weights& w,
                      std::function<double(double, QCD_Variation_Params&)> f)
{
  w.type = Variations_Type::qcd;
  w.names.clear();
  const auto num_variation = s_variations->Size(Variations_Type::qcd);
  w.weights.resize(num_variation + 1,
                   w.weights.empty() ? 1.0 : w.weights.front());
  for (size_t i {1}; i < num_variation + 1; ++i) {
    w.weights[i] = f(w.weights[i], s_variations->Parameters(i - 1));
  }
}

void ATOOLS::Reweight(Weights& w,
                      std::function<double(double, size_t varindex, QCD_Variation_Params&)> f)
{
  w.type = Variations_Type::qcd;
  w.names.clear();
  const auto num_variation = s_variations->Size(Variations_Type::qcd);
  w.weights.resize(num_variation + 1,
                   w.weights.empty() ? 1.0 : w.weights.front());
  for (size_t i {1}; i < num_variation + 1; ++i) {
    w.weights[i] = f(w.weights[i], i - 1, s_variations->Parameters(i - 1));
  }
}

void ATOOLS::ReweightAll(Weights& w,
                      std::function<double(double, size_t varindex, QCD_Variation_Params*)> f)
{
  w.type = Variations_Type::qcd;
  w.names.clear();
  const auto num_variation = s_variations->Size(Variations_Type::qcd);
  w.weights.resize(num_variation + 1,
                   w.weights.empty() ? 1.0 : w.weights.front());
  for (size_t i {0}; i < num_variation + 1; ++i) {
    w.weights[i] = f(w.weights[i], i, (i == 0) ? nullptr : &s_variations->Parameters(i - 1));
  }
}

void ATOOLS::Reweight(Weights& w,
                      std::function<double(double, Qcut_Variation_Params&)> f)
{
  w.type = Variations_Type::qcut;
  w.names.clear();
  const auto num_variation = s_variations->Size(Variations_Type::qcd);
  w.weights.resize(num_variation + 1,
                   w.weights.empty() ? 1.0 : w.weights.front());
  for (size_t i {1}; i < num_variation + 1; ++i) {
    w.weights[i] = f(w.weights[i], s_variations->Qcut_Parameters(i - 1));
  }
}

void ATOOLS::ReweightAll(
    Weights& w,
    std::function<double(double, size_t varindex, Qcut_Variation_Params*)> f)
{
  w.type = Variations_Type::qcut;
  w.names.clear();
  const auto num_variation = s_variations->Size(Variations_Type::qcut);
  w.weights.resize(num_variation + 1,
                   w.weights.empty() ? 1.0 : w.weights.front());
  for (size_t i {0}; i < num_variation + 1; ++i) {
    w.weights[i] = f(w.weights[i], i, (i == 0) ? nullptr : &s_variations->Qcut_Parameters(i - 1));
  }
}

std::ostream& ATOOLS::operator<<(std::ostream& out, const Weights& w)
{
  for (size_t i {0}; i < w.weights.size(); ++i) {
    out << w.Name(i) << '=' << w.weights[i] << '\n';
  }
  return out;
}

bool Weights::IsUnnamedScalar() const
{
  return (weights.size() == 1 && type == Variations_Type::custom);
}

void Weights_Map::Clear()
{
  clear();
  base_weight = 1.0;
}

bool Weights_Map::HasVariations() const
{
  for (const auto& kv : *this)
    if (kv.second.HasVariations())
      return true;
  return false;
}

bool Weights_Map::IsZero() const
{
  if (base_weight == 0.0)
    return true;
  if (empty())
    return false; // empty is interpreted as a unity weight, i.e. 1.0
  for (const auto& kv : *this) {
    if (kv.second.IsZero())
      return true;
  }
  return false;
}

double Weights_Map::Nominal() const
{
  double w {base_weight};
  for (const auto& kv : *this) {
    w *= kv.second.Nominal();
  }
  return w;
}

double Weights_Map::NominalIgnoringVariationType(Variations_Type type) const
{
  double w {base_weight};
  for (const auto& kv : *this) {
    if (kv.second.type != type)
      w *= kv.second.Nominal();
  }
  return w;
}

Weights Weights_Map::Combine(Variations_Type type) const
{
  auto w = Weights {type};
  for (const auto& kv : *this) {
    if (kv.second.type == type)
      w *= kv.second;
  }
  return w;
}

Weights_Map& Weights_Map::operator*=(const Weights_Map& rhs)
{
  base_weight *= rhs.base_weight;
  for (const auto& kv : rhs) {
    auto it = find(kv.first);
    if (it != end()) {
      // both lhs and rhs have this key, hence we use the product
      it->second *= kv.second;
    } else {
      // lhs does not have this key, so we can just copy the values from rhs
      insert(kv);
    }
  }
  return *this;
}

Weights_Map& Weights_Map::operator+=(const Weights_Map& rhs)
{
  if (empty() && rhs.empty()) {
    base_weight += rhs.base_weight;
    return *this;
  }
  if (rhs.IsZero()) {
    return *this;
  }
  if (IsZero()) {
    *this = rhs;
    return *this;
  }

  // we treat the "ME" entries as absolute values, and everything else as
  // relative (pre)factors to those; therefore we can only add weight maps
  // that both have "ME" entries
  assert(this->find("ME") != this->end() && rhs.find("ME") != rhs.end());

  // make sure that lhs has all keys that are present in rhs, with default 1.0
  for (const auto& kv : rhs) {
    auto ret = this->insert(kv);
    if (ret.second)
      ret.first->second = 1.0;
  }

  // auxiliary construct in order to be able to iterate over {lhs, rhs} below
  std::array<const Weights_Map*, 2> operands = {this, &rhs};

  // now do the relative addition of lhs and rhs (ignoring "ME")
  for (auto& kv : *this) {
    if (kv.first == "ME")
      continue;
    const bool correlated_with_me_vars = (kv.second.type == (*this)["ME"].type);
    const auto num_weights = kv.second.weights.size();
    for (int i {0}; i < num_weights; ++i) {
      const int me_idx {correlated_with_me_vars ? i : 0};
      double numerator = 0.0, denominator = 0.0;
      for (auto op : operands) {
        const auto wgts_it = op->find(kv.first); // rhs might not have an entry
        const Weights& meweights = op->find("ME")->second;
        double w = 1.0;
        if (wgts_it != op->end()) {
          w = (wgts_it->second)[i];
        }
        const auto me = op->base_weight * meweights[me_idx];
        numerator += me * w;
        denominator += me;
      }
      kv.second[i] = numerator / denominator;
    }
  }

  // now add the "ME" entries
  (*this)["ME"] *= base_weight;
  base_weight = 1.0;
  (*this)["ME"] += rhs.base_weight * rhs.find("ME")->second;

  // finally, re-scale relative contributions to 1
  for (auto& kv : *this) {
    if (kv.first == "ME")
      continue;
    const auto relfac = kv.second.Nominal();
    if (relfac != 0.0) {
      kv.second /= relfac;
      (*this)["ME"].Nominal() *= relfac;
    }
  }

  return *this;
}

Weights_Map& Weights_Map::operator-=(const Weights_Map& rhs)
{
  Weights_Map negative_rhs = rhs;
  negative_rhs.base_weight = -negative_rhs.base_weight;
  return operator+=(negative_rhs);
}

std::ostream& ATOOLS::operator<<(std::ostream& out, const Weights_Map& w)
{
  out << w.base_weight << ":\n";
  for (const auto& e : w)
    out << e.first << "\n" << e.second << '\n';
  return out;
}

namespace ATOOLS {

  template <> Blob_Data<Weights_Map>::~Blob_Data() {}
  template class Blob_Data<Weights_Map>;
  template Weights_Map& Blob_Data_Base::Get<Weights_Map>();

}
