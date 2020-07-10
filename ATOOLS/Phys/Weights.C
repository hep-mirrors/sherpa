#include "ATOOLS/Phys/Weights.H"
#include "ATOOLS/Phys/Blob.H"

using namespace ATOOLS;

Weights& Weights::operator*=(const Weights& rhs)
{
  if (this->IsUnnamedScalar()) {
    const auto w = weights[0];
    *this = rhs;
    *this *= w;
    return *this;
  }
  if (rhs.IsUnnamedScalar()) {
    const auto w = weights[0];
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

Weights ATOOLS::MakeWeights(Variations_Type t)
{
  Weights w;
  w.type = t;
  if (t == Variations_Type::custom)
    return w;
  const auto num_variation = s_variations->Size(t);
  w.weights.resize(1 + num_variation, 1.0);
  return w;
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
        const auto wgts_it = op->find(kv.first);  // rhs might not have an entry
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
  if (empty() && rhs.empty()) {
    base_weight -= rhs.base_weight;
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
  // relative (pre)factors to those; therefore we can only subtract weight maps
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

  // now do the relative subtraction of lhs and rhs (ignoring "ME")
  for (auto& kv : *this) {
    if (kv.first == "ME")
      continue;
    const bool correlated_with_me_vars = (kv.second.type == (*this)["ME"].type);
    const auto num_weights = kv.second.weights.size();
    for (int i {0}; i < num_weights; ++i) {
      const int me_idx {correlated_with_me_vars ? i : 0};
      double numerator = 0.0, denominator = 0.0;
      for (auto op : operands) {
        const auto wgts_it = op->find(kv.first);  // rhs might not have an entry
        const Weights& meweights = op->find("ME")->second;
        double w = 1.0;
        if (wgts_it != op->end()) {
          w = (wgts_it->second)[i];
        }
        auto me = op->base_weight * meweights[me_idx];
        if (op != this)
          me = -me;
        numerator += me * w;
        denominator += me;
      }
      kv.second[i] = numerator / denominator;
    }
  }

  // now subtract the "ME" entries
  (*this)["ME"] *= base_weight;
  base_weight = 1.0;
  (*this)["ME"] -= rhs.base_weight * rhs.find("ME")->second;

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

namespace ATOOLS {

  template <> Blob_Data<Weights_Map>::~Blob_Data() {}
  template class Blob_Data<Weights_Map>;
  template Weights_Map& Blob_Data_Base::Get<Weights_Map>();

}
