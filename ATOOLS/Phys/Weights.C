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
  nominals_prefactor = 1.0;
  is_absolute = false;
}

double Weights_Map::Get(const std::string& k, size_t i) const
{
  if (i == 0) {
    return Nominal();
  }
  auto it = find(k);
  if (it == end()) {
    return Nominal();
  }
  if (is_absolute) {
    return it->second[i];
  } else {
    return NominalIgnoringPrefactor() * it->second[i] / it->second[0];
  }
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
  if (base_weight == 0.0) {
    return true;
  }
  if (empty()) {
    return false; // empty is interpreted as a unity weight, i.e. 1.0
  }
  for (const auto& kv : *this) {
    if (kv.second.IsZero()) {
      return true;
    }
  }
  return false;
}

double Weights_Map::Nominal() const
{
  if (is_absolute) {
    if (empty()) {
      return base_weight;
    } else {
      return begin()->second.Nominal();
    }
  } else {
    double w {base_weight};
    for (const auto& kv : *this) {
      w *= kv.second.Nominal();
    }
    return nominals_prefactor * w;
  }
}

double Weights_Map::Nominal(const std::string& k) const
{
  const auto res = find(k);
  if (res == this->end())
    THROW(fatal_error, "Weights map does not have an entry for `" + k + "`.");
  return base_weight * res->second.Nominal();
}

double Weights_Map::NominalIgnoringVariationType(Variations_Type type) const
{
  assert(!is_absolute);
  double w {base_weight};
  for (const auto& kv : *this) {
    if (kv.second.type != type)
      w *= kv.second.Nominal();
  }
  return w;
}

Weights Weights_Map::RelativeValues(const std::string& k) const
{
  assert(!is_absolute);
  const auto it = this->find(k);
  if (it != this->end()) {
    auto ret = it->second;
    ret[0] *= nominals_prefactor;
    return ret;
  }
  return 1.0;
}

Weights Weights_Map::Combine(Variations_Type type) const
{
  assert(!is_absolute);
  auto w = Weights {type};
  for (const auto& kv : *this) {
    if (kv.second.type == type)
      w *= kv.second;
  }
  w[0] *= nominals_prefactor;
  return w;
}

void Weights_Map::SetZeroIfCloseToZero(double tolerance)
{
  MakeAbsolute();
  for (auto& kv : *this)
    for (auto& w : kv.second.weights)
      if (IsEqual(w, tolerance))
        w = 0.0;
  MakeRelative();
}

Weights_Map ATOOLS::operator*(Weights_Map lhs, double rhs)
{
  lhs *= rhs;
  return lhs;
}

Weights_Map ATOOLS::operator/(Weights_Map lhs, double rhs)
{
  lhs /= rhs;
  return lhs;
}

Weights_Map& Weights_Map::operator*=(const Weights_Map& rhs)
{
  assert(!is_absolute);
  base_weight *= rhs.base_weight;
  nominals_prefactor *= rhs.nominals_prefactor;
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

Weights_Map ATOOLS::operator*(Weights_Map lhs, const Weights_Map& rhs)
{
  lhs *= rhs;
  return lhs;
}

Weights_Map& Weights_Map::operator+=(const Weights_Map& rhs)
{
  assert(!is_absolute);
  if (empty() && rhs.empty()) {
    base_weight += rhs.base_weight;
    nominals_prefactor *= rhs.nominals_prefactor;
    return *this;
  }
  if (rhs.IsZero()) {
    return *this;
  }
  if (IsZero()) {
    *this = rhs;
    return *this;
  }

  // transform both sides into absolute storage instead of the default relative
  // storage
  MakeAbsolute();
  auto rhs_abs = rhs;
  rhs_abs.MakeAbsolute();

  // now addition is trivial
  for (auto& kv : rhs_abs) {
    (*this)[kv.first] += kv.second;
  }

  // transform back to relative storage
  MakeRelative();

  return *this;
}

Weights_Map& Weights_Map::operator-=(const Weights_Map& rhs)
{
  assert(!is_absolute);
  Weights_Map negative_rhs = rhs;
  negative_rhs.base_weight = -negative_rhs.base_weight;
  return operator+=(negative_rhs);
}

Weights_Map ATOOLS::operator+(Weights_Map lhs, const Weights_Map& rhs)
{
  lhs += rhs;
  return lhs;
}

Weights_Map ATOOLS::operator-(Weights_Map lhs, const Weights_Map& rhs)
{
  lhs -= rhs;
  return lhs;
}

void Weights_Map::MakeRelative()
{
  assert(is_absolute);

  // find any non-zero entry for normalisation, first check nominal entries
  double norm = 0.0;
  for (const auto& kv : *this) {
    norm = kv.second.Nominal();
    if (norm != 0.0) {
      break;
    }
  }
  if (norm == 0.0) {
    // all nominals are zero, we will reset them to 1.0 below and instead
    // account for the zero in the overall nominals prefactor
    nominals_prefactor = 0.0;
    for (const auto& kv : *this) {
      size_t num_vars = kv.second.Size() - 1;
      for (size_t i {0}; i < num_vars; ++i) {
        if (kv.second[i + 1] != 0.0) {
          // found a variation that is non-zero, we can use this as our new
          // overall normalisation factor
          norm = kv.second[i + 1];
          break;
        }
      }
    }
  } else {
    nominals_prefactor = 1.0;
  }
  if (norm == 0.0) {
    // Everything is zero, represent this with base_weight et to 0.0 and
    // everything else set to 1.0 and return.
    base_weight = 0.0;
    nominals_prefactor = 1.0;
    for (auto& kv : *this)
      kv.second = 1.0;
    is_absolute = false;
    return;
  }

  // apply normalisation
  for (auto& kv : *this) {
    kv.second /= norm;
  }
  base_weight = norm;

  // if all nominals are zero, we reset them here to 1.0, the zero gets stored
  // in the overall prefactor; this allows us to have non-zero variations even
  // when all nominals are actually zero (because those nominals get multiplied
  // as a prefactor to the variation)
  if (nominals_prefactor == 0.0) {
    for (auto& kv : *this) {
      kv.second.Nominal() = 1.0;
    }
  }

  is_absolute = false;
}

void Weights_Map::MakeAbsolute()
{
  assert(!is_absolute);

  // without variations, we only get rid of nominals_prefactor and return
  if (this->empty()) {
    base_weight *= nominals_prefactor;
    nominals_prefactor = 1.0;
    is_absolute = true;
    return;
  }

  // store all nominals
  std::map<std::string, double> nominals;
  for (const auto& kv : *this) {
    nominals[kv.first] = kv.second.Nominal();
  }

  // apply nominals of each Weights entry to all other Weights entries
  for (const auto& key_nom : nominals) {
    for (auto& kv : *this) {
      if (kv.first != key_nom.first) {
        kv.second *= key_nom.second;
      }
    }
  }

  // apply base_weight to all entries
  for (auto& kv : *this) {
    kv.second *= base_weight;
  }
  base_weight = 1.0;

  // apply nominals_prefactor to all nominals
  for (auto& kv : *this) {
    kv.second.Nominal() *= nominals_prefactor;
  }
  nominals_prefactor = 1.0;

  is_absolute = true;
}

#ifdef USING__MPI
void Weights_Map::MPI_Allreduce()
{
  int size=mpi->Size();
  if (size>1) {
    assert(!is_absolute);

    if (empty()) {
      mpi->Allreduce(&base_weight, 1, MPI_DOUBLE, MPI_SUM);
      mpi->Allreduce(&nominals_prefactor, 1, MPI_DOUBLE, MPI_PROD);
      return;
    }

    MakeAbsolute();
    for (auto& kv : *this) {
      const size_t n_wgts = kv.second.Size();
      mpi->Allreduce(&kv.second[0], n_wgts, MPI_DOUBLE, MPI_SUM);
    }
    MakeRelative();
  }
}
#endif

double Weights_Map::NominalIgnoringPrefactor() const
{
  if (is_absolute) {
    if (empty()) {
      return base_weight;
    } else {
      return begin()->second.Nominal();
    }
  } else {
    double w {base_weight};
    for (const auto& kv : *this) {
      w *= kv.second.Nominal();
    }
    return w;
  }
}

std::ostream& ATOOLS::operator<<(std::ostream& out, const Weights_Map& w)
{
  if (!w.is_absolute) {
    out << w.base_weight << " (nominals prefactor = " << w.nominals_prefactor
        << "):\n";
  }
  for (const auto& e : w) {
    out << e.first << "\n" << e.second << '\n';
  }
  return out;
}

namespace ATOOLS {

  Weights_Map sqrt(const Weights_Map& w)
  {
    auto root = w;
    root.base_weight = std::sqrt(root.base_weight);
    root.nominals_prefactor = std::sqrt(root.nominals_prefactor);
    for (auto& kv : root) {
      const size_t n_wgts = kv.second.Size();
      for (size_t i {0}; i < n_wgts; ++i) {
        kv.second[i] = std::sqrt(kv.second[i]);
      }
    }
    return root;
  }

  template <> Blob_Data<Weights_Map>::~Blob_Data() {}
  template class Blob_Data<Weights_Map>;
  template Weights_Map& Blob_Data_Base::Get<Weights_Map>();

}
