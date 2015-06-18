#include "DIRE/Shower/Gauge.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE DIRE::Kernel_Key
#define OBJECT_TYPE DIRE::Gauge
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

using namespace DIRE;

Gauge::Gauge(const Kernel_Key &k):
  p_sk(k.p_k), m_type(k.m_type)
{
}

Gauge::~Gauge()
{
}

double Gauge::Value(const Splitting &s) const
{
  return Coupling(s)*Weight(s);
}

double Gauge::Estimate(const Splitting &s) const
{
  return CplMax(s)*Charge(s);
}

double Gauge::K(const Splitting &s) const
{
  THROW(not_implemented,"Invalid call");
}

double Gauge::KMax(const Splitting &s) const
{
  THROW(not_implemented,"Invalid call");
}

double Gauge::Nf(const Splitting &s) const
{
  THROW(not_implemented,"Invalid call");
}
