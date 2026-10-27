#include "CSSHOWER++/Showers/Shower_Enhancement.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MathTools.H"


using namespace CSSHOWER;
using namespace ATOOLS;

namespace {
  const char *variables[] = {"z", "y", "x", "t", "Q2"};
}

Shower_Enhancement::Shower_Enhancement(const std::string &expression)
  : m_values{.5,.5,.5,1.,1.},
    m_min(std::numeric_limits<double>::lowest()),
    m_max(std::numeric_limits<double>::max()), m_hasbound(false)
{
  // Split only after VAR{...}, so logical OR inside the expression is retained.
  const size_t close = expression.rfind('}');
  std::string formula;
  if (expression == "1") formula = "1";
  else {
    if (expression.find("VAR{") != 0 || close == std::string::npos)
      THROW(fatal_error, "SHOWER: ENHANCE_FUNCTION requires VAR{expression}[|minimum|maximum]");
    formula = expression.substr(4, close-4);
    if (close+1 < expression.size()) {
      const std::string suffix = expression.substr(close+1);
      if (suffix[0] != '|')
        THROW(fatal_error, "Invalid SHOWER: ENHANCE_FUNCTION range");
      std::stringstream stream(suffix.substr(1));
      std::string lo, hi, extra;
      if (!std::getline(stream, lo, '|') || !std::getline(stream, hi, '|') ||
          lo.empty() || hi.empty() || std::getline(stream, extra, '|') || suffix.back()=='|')
        THROW(fatal_error, "SHOWER: ENHANCE_FUNCTION requires both range limits");
      m_min = ToType<double>(lo);
      m_max = ToType<double>(hi);
      m_hasbound = true;
      if (IsBad(m_min) || IsBad(m_max) || m_min > m_max || m_min <= 0.)
        THROW(fatal_error, "Invalid SHOWER: ENHANCE_FUNCTION range: require 0 < minimum <= maximum");
    }
  }
  m_calc.SetTagReplacer(this);
  m_calc.AddTag("M_PI", ToString(M_PI,17));
  m_calc.AddTag("M_E", ToString(exp(1.),17));
  for (const auto variable : variables) m_calc.AddTag(variable, "0.5");
  m_calc.Interprete(formula);
  if (m_calc.Calculate()->Type()!='D')
    THROW(fatal_error, "SHOWER: ENHANCE_FUNCTION must evaluate to a real scalar");
}

std::string Shower_Enhancement::ReplaceTags(std::string &expression) const
{
  return m_calc.ReplaceTags(expression);
}

void Shower_Enhancement::AssignId(Term *term)
{
  for (size_t i=0; i<5; ++i)
    if (term->Tag() == variables[i]) { term->SetId(i); return; }
  if (term->Tag()=="M_PI" || term->Tag()=="M_E") { term->SetId(5); return; }
  THROW(fatal_error, "Unknown SHOWER: ENHANCE_FUNCTION variable: "+term->Tag());
}

Term *Shower_Enhancement::ReplaceTags(Term *term) const
{
  if (term->Id()<5) term->Set(m_values[term->Id()]);
  return term;
}

double Shower_Enhancement::Evaluate(double z, double y, double x, double t, double Q2)
{
  m_values[0]=z; m_values[1]=y; m_values[2]=x; m_values[3]=t; m_values[4]=Q2;
  const Term *result = m_calc.Calculate();
  if (result->Type()!='D')
    THROW(fatal_error, "SHOWER: ENHANCE_FUNCTION must evaluate to a real scalar");
  const double value = result->Get<double>();
  if (IsBad(value))
    THROW(fatal_error, "Non-finite SHOWER: ENHANCE_FUNCTION value");
  return std::max(m_min, std::min(m_max, value));
}
