#ifndef CSSHOWER_Showers_Shower_Enhancement_H
#define CSSHOWER_Showers_Shower_Enhancement_H

#include "ATOOLS/Math/Algebra_Interpreter.H"

namespace CSSHOWER {

  // VAR{expression}[|minimum|maximum], evaluated on a shower trial.
  class Shower_Enhancement : public ATOOLS::Tag_Replacer {
    ATOOLS::Algebra_Interpreter m_calc;
    double m_values[5];
    double m_min, m_max;
    bool m_hasbound;
  public:
    explicit Shower_Enhancement(const std::string &expression);
    double Evaluate(double z, double y, double x, double t, double Q2);
    bool HasBound() const { return m_hasbound; }
    double Bound() const { return m_max; }
    std::string ReplaceTags(std::string &expression) const;
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;
    void AssignId(ATOOLS::Term *term);
  };
}
#endif
