#ifndef DIRE__Shower__Gauge_H
#define DIRE__Shower__Gauge_H

#include "ATOOLS/Org/Getter_Function.H"
#include "DIRE/Tools/Kernel_Key.H"
#include "DIRE/Tools/Splitting.H"

namespace DIRE {

  class Kernel;
  class Vertex;

  class Gauge {
  protected:

    Kernel *p_sk;

    int m_type;

  public:

    Gauge(const Kernel_Key &k);

    virtual ~Gauge();

    double Value(const Splitting &s) const;
    double Estimate(const Splitting &s) const;

    virtual double K(const Splitting &s) const;
    virtual double KMax(const Splitting &s) const;

    virtual double Nf(const Splitting &s) const;
    virtual double Scale(const Splitting &s) const = 0;

    virtual double Weight(const Splitting &s) const = 0;
    virtual double Charge(const Splitting &s) const = 0;

    virtual double Coupling(const Splitting &s) const = 0;
    virtual double CplMax(const Splitting &s) const = 0;

    virtual bool GeneratePoint(Splitting &s) const = 0;
    virtual int Construct(Splitting &s) const = 0;
    
    virtual bool Allowed(const Splitting &s) const = 0;

    virtual void SetLimits() = 0;

    virtual double Solve(const double &cpl) const = 0;

  };// end of class Gauge

  typedef ATOOLS::Getter_Function
  <Gauge,Kernel_Key,std::less<std::string> > Gauge_Getter;

}// end of namepsace DIRE

#endif
