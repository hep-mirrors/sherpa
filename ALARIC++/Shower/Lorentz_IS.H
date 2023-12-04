#ifndef ALARIC__Shower__Lorentz_IS_H
#define ALARIC__Shower__Lorentz_IS_H

#include "ALARIC++/Shower/Lorentz.H"

namespace ALARIC {

  class Lorentz_IS: public Lorentz {
  protected:

    double Jacobian(const Splitting &s) const;

    double PDFEstimate(const Splitting &s) const;

  public:

    Lorentz_IS(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Cluster(Splitting &s,const int mode,
		 const int i, const int j,const int k,
		 const ATOOLS::Cluster_Amplitude *ampl) const;

  };// end of class Lorentz_IS

}// end of namepsace ALARIC

#endif
