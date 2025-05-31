#ifndef MALARIC__Shower__Lorentz_FS_H
#define MALARIC__Shower__Lorentz_FS_H

#include "MALARIC/Shower/Lorentz.H"

namespace MALARIC {

  class Lorentz_FS_Rad: public Lorentz {
  protected:

    double Jacobian(const Splitting &s) const;

  public:

    Lorentz_FS_Rad(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Cluster(Splitting &s,const int mode,
		 const int i, const int j,const int k,
		 const ATOOLS::Cluster_Amplitude *ampl) const;

    double MEPSWeight(const Splitting &s) const;

  };// end of class Lorentz_FS_Rad

  class Lorentz_FS_Split: public Lorentz {
  protected:

    double Jacobian(const Splitting &s) const;

  public:

    Lorentz_FS_Split(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Cluster(Splitting &s,const int mode,
		 const int i, const int j,const int k,
		 const ATOOLS::Cluster_Amplitude *ampl) const;

    double MEPSWeight(const Splitting &s) const;

  };// end of class Lorentz_FS_Split

}// end of namepsace MALARIC

#endif
