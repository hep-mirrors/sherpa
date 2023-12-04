#ifndef ALARIC__Shower__Lorentz_FS_H
#define ALARIC__Shower__Lorentz_FS_H

#include "ALARIC++/Shower/Lorentz.H"

namespace ALARIC {

  class Lorentz_FS_Rad: public Lorentz {
  protected:

    double Jacobian(const Splitting &s) const;

  public:

    Lorentz_FS_Rad(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Cluster(Splitting &s,const int mode,
		 const int i, const int j,const int k,
		 const ATOOLS::Cluster_Amplitude *ampl) const;

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

  };// end of class Lorentz_FS_Split

}// end of namepsace ALARIC

#endif
