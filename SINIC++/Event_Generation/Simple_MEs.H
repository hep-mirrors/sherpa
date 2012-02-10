#ifndef SINIC_Main_Simple_MEs_H
#define SINIC_Main_Simple_MEs_H

#include "ATOOLS/Phys/Flavour.H"

namespace SINIC {
  class Simple_MEs {
  private:
    double q1q2_q1q2(const double & hats,const double & hatt,const double & hatu,
		     const double & mu2=0.);
    double qq_qq(const double & hats,const double & hatt,const double & hatu,
		 const double & mu2=0.);
    double qqb_qqb(const double & hats,const double & hatt,const double & hatu,
		   const double & mu2=0.);
    double gq_gq(const double & hats,const double & hatt,const double & hatu,
		 const double & mu2=0.);
    double gg_gg(const double & hats,const double & hatt,const double & hatu,
		 const double & mu2=0.);
    double gg_qqb(const double & hats,const double & hatt,const double & hatu,
		  const double & mu2=0.);
  public:
    Simple_MEs();
    ~Simple_MEs();
    
    double operator()(const ATOOLS::Flavour & in1,const ATOOLS::Flavour & in2,
		      const ATOOLS::Flavour & out1,const ATOOLS::Flavour & out2,
		      const double & hats,const double & hatt,const double & hatu,
		      const double & mu2=0.);

  };
}

#endif
