#ifndef Extra_XS_h
#define Extra_XS_h

#include "Flavour.H"
#include "Vector.H"

namespace AMEGIC {
  class Extra_XS {
    double  s,MZ2,DZZ,ReDZZ,P_PL,P_MI;
    int     xs;
    ATOOLS::Flavour flout[2];
    double Lambda(double,double,double);

    double Pqr(double,double);
    double DiffBhabhaScattering(ATOOLS::Vec4D*);
    double DiffMollerScattering(ATOOLS::Vec4D*);
    double NeutralinoPairProduction();
    double SelectronPairProduction();
    double CharginoPairProduction();
    double HiggsPairProduction();
    double HiggsBremsstrahlung(); 
    double WPairProduction();
    double SneutrinoPairProduction();
    double ZPairProduction();   
    double sFermionPairProduction();
 public: 
    void Init(int,int,ATOOLS::Flavour*,ATOOLS::Flavour*);
    double Calculate();
    
    void   InitDiff(int,int,ATOOLS::Flavour*,ATOOLS::Flavour*);
    double Differential(ATOOLS::Vec4D*);
  };
}
#endif
