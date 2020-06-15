#ifndef CFPSHOWER_Calculators_Kinematics_FF2_H
#define CFPSHOWER_Calculators_Kinematics_FF2_H

#include "CFPSHOWER++/Calculators/Kinematics_Base.H"

namespace CFPSHOWER {
  class Kinematics_FF2_Coll : public Kinematics_Base {
  protected:
    double CalculateY(Splitting & split);
    bool   ConstructSystem(Splitting & split);
    bool   KinCheck(Splitting & split);
  public:
    Kinematics_FF2_Coll();
    ~Kinematics_FF2_Coll();

    bool operator()(Splitting & split,Configuration & config);
    bool UpdateSystem(Splitting & split,Configuration & config);
    void CalculateWeight(Splitting & split);
  };

}

#endif
