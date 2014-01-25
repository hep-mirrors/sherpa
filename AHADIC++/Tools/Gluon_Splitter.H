#ifndef AHADIC_Tools_Gluon_Splitter_H
#define AHADIC_Tools_Gluon_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"
#include "AHADIC++/Tools/Dipole.H"

namespace AHADIC {
  class Gluon_Splitter : public Splitter_Base {
  private:
    double m_etay, m_etay_lead;
    Proto_Particle * p_out1, * p_out2;

    bool   SelectSplitter(Proto_Particle * part1,Proto_Particle * part2);
    bool   PoppedMassPossible(const double & m2);
    bool   ConstructSystem();
    void   ConstructKinematics(const double & etay);  
    bool   AcceptSystem(const double & pt2max);
    void   MakeKinematics();
    void   MakeParticles();
    double FixExponent();
  public:
    Gluon_Splitter();

    bool operator()(Dipole * dip,const bool & first,const bool & vetodq=false);
    void GetNewParticles(Proto_Particle *& out1,Proto_Particle *& out2) const;
  };

  inline void Gluon_Splitter::
  GetNewParticles(Proto_Particle *& out1,Proto_Particle *& out2) const {
    out1 = p_out1; out2 = p_out2;
  }
}

#endif
