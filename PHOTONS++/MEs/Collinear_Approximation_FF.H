#ifndef PHOTONS_MEs_Collinear_Approximation_FF_H
#define PHOTONS_MEs_Collinear_Approximation_FF_H

#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "PHOTONS++/Tools/Dipole_FF.H"
#include "METOOLS/Loops/Divergence_Array.H"

namespace PHOTONS {
  class Collinear_Approximation_FF: public PHOTONS_ME_Base, public Dipole_FF {
    private:
    int m_switch;
    virtual void    BoostOriginalPVVToMultipoleCMS();
    virtual void    FillMomentumArrays(const Particle_Vector_Vector&);


      virtual Complex InfraredSubtractedME_0_0();
      virtual Complex InfraredSubtractedME_0_1(const int& = 1);
      virtual Complex InfraredSubtractedME_0_2();

      virtual Complex InfraredSubtractedME_1_05(unsigned int, const double& = 1.);
      virtual Complex InfraredSubtractedME_1_15(unsigned int, const double& = 1.);

      virtual Complex InfraredSubtractedME_2_1(unsigned int, unsigned int);

    public:
      Collinear_Approximation_FF(const Particle_Vector_Vector&);
      virtual ~Collinear_Approximation_FF();

      virtual double GetBeta_0_0();
      virtual double GetBeta_0_1(const int& = 1);
      virtual double GetBeta_0_2();

      virtual double GetBeta_1_1(unsigned int);
      virtual double GetBeta_1_2(unsigned int);

      virtual double GetBeta_2_2(unsigned int, unsigned int);

      virtual double SmodFull(unsigned int);
      virtual double Smod(unsigned int);
      virtual double Smod(unsigned int, unsigned int, unsigned int);
      double Dmod(unsigned int, unsigned int, unsigned int);
      double Kallen(double, double, double);
      virtual void Print_Info();
  };
}
#endif
