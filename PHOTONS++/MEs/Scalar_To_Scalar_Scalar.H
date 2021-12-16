#ifndef PHOTONS_MEs_Scalar_To_Scalar_Scalar_H
#define PHOTONS_MEs_Scalar_To_Scalar_Scalar_H

#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "PHOTONS++/Tools/Dipole_FF.H"

namespace PHOTONS {
  class Scalar_To_Scalar_Scalar: public PHOTONS_ME_Base, public Dipole_FF {
    private:
      bool      m_switch;
      Complex   m_Gamma;

      virtual void    BoostOriginalPVVToMultipoleCMS();

      virtual Complex InfraredSubtractedME_0_0();
      virtual Complex InfraredSubtractedME_0_1(const int& = 1);
      virtual Complex InfraredSubtractedME_0_2();

      virtual Complex InfraredSubtractedME_1_05(unsigned int, const double& = 1.);
      virtual Complex InfraredSubtractedME_1_15(unsigned int, const double& = 1.);

      virtual Complex InfraredSubtractedME_2_1(unsigned int, unsigned int);

    public:
      Scalar_To_Scalar_Scalar(const Particle_Vector_Vector&);
      virtual ~Scalar_To_Scalar_Scalar();

      virtual double GetBeta_0_0();
      virtual double GetBeta_0_1(const int& = 1);
      virtual double GetBeta_0_2();

      virtual double GetBeta_1_1(unsigned int);
      virtual double GetBeta_1_2(unsigned int);

      virtual double GetBeta_2_2(unsigned int, unsigned int);

      virtual void   FillMomentumArrays(const Particle_Vector_Vector&);
      virtual double SmodFull(unsigned int);
      virtual double Smod(unsigned int);
      virtual double Smod(unsigned int, unsigned int, unsigned int);
      virtual void Print_Info();
  };
}
#endif
