#ifndef PHOTONS_MEs_Vector_To_Scalar_Scalar_H
#define PHOTONS_MEs_Vector_To_Scalar_Scalar_H

#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "PHOTONS++/Tools/Dipole_FF.H"

using namespace ATOOLS;

namespace PHOTONS {
  class Vector_To_Scalar_Scalar: public PHOTONS_ME_Base, public Dipole_FF {
    private:
      bool      m_switch;
      double    m_Gamma;

      virtual void    BoostOriginalPVVToMultipoleCMS();
      virtual void    FillMomentumArrays(const Particle_Vector_Vector&);


      virtual Complex InfraredSubtractedME_0_0();
      virtual Complex InfraredSubtractedME_0_1();
      virtual Complex InfraredSubtractedME_0_2();

      virtual Complex InfraredSubtractedME_1_05(unsigned int);
      virtual Complex InfraredSubtractedME_1_15(unsigned int);

      virtual Complex InfraredSubtractedME_2_1(unsigned int, unsigned int);

    public:
      Vector_To_Scalar_Scalar(const Particle_Vector_Vector&);
      virtual ~Vector_To_Scalar_Scalar();

      virtual double GetBeta_0_0();
      virtual double GetBeta_0_1();
      virtual double GetBeta_0_2();

      virtual double GetBeta_1_1(unsigned int);
      virtual double GetBeta_1_2(unsigned int);

      virtual double GetBeta_2_2(unsigned int, unsigned int);

      virtual double Smod(unsigned int);
  };
}
#endif
