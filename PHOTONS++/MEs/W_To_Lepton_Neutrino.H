#ifndef PHOTONS_MEs_W_To_Lepton_Neutrino_H
#define PHOTONS_MEs_W_To_Lepton_Neutrino_H

#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "PHOTONS++/Tools/Dipole_FI.H"

namespace PHOTONS {
  class W_To_Lepton_Neutrino: public PHOTONS_ME_Base, public Dipole_FI {
    private:
      bool      m_switch;

      virtual void    BoostOriginalPVVToMultipoleCMS();
      virtual void    FillMomentumArrays(const Particle_Vector_Vector&);


      virtual Complex InfraredSubtractedME_0_0();
      virtual Complex InfraredSubtractedME_0_1();
      virtual Complex InfraredSubtractedME_0_2();

      virtual Complex InfraredSubtractedME_1_05(unsigned int);
      virtual Complex InfraredSubtractedME_1_15(unsigned int);

      virtual Complex InfraredSubtractedME_2_1(unsigned int, unsigned int);

    public:
      W_To_Lepton_Neutrino(const Particle_Vector_Vector&);
      virtual ~W_To_Lepton_Neutrino();

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
