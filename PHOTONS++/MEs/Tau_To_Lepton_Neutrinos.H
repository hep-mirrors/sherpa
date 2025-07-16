#ifndef PHOTONS_MEs_Tau_To_Lepton_Neutrinos_H
#define PHOTONS_MEs_Tau_To_Lepton_Neutrinos_H

#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "PHOTONS++/Tools/Dipole_FI.H"

namespace PHOTONS {
  class Tau_To_Lepton_Neutrinos: public PHOTONS_ME_Base, public Dipole_FI {
    private:
      bool      m_switch;
      bool      m_fermi;

      virtual void    BoostOriginalPVVToMultipoleCMS();
      virtual void    FillMomentumArrays(const Particle_Vector_Vector&);


      virtual Complex InfraredSubtractedME_0_0();
      virtual Complex InfraredSubtractedME_0_1();
      virtual Complex InfraredSubtractedME_0_2();

      virtual Complex InfraredSubtractedME_1_05(unsigned int);
      virtual Complex InfraredSubtractedME_1_15(unsigned int);

      virtual Complex InfraredSubtractedME_2_1(unsigned int, unsigned int);

    public:
      Tau_To_Lepton_Neutrinos(const Particle_Vector_Vector&);
      virtual ~Tau_To_Lepton_Neutrinos();

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
