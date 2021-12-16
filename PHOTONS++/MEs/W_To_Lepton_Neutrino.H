#ifndef PHOTONS_MEs_W_To_Lepton_Neutrino_H
#define PHOTONS_MEs_W_To_Lepton_Neutrino_H

#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "PHOTONS++/Tools/Dipole_FI.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "AddOns/OpenLoops/OpenLoops_Interface.H"
#include "AddOns/OpenLoops/OpenLoops_Virtual.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"


namespace PHOTONS {
  class W_To_Lepton_Neutrino: public PHOTONS_ME_Base, public Dipole_FI {
    private:
      bool      m_switch;
      METOOLS::DivArrC m_res, m2_loop;
      ATOOLS::Vec4D      m_csmoms[2][9][9];

      virtual void    BoostOriginalPVVToMultipoleCMS();
      virtual void    FillMomentumArrays(const Particle_Vector_Vector&);


      virtual Complex InfraredSubtractedME_0_0();
      virtual Complex InfraredSubtractedME_0_1(const int& = 1);
      virtual Complex InfraredSubtractedME_0_2();

      virtual Complex InfraredSubtractedME_1_05(unsigned int, const double& = 1.);
      virtual Complex InfraredSubtractedME_1_15(unsigned int, const double& = 1.);

      virtual Complex InfraredSubtractedME_2_1(unsigned int, unsigned int);

    public:
      W_To_Lepton_Neutrino(const Particle_Vector_Vector&);
      virtual ~W_To_Lepton_Neutrino();

      virtual double GetBeta_0_0();
      virtual double GetBeta_0_1(const int& = 1);
      virtual double GetBeta_0_2();

      virtual double GetBeta_1_1(unsigned int);
      virtual double GetBeta_1_2(unsigned int);

      virtual double GetBeta_2_2(unsigned int, unsigned int);

      virtual double SmodFull(unsigned int);
      virtual double Smod(unsigned int);
      virtual double Smod(unsigned int, unsigned int, unsigned int);
      virtual void Print_Info();
  };
}

#endif
