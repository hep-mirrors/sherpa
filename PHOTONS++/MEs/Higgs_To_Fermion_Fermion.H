#ifndef PHOTONS_MEs_Higgs_To_Fermion_Fermion_H
#define PHOTONS_MEs_Higgs_To_Fermion_Fermion_H

#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "PHOTONS++/Tools/Dipole_FF.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__YFS_NNLO
#include "PHOTONS++/MEs/RVTools/Higgs_Decay_Two_Loop_Functions.H"
#endif

namespace PHOTONS {
  class Higgs_To_Fermion_Fermion: public PHOTONS_ME_Base, public Dipole_FF {
    private:
      bool      m_switch;
      double m2, m2p, Qf, Qfp, If, Ifp, m_smod;
      Complex m_x, m_r_res, m_r2_res, m_r2s_res;
#ifdef USING__YFS_NNLO
      Higgs_Decay_Two_Loop_Functions * p_TL;
      Higgs_Decay_Two_Loop_Functions * p_TL_scaled;
#endif
      METOOLS::DivArrC m_res,m_nlo_res,m_rv_res,m_rv_res_scaled,m2_rv,m2_rvfull,m_nnlo_res,m2_loop;


      virtual void    BoostOriginalPVVToMultipoleCMS();


      virtual Complex InfraredSubtractedME_0_0();
      virtual Complex InfraredSubtractedME_0_1(const int& = 1);
      virtual Complex InfraredSubtractedME_0_2();

      virtual Complex InfraredSubtractedME_1_05(unsigned int, const double& =1.);
      virtual Complex InfraredSubtractedME_1_15(unsigned int, const double& =1.);

      virtual Complex InfraredSubtractedME_2_1(unsigned int, unsigned int);

    public:
      Higgs_To_Fermion_Fermion(const Particle_Vector_Vector&);
      virtual ~Higgs_To_Fermion_Fermion();

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
