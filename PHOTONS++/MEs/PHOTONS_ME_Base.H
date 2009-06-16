#ifndef PHOTONS_MEs_PHOTONS_ME_Base_H
#define PHOTONS_MEs_PHOTONS_ME_Base_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/MyComplex.H"
#include "HELICITIES/Main/C_Vector.H"
#include "HELICITIES/Main/XYZFuncs.H"
#include "HELICITIES/Main/Polarization_Tools.H"
#include "PHOTONS++/Main/Dipole_Type.H"
#include "ATOOLS/Org/Run_Parameter.H"


namespace PHOTONS {

  class PHOTONS_ME_Base {
    protected:
      unsigned int  m_nV;
      unsigned int  m_nR;
      double        m_e;
      double        m_alpha;
      double        m_sW;
      double        m_cW;
      Complex       m_i;

      Complex  m_cL;
      Complex  m_cR;

      Poincare * p_boost;
      Poincare * p_rot;

      Flavour  m_flavs[9];
      Vec4D *  m_moms;
      Vec4D    m_moms0[9];
      Vec4D    m_moms1[9][9];
      Vec4D    m_moms2[9][9];
      int      m_spins[9];

      Particle_Vector_Vector m_pvv_zero;
      Particle_Vector_Vector m_pvv_one;

      virtual void    BoostOriginalPVVToMultipoleCMS() = 0;
      virtual void    FillMomentumArrays() = 0;


      virtual Complex InfraredSubtractedME_0_0() = 0;
      virtual Complex InfraredSubtractedME_0_1() = 0;
      virtual Complex InfraredSubtractedME_0_2() = 0;

      virtual Complex InfraredSubtractedME_1_05(unsigned int) = 0;
      virtual Complex InfraredSubtractedME_1_15(unsigned int) = 0;

      virtual Complex InfraredSubtractedME_2_1(unsigned int, unsigned int) = 0;

    public:
      PHOTONS_ME_Base();
      virtual ~PHOTONS_ME_Base();


      virtual double GetBeta_0_0() = 0;
      virtual double GetBeta_0_1() = 0;
      virtual double GetBeta_0_2() = 0;

      virtual double GetBeta_1_1(unsigned int) = 0;
      virtual double GetBeta_1_2(unsigned int) = 0;

      virtual double GetBeta_2_2(unsigned int, unsigned int) = 0;

      virtual double Smod(unsigned int) = 0;
  };

}
#endif
