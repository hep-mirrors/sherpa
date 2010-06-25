#ifndef PHOTONS_MEs_PHOTONS_ME_Base_H
#define PHOTONS_MEs_PHOTONS_ME_Base_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "PHOTONS++/Main/Dipole_Type.H"

namespace ATOOLS {
  class Flavour;
  class Poincare;
}

namespace PHOTONS {

  class PHOTONS_ME_Base {
    protected:
      std::string        m_name;

      unsigned int       m_nV;
      unsigned int       m_nR;
      double             m_alpha;
      double             m_e;
      double             m_sW;
      double             m_cW;
      double             m_GF;
      double             m_sqrt2;
      Complex    m_i;

      Complex    m_cL;
      Complex    m_cR;

      ATOOLS::Poincare * p_boost;
      ATOOLS::Poincare * p_rot;

      ATOOLS::Flavour    m_flavs[9];
      double             m_masses[9];
      ATOOLS::Vec4D *    m_moms;
      ATOOLS::Vec4D      m_moms0[9];
      ATOOLS::Vec4D      m_moms1[9][9];
      ATOOLS::Vec4D      m_moms2[9][9];
      int                m_spins[9];

      const Particle_Vector_Vector& m_pvv_zero;

      virtual void    BoostOriginalPVVToMultipoleCMS() = 0;

      virtual Complex InfraredSubtractedME_0_0() = 0;
      virtual Complex InfraredSubtractedME_0_1() = 0;
      virtual Complex InfraredSubtractedME_0_2() = 0;

      virtual Complex InfraredSubtractedME_1_05(unsigned int) = 0;
      virtual Complex InfraredSubtractedME_1_15(unsigned int) = 0;

      virtual Complex InfraredSubtractedME_2_1(unsigned int,
                                               unsigned int) = 0;

    public:
      PHOTONS_ME_Base(const Particle_Vector_Vector&);
      virtual ~PHOTONS_ME_Base();


      virtual double GetBeta_0_0() = 0;
      virtual double GetBeta_0_1() = 0;
      virtual double GetBeta_0_2() = 0;

      virtual double GetBeta_1_1(unsigned int) = 0;
      virtual double GetBeta_1_2(unsigned int) = 0;

      virtual double GetBeta_2_2(unsigned int, unsigned int) = 0;

      virtual void   FillMomentumArrays(const Particle_Vector_Vector&) = 0;
      virtual double Smod(unsigned int) = 0;

      inline std::string Name() { return m_name; }

      static PHOTONS_ME_Base * GetIRsubtractedME
                               (const Particle_Vector_Vector& pvv);
      static PHOTONS_ME_Base * GetIRsubtractedME
                               (const std::string& tag,
                                const Particle_Vector_Vector& pvv);
  };
}

#define DECLARE_PHOTONS_ME_GETTER(NAME,TAG)                        \
  DECLARE_GETTER(NAME,TAG,PHOTONS::PHOTONS_ME_Base,PHOTONS::Particle_Vector_Vector); \
  void NAME::PrintInfo(std::ostream &str,const size_t width) const \
  {                                                                \
    str<<#TAG;                                                     \
  }

typedef ATOOLS::Getter_Function<PHOTONS::PHOTONS_ME_Base,PHOTONS::Particle_Vector_Vector>
        PHOTONS_ME_Getter;

#endif
