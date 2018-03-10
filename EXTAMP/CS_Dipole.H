#ifndef EXTAMP_Process_CS_Dipole_H
#define EXTAMP_Process_CS_Dipole_H

#include <string>
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

namespace PHASIC {
  class Spin_Color_Correlated_ME2;
}

namespace EXTAMP {

  class Dipole_Wrapper_Process;
  class Spin_Color_Correlated_ME2;

  enum SplittingType { FF,IF,FI,II };
  std::ostream &operator<<(std::ostream &str,const SplittingType& st);


  enum FlavourType { gtogg,gtoqq,qtoqg };
  std::ostream &operator<<(std::ostream &str,const FlavourType& ft);
  

  struct Dipole_Info {
    
    Dipole_Info(const ATOOLS::Flavour_Vector& flavs,
		const size_t& i, const size_t& j, const size_t& k,
		const int& subtrtype, const double& alphamin, const double& alphamax);

    SplittingType          m_split_type;
    FlavourType            m_flav_type;
    ATOOLS::Flavour_Vector m_real_flavs;
    
    /* Indices i,j,k in momentum/flavour vector of corresponding real
       emission process */
    size_t m_real_i, m_real_j, m_real_k;
    
    /* Subtraction type: 0: plain Catani-Seymour
                         1: Dire
			 2: modified Catani-Seymour (Stefan) */
    int m_subtype;

    /* Cuts on alpha as defined in hep-ph/0307268 */
    double m_alphamin, m_alphamax;
    
  };


  std::ostream &operator<<(std::ostream &str,const Dipole_Info& di);


  struct Dipole_Kinematics {

    ATOOLS::Vec4D_Vector m_born_mom;

    /* Alpha parameter as defined in hep-ph/0307268 */
    virtual double Alpha() const = 0;

    /* Get shower-like variables */
    virtual double ShowerX()  const = 0;
    virtual double ShowerY()  const = 0;
    virtual double ShowerQ2() const = 0;

    bool PassesAlphaMin(const double& alphamin) const
    { return (Alpha() > alphamin); }

    /* Check whether alpha lies in the given interval */
    bool PassesAlphaCuts(const double& alphamin,
			 const double& alphamax) const
    {
      const double& alpha = Alpha();
      return (alpha > alphamin) && (alpha < alphamax);
    }
    
  };


  class CS_Dipole {

    friend class EXTAMP::Dipole_Wrapper_Process;

  public:

    CS_Dipole(const Dipole_Info& di);
    virtual ~CS_Dipole() {};

    const Dipole_Info& Info() const { return m_dip_info; }

    /* Calculate contribution to differential cross section: Call
       CalcKinematics first, then get XS from Calc() method.
       CalcKinematics stores results, resulting born momenta can be
       accessed through Momenta() method. */
    virtual void  CalcKinematics(const ATOOLS::Vec4D_Vector& p) = 0;
    virtual const ATOOLS::Vec4D_Vector& Momenta() const = 0;
    double  Calc() const;

    bool PassesAlphaCuts() const;
    bool PassesAlphaMin () const;
    
    /* Indices i,j,k in the real emission flavour config */
    const size_t& I() const { return m_dip_info.m_real_i; }
    const size_t& J() const { return m_dip_info.m_real_j; }
    const size_t& K() const { return m_dip_info.m_real_k; }

    /* Indices (ij) and k in the born flavour config */
    const size_t& BornIJ() const { return Emitter(); }
    size_t BornK () const { return (K()<Emitted()? K() : K()-1); }

    /* Convention: consider min(i,j) the emitter, max(i,j) the
       emitted */
    const size_t& Emitter() const { return std::min(I(),J()); }
    const size_t& Emitted() const { return std::max(I(),J()); }

    const std::vector<size_t>& IDVector() const {return m_id_vector; }

    const FlavourType&   FlavType()  const { return m_dip_info.m_flav_type;  }
    const SplittingType& SplitType() const { return m_dip_info.m_split_type; }

    const ATOOLS::Flavour_Vector& Flavours()     const {return m_born_flavs; }
    const ATOOLS::Flavour_Vector& RealFlavours() const {return m_dip_info.m_real_flavs; }
    const ATOOLS::Flavour& FlavI()  const { return m_dip_info.m_real_flavs[I()]; }
    const ATOOLS::Flavour& FlavJ()  const { return m_dip_info.m_real_flavs[J()]; }
    const ATOOLS::Flavour& FlavIJ() const { return m_born_flavs[BornIJ()]; }

    /* Given two indices i,j, and a flavour vector, construct born a
       flavour configuration by combining partons i and j */
    static ATOOLS::Flavour_Vector ConstructBornFlavours(const size_t& i, const size_t& j,
							const ATOOLS::Flavour_Vector& flavs);

    /* Construct an ID vector encoding the id's of particles as they
       are ordered in the born flavour vector m_born_flavs. This is
       needed for ATOOLS::NLO_Subevents, which own a pointer to such a
       vector. Note: id's are in binary encoding, meaning the i'th
       particle has id 1<<i and the combined particles i and j have
       inded (1<<i|1<<j). */
    static std::vector<size_t> ConstructIDVector(const size_t& i, const size_t& j,
						 const ATOOLS::Flavour_Vector& flavs);

    const int& SubtractionType() const   { return m_dip_info.m_subtype; }
    void SetSubtractionType(int subtype) { m_dip_info.m_subtype = subtype; }

    virtual const Dipole_Kinematics* const LastKinematics() const = 0;

    PHASIC::Spin_Color_Correlated_ME2* CorrelatedME() { return p_corr_me; }

  protected:

    /* Calc() relies on the implementation of all pure virtual
       functions below. They are implemented by child classes
       representing FF, II, FI, and IF dipoles. */
    double CalcCorrelator() const;
    virtual double CalcKinDependentPrefac() const = 0;
    /* Calc \tilde{p}^\mu similar to table 1 of arXiv:0709.2881v1 */
    virtual ATOOLS::Vec4D CalcPtilde() const = 0;
    /* Prefactor of spin- and color-correlated contribution, enters as
       CalcB() * <1,...,m;a,b| ptilde^\mu T_ij T_k ptilde^\nu |b,a;m,...m1> * T_ij^{-2} * ptilde^{-2} */
    virtual double CalcB() const = 0;
    /* Prefactor of purely color-correlated contribution,  enters as
       CalcA() * <1,...,m;a,b| T_ij T_k |b,a;m,...m1> * T_ij^{-2} */
    virtual double CalcA() const = 0;

    /* Determine the mother flavour ij of QCD splitting ij -> i+j */
    static ATOOLS::Flavour CombinedFlavour(const size_t& i, const size_t& j,
					   const ATOOLS::Flavour_Vector& flavs);

    constexpr static double m_CF = 4./3.;
    constexpr static double m_CA = 3.0;
    constexpr static double m_TR = 1./2.;

  private:
    
    PHASIC::Spin_Color_Correlated_ME2* p_corr_me;

    ATOOLS::Flavour_Vector m_born_flavs;

    /* An ID vector encoding the id's of particles as they are ordered
       in the born flavour vector m_born_flavs. This is needed for
       ATOOLS::NLO_Subevents, which own a pointer to such a vector. */
    std::vector<size_t> m_id_vector;

    /* Usual symmetry factor and kinematics-independent factor of
       splitting operators, e.g. prefactor of eq (5.8) of
       hep-ph/9605323v3 for case of FF dipole */
    double m_symfac;
    double m_const_prefac;

    Dipole_Info m_dip_info;

  };

}

#endif
