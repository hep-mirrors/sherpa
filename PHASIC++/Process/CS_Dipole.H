#ifndef PHASIC_Process_CS_Dipole_H
#define PHASIC_Process_CS_Dipole_H

#include <string>
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

namespace MODEL {
  class Coupling_Map;
  class Coupling_Data;
}

namespace EXTAMP {
  class Dipole_Wrapper_Process;
}

namespace PHASIC {


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
    virtual double Calc() const = 0;
    virtual void CalcKinematics(const ATOOLS::Vec4D_Vector& p) = 0;
    virtual const ATOOLS::Vec4D_Vector& Momenta() const = 0;

    bool PassesAlphaCuts() const;
    bool PassesAlphaMin () const;
    
    /* Getter functions instantiating a concrete implementation and
       returning a pointer to it */
    static CS_Dipole* Get(const Dipole_Info& di);
    static CS_Dipole* Get(const std::string& tag,
			  const Dipole_Info& di);

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

    /* Give the dipole an pointers to Coupling_Data instances, where
       they can collect their (pre-calculated) running couplings */
    void SetCouplings(MODEL::Coupling_Data* p_rqcd,
		      MODEL::Coupling_Data* p_rqed);

    /* Need access methods in order to synchronize these couplings
       with wrapper processes */
    MODEL::Coupling_Data* RunningQCD() { return p_aqcd; }
    MODEL::Coupling_Data* RunningQED() { return p_aqed; }

    const int& SubtractionType() const   { return m_dip_info.m_subtype; }
    void SetSubtractionType(int subtype) { m_dip_info.m_subtype = subtype; }

    virtual const Dipole_Kinematics* const LastKinematics() const = 0;

  protected:

    /* Determine the mother flavour ij of QCD splitting ij -> i+j */
    static ATOOLS::Flavour CombinedFlavour(const size_t& i, const size_t& j,
					   const ATOOLS::Flavour_Vector& flavs);

    ATOOLS::Flavour_Vector m_born_flavs;

    /* An ID vector encoding the id's of particles as they are ordered
       in the born flavour vector m_born_flavs. This is needed for
       ATOOLS::NLO_Subevents, which own a pointer to such a vector. */
    std::vector<size_t> m_id_vector;

    double m_symfac;

    Dipole_Info m_dip_info;

    MODEL::Coupling_Data* p_aqcd, * p_aqed;

    constexpr static double m_CF = 4./3.;
    constexpr static double m_CA = 3.0;
    constexpr static double m_TR = 1./2.;

  };


  ///////////////////////////////////////////////////////////////
  ////////// FINAL FINAL ////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  

  /* Kinematic variables of hep-ph/9605323v3 eq. (5.3) - (5.6) */
  struct FF_Dipole_Kinematics : public Dipole_Kinematics {
    ATOOLS::Vec4D m_pi;
    ATOOLS::Vec4D m_pj;
    ATOOLS::Vec4D m_pk;
    ATOOLS::Vec4D m_pij_tilde;
    ATOOLS::Vec4D m_pk_tilde;
    double m_zi;
    double m_zj;
    double m_y;

    double Alpha() const { return m_y; }

    double ShowerX()  const { return m_zi; }
    double ShowerY()  const { return m_y ; }
    double ShowerQ2() const { return (m_pi+m_pj+m_pk).Abs2(); }
  };


  class FF_Dipole : virtual public CS_Dipole {

  public:

    FF_Dipole(const Dipole_Info& di) : CS_Dipole(di) {};

    virtual ~FF_Dipole() {};

    const ATOOLS::Vec4D_Vector& Momenta() const { return m_kin.m_born_mom; }

  protected:

    void CalcKinematics(const ATOOLS::Vec4D_Vector& p);

    const Dipole_Kinematics* const LastKinematics() const { return &m_kin; }

    /* Legacy structure: NLO_subevts only maintain pointers to const
       momenta and flavours. These pointers need to point to
       persistent copies of momenta/flavours that do not run out of
       scope. Therefore have to maintain a copy of the last calculated
       momentum configuration in this class and let the momentum
       pointer of one NLO_subevt point to that instance. */
    FF_Dipole_Kinematics m_kin;

  };


  ///////////////////////////////////////////////////////////////
  ////////// FINAL INITIAL //////////////////////////////////////
  ///////////////////////////////////////////////////////////////


  /* Kinematic variables of hep-ph/9605323v3 eq. (5.3) - (5.6) */
  struct FI_Dipole_Kinematics : public Dipole_Kinematics {
    ATOOLS::Vec4D m_pi;
    ATOOLS::Vec4D m_pj;
    ATOOLS::Vec4D m_pa;
    ATOOLS::Vec4D m_pij_tilde;
    ATOOLS::Vec4D m_pa_tilde;
    double m_x;
    double m_zi;
    double m_zj;

    double Alpha() const { return 1.-m_x; }

    double ShowerX()  const { return m_zi; }
    double ShowerY()  const { return m_x ; }
    double ShowerQ2() const { return (m_pi+m_pj-m_pa).Abs2(); }
  };


  class FI_Dipole : virtual public CS_Dipole {

  public:

    FI_Dipole(const Dipole_Info& di) : CS_Dipole(di) {};

    virtual ~FI_Dipole() {};

    const ATOOLS::Vec4D_Vector& Momenta() const { return m_kin.m_born_mom; }

  protected:

    void CalcKinematics(const ATOOLS::Vec4D_Vector& p);

    const Dipole_Kinematics* const LastKinematics() const { return &m_kin; }

    FI_Dipole_Kinematics m_kin;

  };


  ///////////////////////////////////////////////////////////////
  ////////// INITIAL FINAL //////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  

  /* Kinematic variables of hep-ph/9605323v3 eq. (5.3) - (5.6) */
  struct IF_Dipole_Kinematics : public Dipole_Kinematics {
    ATOOLS::Vec4D m_pa;
    ATOOLS::Vec4D m_pi;
    ATOOLS::Vec4D m_pk;
    ATOOLS::Vec4D m_pai_tilde;
    ATOOLS::Vec4D m_pk_tilde;
    double m_ui;
    double m_x;

    double Alpha() const { return m_ui; }

    double ShowerX()  const { return m_x; }
    double ShowerY()  const { return m_ui; }
    double ShowerQ2() const { return (-m_pa+m_pi+m_pk).Abs2(); }
  };


  class IF_Dipole : virtual public CS_Dipole {

  public:

    IF_Dipole(const Dipole_Info& di) : CS_Dipole(di) {};

    virtual ~IF_Dipole() {};

    const ATOOLS::Vec4D_Vector& Momenta() const { return m_kin.m_born_mom; }

  protected:

    void CalcKinematics(const ATOOLS::Vec4D_Vector& p);

    const Dipole_Kinematics* const LastKinematics() const { return &m_kin; }

    IF_Dipole_Kinematics m_kin;

  };


  ///////////////////////////////////////////////////////////////
  ////////// INITIAL INITIAL ////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  

  /* Kinematic variables of hep-ph/9605323v3 eq. (5.3) - (5.6) */
  struct II_Dipole_Kinematics : public Dipole_Kinematics {
    ATOOLS::Vec4D m_pa;
    ATOOLS::Vec4D m_pi;
    ATOOLS::Vec4D m_pb;
    ATOOLS::Vec4D m_pai_tilde;
    ATOOLS::Vec4D m_pb_tilde;
    double m_ui;
    double m_x;
    double m_v;

    double Alpha() const { return m_v; }

    double ShowerX()  const { return m_x; }
    double ShowerY()  const { return m_v; }
    double ShowerQ2() const { return (-m_pa+m_pi-m_pb).Abs2(); }
  };


  class II_Dipole : virtual public CS_Dipole {

  public:

    II_Dipole(const Dipole_Info& di) : CS_Dipole(di) {};

    virtual ~II_Dipole() {};

    const ATOOLS::Vec4D_Vector& Momenta() const { return m_kin.m_born_mom; }

  protected:

    void CalcKinematics(const ATOOLS::Vec4D_Vector& p);

    const Dipole_Kinematics* const LastKinematics() const { return &m_kin; }
    
    II_Dipole_Kinematics m_kin;

  };

}

#define DECLARE_CSDIPOLE_GETTER(NAME,TAG)			    \
  DECLARE_GETTER(NAME,TAG,PHASIC::CS_Dipole,PHASIC::Dipole_Info);   \
  void ATOOLS::Getter<PHASIC::CS_Dipole,PHASIC::Dipole_Info,NAME>:: \
  PrintInfo(std::ostream &str,const size_t width) const		    \
  {								    \
    str<<#TAG;							    \
  }

#endif
