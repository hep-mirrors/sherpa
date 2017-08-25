#ifndef AddOns_OpenLoops_OpenLoops_Dipoles_H
#define AddOns_OpenLoops_OpenLoops_Dipoles_H

#include "PHASIC++/Process/CS_Dipole.H"

namespace PHASIC {
  class Dipole_Info;
  class FF_Dipole_Kinematics;
}

namespace OpenLoops {
  
  class Dipole : virtual public PHASIC::CS_Dipole {

  public:

    Dipole(const PHASIC::Dipole_Info& di, const int& id,
	   const AmplitudeType& type);
    
    ~Dipole() {}
    
  protected:

    int m_ol_id;
    double m_prefac;
    AmplitudeType m_amptype;
    const AmplitudeType AmpType() const { return m_amptype; }

    /* Calculate fill correlator using methods below */
    double CalcCorrelator() const;

    /* Calc \tilde{p}^\mu similar to table 1 of arXiv:0709.2881v1 */
    virtual ATOOLS::Vec4D CalcPtilde() const = 0;

    /* Prefactor of spin- and color-correlated contribution, enters as
       CalcB() * <1,...,m;a,b| ptilde^\mu T_ij T_k ptilde^\nu |b,a;m,...m1> * T_ij^{-2} * ptilde^{-2} */
    virtual double CalcB() const = 0;

    /* Prefactor of purely color-correlated contribution,  enters as
       CalcA() * <1,...,m;a,b| T_ij T_k |b,a;m,...m1> * T_ij^{-2} */
    virtual double CalcA() const = 0;

  };

  
  ////////////////////////////// // ///////////////////////////////////////
  ////////////////////////////// FF ///////////////////////////////////////
  ////////////////////////////// // ///////////////////////////////////////
  

  class FF_Dipole: public OpenLoops::Dipole, public PHASIC::FF_Dipole {

  public:

    FF_Dipole(const PHASIC::Dipole_Info& pi,
	      const int& id, const AmplitudeType& type) :
      PHASIC::CS_Dipole(pi), PHASIC::FF_Dipole(pi), Dipole(pi,id,type) { }
    
    ~FF_Dipole() {}

    double Calc() const;

  private:

    double CalcB() const;
    double CalcA() const;
    ATOOLS::Vec4D CalcPtilde() const;

  };


  ////////////////////////////// // ///////////////////////////////////////
  ////////////////////////////// FI ///////////////////////////////////////
  ////////////////////////////// // ///////////////////////////////////////
  

  class FI_Dipole: public OpenLoops::Dipole, public PHASIC::FI_Dipole {

  public:

    FI_Dipole(const PHASIC::Dipole_Info& pi,
	      const int& id, const AmplitudeType& type) :
      PHASIC::CS_Dipole(pi), PHASIC::FI_Dipole(pi), Dipole(pi,id,type) { }

    ~FI_Dipole() {}

    double Calc() const;

  private:

    double CalcB() const;
    double CalcA() const;
    ATOOLS::Vec4D CalcPtilde() const;

  };


  ////////////////////////////// // ///////////////////////////////////////
  ////////////////////////////// IF ///////////////////////////////////////
  ////////////////////////////// // ///////////////////////////////////////
  

  class IF_Dipole: public OpenLoops::Dipole, public PHASIC::IF_Dipole {

  public:

    IF_Dipole(const PHASIC::Dipole_Info& pi,
	      const int& id, const AmplitudeType& type) :
      PHASIC::CS_Dipole(pi), PHASIC::IF_Dipole(pi), Dipole(pi,id,type) { }

    ~IF_Dipole() {}

    double Calc() const;

  private:

    double CalcB() const;
    double CalcA() const;
    ATOOLS::Vec4D CalcPtilde() const;

  };


  ////////////////////////////// // ///////////////////////////////////////
  ////////////////////////////// II ///////////////////////////////////////
  ////////////////////////////// // ///////////////////////////////////////
  

  class II_Dipole: public OpenLoops::Dipole, public PHASIC::II_Dipole {

  public:

    II_Dipole(const PHASIC::Dipole_Info& pi,
	      const int& id, const AmplitudeType& type) :
      PHASIC::CS_Dipole(pi), PHASIC::II_Dipole(pi), Dipole(pi,id,type) { }

    ~II_Dipole() {}

    double Calc() const;

  private:

    double CalcB() const;
    double CalcA() const;
    ATOOLS::Vec4D CalcPtilde() const;
    
  };

}

#endif
