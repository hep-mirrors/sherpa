//bof
//Version: 2 ADICIC++-0.0/2004/08/06    (surely preliminarily)

//Classes are provided, which manage the calculation of the kinematics
//according to given recoil strategies.



#ifndef _Recoil_Calculator_H_
#define _Recoil_Calculator_H_ _Recoil_Calculator_H_


#include <enumextra>
#include "Vector.H"


#define RECOIL_CALCULATOR_OUTPUT RECOIL_CALCULATOR_OUTPUT
#undef  RECOIL_CALCULATOR_OUTPUT





namespace ADICIC {



  struct Recoil_Setup;





  class Recoil_Calculator {

  private:

    Recoil_Calculator(const Recoil_Calculator&);
    Recoil_Calculator& operator=(const Recoil_Calculator&);
    const bool operator==(const Recoil_Calculator&) const;

  private:

    static int s_count;
    static const ATOOLS::Vec4D s_zaxis;

  protected:

    Trio f_recoil;

    ATOOLS::Vec4D m_p1;
    ATOOLS::Vec4D m_p3;

 public:

    inline Recoil_Calculator();
    virtual ~Recoil_Calculator();

    static const int&           InStore;
    static const ATOOLS::Vec4D& ZAxis;

    inline void GetResult(Trio&, ATOOLS::Vec4D&, ATOOLS::Vec4D&) const;

    virtual void Which() const=0;
    virtual const bool GenerateCmsMomenta(const Recoil_Setup&,
					  const ATOOLS::Vec4D&)=0;

  };





  template<class _Strategy> class Recoil : public Recoil_Calculator {

  private:

    //Blocked methods:

    Recoil(const Recoil<_Strategy>&);
    Recoil<_Strategy>& operator=(const Recoil<_Strategy>&);
    const bool operator==(const Recoil<_Strategy>&) const;

  private:

    double m_costheta;
    double m_sintheta;
    double m_phi;

    const Recoil_Setup* p_ini;

    ATOOLS::Vec4D m_cmsaxis;

    //Private methods:

    const bool Initialize();
    void RotateOnto(const ATOOLS::Vec4D&);
    void Rotate(const ATOOLS::Vec4D&, const ATOOLS::Vec4D&);

    const bool Calculate();

  public:

    Recoil();
    ~Recoil() {}

    //Public (inline) methods:

    inline void Which() const;
    inline const bool GenerateCmsMomenta(const Recoil_Setup&,
					 const ATOOLS::Vec4D&);

  };





  class Dipole_Handler;





  struct Recoil_Setup {
  private:
    friend class Dipole_Handler;
    //friend class Recoil_Calculator;
    //Derived classes are not instantly set to be friends as well.
  private:
    double E1, E2, E3;
    double Q1, Q2, Q3;
    //The Q numbers will be needed as recently as the massive case
    //has been implemented.
  public:
    inline const double GetE1() const;
    inline const double GetE2() const;
    inline const double GetE3() const;
    inline const double GetQ1() const;
    inline const double GetQ2() const;
    inline const double GetQ3() const;
  };



}    //eo namespace ADICIC





#include "Recoil_Calculator.inl.hh"


#endif    //eo _Recoil_Calculator_H_



//eof
