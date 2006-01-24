//bof
//Version: 3 ADICIC++-0.0/2005/07/25    (surely preliminarily)

//Classes are provided, which manage the calculation of the kinematics
//according to given recoil strategies.



#ifndef _Recoil_Calculator_H_
#define _Recoil_Calculator_H_ _Recoil_Calculator_H_


#include <vector>
#include <enumextra>
#include "Vector.H"
#include "Recoil_Strategy.hpp"
#include "Dipole_Handler.H"


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

    const Dipole*         p_dip;
    const Sudakov_Result* p_sur;

    Recoil_Result* p_rer;

 public:

    inline Recoil_Calculator();
    virtual ~Recoil_Calculator();

    static const int&           InStore;
    static const ATOOLS::Vec4D& ZAxis;

    virtual const Recoil_Strategy::Type IsType() const=0;
    virtual void Which() const=0;
    virtual void TestKey(Dipole_Handler::Key) const=0;
    virtual const bool GenerateMomenta(const Dipole&, const Sudakov_Result&,
				       Recoil_Result&)=0;

  };





  template<Recoil_Strategy::Type _Strategy>
  class Recoil : public Recoil_Calculator {

  private:

    //Blocked methods:
    Recoil(const Recoil<_Strategy>&);
    Recoil<_Strategy>& operator=(const Recoil<_Strategy>&);
    const bool operator==(const Recoil<_Strategy>&) const;

  private:

    //The following works, but is not needed.
    //static const Recoil_Strategy::Type s_type=_Strategy;

    double m_costheta;
    double m_sintheta;
    double m_phi;

    Multidouble m_e;

    //Private methods:
    inline const bool Generate();    //Specific!

    void CrossProductTest(const ATOOLS::Vec4D&) const;
    void RotateOnto(const ATOOLS::Vec4D&);
    void Rotate(const ATOOLS::Vec4D&, const ATOOLS::Vec4D&);

    const bool CmsMomenta();
    const bool LabMomenta();

    const bool Calculate();    //Specific!

  public:

    Recoil();
    ~Recoil();

    //Public (inline) methods:
    inline const Recoil_Strategy::Type IsType() const;
    inline void Which() const;    //Specific!
    inline void TestKey(Dipole_Handler::Key) const;    //Specific!
    inline const bool GenerateMomenta(const Dipole&, const Sudakov_Result&,
				      Recoil_Result&);

  };





  //Function
  void MakeRecos(const std::vector<bool>&, std::vector<Recoil_Calculator*>&);



}    //eo namespace ADICIC





#include "Recoil_Calculator.inl.hh"


#endif    //eo _Recoil_Calculator_H_



//eof
