#ifndef QCD_Currents_H
#define QCD_Currents_H

#include "Current.H"
#include "C_Vector.H"
#include "C_Tensor.H"
#include "C_Spinor.H"
#include "Spinor.H"

namespace EXTRAXS {

  typedef std::vector<ATOOLS::CVec4D>  CVec_Vector;
  typedef std::vector<ATOOLS::CAsT4D>  CAsT_Vector;
  typedef std::vector<ATOOLS::CSpinor> CSpinor_Vector;

  class QCD_G: public Current<ATOOLS::CVec4D> {
  protected:

    std::string CLabel() const;

  private:

    ATOOLS::Vec4D  m_k;
    ATOOLS::Spinor m_kp, m_km;

    ATOOLS::CVec4D EM(const ATOOLS::Vec4D &p,const int cr,const int ca);
    ATOOLS::CVec4D EP(const ATOOLS::Vec4D &p,const int cr,const int ca);

  public:

    inline QCD_G(): Current<ATOOLS::CVec4D>(ATOOLS::kf::gluon) {}

    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca);
    void SetGauge(const ATOOLS::Vec4D &k);

    void AddPropagator();

    Complex Contract(const Current_Base &c,
		     const long unsigned int &hm,
		     const long unsigned int &hp) const;

    char Type() const;    

  };// end of class QCD_G

  class QCD_QGC: public Current<ATOOLS::CAsT4D> {
  protected:

    std::string CLabel() const;

  public:

    inline QCD_QGC(): Current<ATOOLS::CAsT4D>(ATOOLS::kf::gluonqgc) {}

    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca);
    void SetGauge(const ATOOLS::Vec4D &k);

    void AddPropagator();

    Complex Contract(const Current_Base &c,
		     const long unsigned int &hm,
		     const long unsigned int &hp) const;

    char Type() const;    

  };// end of class QCD_QGC

  class QCD_Q: public Current<ATOOLS::CSpinor> {
  protected:

    std::string CLabel() const;

  public:

    inline QCD_Q(const ATOOLS::Flavour &fl): 
      Current<ATOOLS::CSpinor>(fl) {}

    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca);
    void SetGauge(const ATOOLS::Vec4D &k);

    void AddPropagator();

    Complex Contract(const Current_Base &c,
		     const long unsigned int &hm,
		     const long unsigned int &hp) const;

    char Type() const;    

  };// end of class QCD_Q

}// end of namespace EXTRAXS

#endif
