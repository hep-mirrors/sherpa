#ifndef SHRIMPS_Eikonals_Form_Factors_H
#define SHRIMPS_Eikonals_Form_Factors_H

#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Function_Base.H"
#include <vector>

namespace SHRIMPS {
  /*!
    \class Form_Factor
    \brief The form factor and its Fourier transform, yielding the boundary 
    conditions of the differential equations defining the single channel 
    eikonals.  

    The two following forms are implemented (and tested):
    - Gaussian
      \f[
      F(q) = \beta_0^2\,(1+\kappa)\,\exp\left(-\frac{(1+\kappa)q^2}{\Lambda^2}\right)
      \f]
    - dipole
      \f[
      F(q) = \beta_0^2\,(1+\kappa)\,
             \frac{\exp\left[-\xi\frac{(1+\kappa)q^2}{\Lambda^2}\right]}
	          {\left[1+\frac{(1+\kappa)q^2}{\Lambda^2}\right]^2}\,.
      \f]
    
    The most important method here is Form_Factor::FourierTransform(const double & b) 
    yielding the Fourier transform of the form factor at a given impact parameter 
    \f$b\f$.  It has been tested to be able to achieve practically arbitrary
    precision.  
  */
  class Form_Factor : public ATOOLS::Function_Base {
  private:
    class Norm_Argument : public ATOOLS::Function_Base {
    private:
      Form_Factor * p_ff;
    public:
      Norm_Argument(Form_Factor * ff) : p_ff(ff) {}
      double operator()(double q);
    };

    class FT_Argument : public ATOOLS::Function_Base {
    private:
      Form_Factor * p_ff;
      double        m_b;
    public:
      FT_Argument() {}
      FT_Argument(Form_Factor * ff) : p_ff(ff), m_b(0.) {}
      void   SetB(const double & b) { m_b = b; }
      double operator()(double q);
    };
  private:
    FT_Argument m_ftarg;

    ff_form::code m_form;
    int           m_number;
    double        m_Lambda2, m_beta, m_kappa, m_xi;
    double        m_prefactor, m_ftnorm, m_norm;
    double        m_bmax, m_deltab;
    size_t        m_bsteps;
    double        m_ffmin, m_ffmax;
    double        m_accu;
    int           m_test;


    std::vector<double> m_values;

    double Norm();
    double NormAnalytical();
    double CalculateFourierTransform(const double & b);
    double AnalyticalFourierTransform(const double & b);
    void   FillFTGrid();
    void   PrintFFGrids(const int & mode);

    void   TestBessel() const;
    void   TestFormFactor();
  public:
    Form_Factor(const int & number,const int & test=0);
    Form_Factor(const ff_form::code & fff,const int & number,
		const int & test=0);
    inline ~Form_Factor() {}
    
    void   Initialise(const std::vector<double> & params);
    double operator()(const double q2);
    double operator()() { return 1.; }
    double FourierTransform(const double & b) const;
    double ImpactParameter(const double & val) const;
    double SelectQT2(const double & qt2max,const double & qt2min=-1.) const;

    const double & Bmax() const      { return m_bmax;      }
    const double & DeltaB() const    { return m_deltab;    }
    const size_t & Bbins() const     { return m_bsteps;    }
    const double & FFmin() const     { return m_ffmin;     }
    const double & FFmax() const     { return m_ffmax;     }
    const double & Prefactor() const { return m_prefactor; }
    const double & Beta0() const     { return m_beta;      }
    const double & Kappa() const     { return m_kappa;     }
    const double & Lambda2() const   { return m_Lambda2;   }
    const double & Xi() const        { return m_xi;   }
    const int    & Number() const    { return m_number;    }
  };
}

#endif
