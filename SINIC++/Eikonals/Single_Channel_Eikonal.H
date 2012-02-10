#ifndef SINIC_Eikonals_Single_Channel_Eikonal_H
#define SINIC_Eikonals_Single_Channel_Eikonal_H

#include "SINIC++/Eikonals/Form_Factors.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

namespace SINIC {
  struct deqmode {
    enum code {
      RungeKutta4Transformed = -4,
      RungeKutta4            = 4,
      RungeKutta2            = 2
    };
  };

  class Single_Channel_Eikonal : public ATOOLS::Function_Base {
    class Convolution2D : public ATOOLS::Function_Base {
      class Convolution1D : public ATOOLS::Function_Base {
      private:
	Single_Channel_Eikonal * p_eikonal;
	double                   m_b,m_b1,m_Y,m_y;
	double                   m_errmax12, m_errmax21;
	int                      m_test;
      public:
	Convolution1D(Single_Channel_Eikonal * eikonal, const double & Y,const int & test) : 
	  p_eikonal(eikonal), m_Y(Y), 
	  m_errmax12(0.), m_errmax21(0.), m_test(test) {}
	void   SetY(const double & y)   { m_y  = y;  }
	void   SetB(const double & b)   { m_b  = b;  }
	void   Setb1(const double & b1) { m_b1 = b1; }
	double operator()(double theta1);
	void   PrintErrors() { 
	  msg_Out()<<"Maximal errors in evaluating product of single terms: "<<std::endl<<"    "
		   <<"delta_max{Omega_12} = "<<m_errmax12<<", "
		   <<"delta_max{Omega_21} = "<<m_errmax12<<"."<<std::endl;
	}
      };
    private:
      Single_Channel_Eikonal   * p_eikonal;      
      double                     m_b,m_y,m_Y,m_accu;
      Convolution1D            * p_convolution1D;
      ATOOLS::Gauss_Integrator * p_integrator;
      int                        m_test;
    public:
      Convolution2D(Single_Channel_Eikonal * eikonal,const double & Y,const int & test) : 
	p_eikonal(eikonal), m_Y(Y), m_accu(1.e-3), m_test(test)
      {
	p_convolution1D = new Convolution1D(p_eikonal,m_Y,m_test);
	p_integrator    = new ATOOLS::Gauss_Integrator(p_convolution1D);
      }
      ~Convolution2D() {
	if (p_convolution1D) { delete p_convolution1D; p_convolution1D = NULL; }
	if (p_integrator)    { delete p_integrator;    p_integrator    = NULL; }
      }
      void   SetY(const double & y) { 
	m_y = y; 
	p_convolution1D->SetY(m_y); 
      }
      void   SetB(const double & b) { 
	m_b = b; 
	p_convolution1D->SetB(m_b); 
      }
      double operator()(double b1);
      void   PrintErrors() { p_convolution1D->PrintErrors(); }
    };
  private:
    int           m_test;
    deqmode::code m_deqmode;
    Form_Factor * p_ff1, * p_ff2;

    double m_Y, m_ycutoff, m_yshift;
    int    m_ybins;
    double m_deltay;
    double m_beta2, m_lambda, m_alpha, m_expfactor;
    double m_b1max, m_b2max, m_ff1max, m_ff2max, m_Bmax;
    int    m_ff1bins, m_ff2bins, m_Bbins;
    double m_deltaff1, m_deltaff2, m_deltaB;
    double m_accu, m_maxconv;
    bool   m_inelastic;

    std::vector<std::vector<std::vector<double> > >   m_grid1, m_grid2;
    std::vector<double>                               m_gridB;


    Convolution2D            * p_convolution2D;
    ATOOLS::Gauss_Integrator * p_integrator;


    void   ProduceInitialGrids();  
    void   InitialiseBoundaries(const int & i,const int & j,double & val1,double & val2);
    void   SolveSystem(const int & i,const int & j,double & val1,double & val2,const int & steps);
    void   RungeKutta2(const int & i, const int & j,double & val1,double & val2,const int & steps);
    void   RungeKutta4(const int & i, const int & j,double & val1,double & val2,const int & steps);
    void   RungeKutta4Transformed(const int & i, const int & j,
				  double & val1,double & val2,const int & steps);
    int    AdjustGrid(const int & i, const int & j,double & val1,double & val2);
    bool   CheckAccuracy(const int & i,const int & j,const int & ysteps,
			 const std::vector<double> & comp1,const std::vector<double> & comp2);

    void   PrintOmega_ik();
    void   TestSingleEikonal(const double & b1=0.,const double & b2=0.);
    void   TestDEQSolvers();
  public:
    Single_Channel_Eikonal(const deqmode::code & deq=deqmode::RungeKutta4,const int & m_test=0);
    ~Single_Channel_Eikonal();    

    void   Initialise(Form_Factor * ff1,Form_Factor * ff2,
		      const double & lambda,const double & alpha,
		      const double & Y,const double & ycutoff=0.);

    bool   GeneratePositions(const double & B,double & b1,double & theta1);

    double Omega12(const double & b1,const double & b2,const double & y,const bool & plot=false) const;
    double Omega21(const double & b1,const double & b2,const double & y,const bool & plot=false) const;
    double operator()(double B);

    void   ProduceImpactParameterGrid(const double & y);  
    double IntegrateOutImpactParameters(const double & B,const double & Y=0.);

    void SetInelastic(const bool & inelastic) { m_inelastic = inelastic; }
    void SetMaxConv(const double & maxconv)   { m_maxconv = maxconv; }

    const bool & Inelastic() const { return m_inelastic; }
    const double & MaxConv() const { return m_maxconv; }
    const double Prefactor() const { return ATOOLS::sqr(p_ff1->Prefactor()*p_ff2->Prefactor()); }

    const double & Beta2() const   { return m_beta2; }
    const double & Lambda2() const { return p_ff1->Lambda2(); }
    const double & Delta() const   { return m_alpha; }
    const double & Kappa_i() const { return p_ff1->Kappa(); }
    const double & Kappa_k() const { return p_ff2->Kappa(); }

    const Form_Factor * FF1() const { return p_ff1; }
    const Form_Factor * FF2() const { return p_ff1; }
  };
}

#endif
