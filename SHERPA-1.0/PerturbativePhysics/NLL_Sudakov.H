#ifndef NLL_Sudakov_H
#define NLL_Sudakov_H

//#include "Sudakov_Tools.H"
#include "Running_AlphaS.H"
#include "Fast_Function.H"
#include "Gauss_Integrator.H"
#include "Kabbala.H"

namespace SHERPA {
  class NLL {
  public:
    enum code { G_g_G_f    = 0,
		G_q        = 1,
		G_g        = 2,
                G_f        = 3,
		//	what about Delta_q, Delta_g, Delta_f
		// direct versions - LLA  2/(2 pi q) * Int( p_?_lla)
		G_g_lla    = 10,     
		G_q_lla    = 11,
		// direct versions - MLLA  2/(2 pi q) * Int( p_?_mlla)
		G_g_mlla   = 20,
		G_q_mlla   = 21,

                three_jet   = 30,   //  G_q(Qmax,Q)*Delta_g(Q) 

		four_jet_1 = 41,  //  G_q(Qmax,Q)*Delta_g(Q) * Int( NLL::four_jet_2)
		four_jet_2 = 42,  //  G_g*Delta_g   and/or   G_f*Delta_f

		five_jet_1 = 51,  //  G_q * Delta_g * ((Int( NLL::five_jet_2 ))^2 
		                  // + Int( NLL::five_jet_3) [ + Int( NLL::five_jet_4) ])
		five_jet_2 = 52,  //  ( G_g * Delta_g )
                five_jet_3 = 53,  //  G_g * Delta_g * (Int( NLL::five_jet_5))
		five_jet_4 = 54,  //  2*G_f*Dleta_f*(Int( NLL::five_jet_6))
                five_jet_5 = 55,  //  2*G_g*Delta_g  and/or  4*G_f*Delta_f
                five_jet_6 = 56,  //  2*(G_q - G_g + Gq)*Delta_g

                P_g_lla    = 200, //   P_g_LLA(qstat,tt)     tt=z
		P_q_lla    = 201, //   P_q_LLA(qstat,tt)     tt=z
		P_g_mlla   = 210, //   P_g_MLLA(qstat,tt);   tt=z
		P_q_mlla   = 211, //   P_q_MLLA(qstat,tt);   tt=z

		Unknown    = 99,

		int_G_g_G_f= 300,
		int_G_q    = 301,
		int_G_g    = 302,
		int_G_f    = 303,

		// select methode for Sudakov determination
		use_table   = 900,
		inte_self   = 901,
                analytic    = 902,
                alphas_repl = 903
    };
  };

  class NLL_Sudakov : public AMATOOLS::Function_Base {
    double tmax,tmin;
    NLL::code            inte;    // specifies integration algorithm used for Deltas
    AMATOOLS::Gauss_Integrator gauss;    // integrator routine
    double qmax,qmin,qstat;
    int    Nc,Nf;
    double CA,CF,TR;
    double K;
    double beta0,beta1,lambda2,mu2,asmu;
    bool   nlo,run_as;
    
    void   FixLambda2();

    NLL::code            type;    // specifies return value of operator()

    // logs of sudakovs
    AMATOOLS::Fast_Function log_delta_q, log_delta_g;

    // branching probabilities (branch at q, with start scale Q)
    double GammaQ(double q, double Q);
    double GammaG(double q, double Q);
    double GammaF(double q);

    // integrated branching probabilities 
    double IntGammaI(double q, double i0, double i1, double i2, double i3);

    double IntGammaQ(double Q0, double Q);
    double IntGammaG(double Q0, double Q);
    double IntGammaF(double Q0, double Q);

    // splitting functions
    double P_g_LLA(double q, double z);
    double P_q_LLA(double q, double z);
    double P_g_MLLA(double q, double z);
    double P_q_MLLA(double q, double z);

  public:
    NLL_Sudakov(double,double);
    // Sudakovs (default: relying on fast interpolatoin of log sudakov)
    // no branching between Q (high scale) and q (low scale)
    double DeltaQ(double Q, double q);
    double DeltaG(double Q, double q);
    double DeltaF(double Q, double q);

    double R2(double Q0, double Q); 
    double R3(double Q0, double Q); 
    double R4(double Q0, double Q); 
    double R5(double Q0, double Q); 

    // Kabbala-Sudakovs, do not only return values but also a string
    AMEGIC::Kabbala DeltaQ(int,double Q, double q);
    AMEGIC::Kabbala DeltaG(int,double Q, double q);
    AMEGIC::Kabbala DeltaF(int,double Q, double q);

    // LL AlphaS routine
    double AlphaS(double);

    // general function, depending on specified type.
    double operator()(double);

    // consitency checks
    void CheckSudakovs();
  };


}
#endif
