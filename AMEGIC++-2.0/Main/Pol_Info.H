#ifndef Pol_Info_h
#define Pol_Info_h

#include "Flavour.H"
#include "MathTools.H"
#define Explicit_Pols

//syntax for polarised particles: (in processes.dat)
//
//    degree of polarisation (0.=unpolarized,1.=fully polarized) (only 1 or 0 for outgoing)
//                        | l for linear
//                        | | polarization angle
//                        | | |
//   Bosons:   linear:   (.8l45)
//                        |
//                        | +/-/0(for massive bosons)
//                        |  |
//             circular: (.8c+)
//                        |
//                        | +/-
//                        |  |
//   Fermion helicities: (.8h+) (up to now only ok for massless)

namespace AMEGIC {

class mt {
public:
  enum momtype{  mom  = 10,  //external momentum
		 prop = 11,  //propagator
		 p_p  = 1,   //polarisation +1 (circular)
		 p_m  =-1,   //             -1 (circular)
		 p_l  = 0,   //             longitudinal
		 p_s  = 2,   //             scalar
		 p_si = 7,   //             massless scalar
		 p_l0 = 3,   //             linear with direction
		 p_l1 = 4,   //             linear, perpendicular to p_l0
	  	 p_lh = 5,   //             horizontaly (linear)
		 p_lp = 6,    //             perpendicular (linear)
		 p_none = 8
  };
};

  class Pol_Info{
  public:
    Pol_Info& operator=(const Pol_Info& p){
      if (this!=&p) {
	num    = p.num;
	p_type = p.p_type;
	angle  = p.angle;
	if(type)   delete[] type;
	if(factor) delete[] factor;
	type   = new int[num];
	factor = new double[num];
	for(int i=0;i<num;i++){
	  type[i]   = p.type[i];
	  factor[i] = p.factor[i];
	}
      }
      return *this;
    }
    Pol_Info() { num=0; type=0; factor=0; p_type=' '; angle=0.;}
    Pol_Info(const APHYTOOLS::Flavour& fl)
    {
      int dof = 1;
      if(fl.isfermion())                                 { dof = 2;p_type='h';};
      if(fl.isvector() &&  AMATOOLS::IsZero(fl.mass()))  { dof = 2;p_type='c';}
      if(fl.isvector() && !AMATOOLS::IsZero(fl.mass()))  {

#ifdef Explicit_Pols
	dof=3;
#else
        dof=1;
#endif
	p_type='c';
      }
      init((int)dof);
      int tf[3]  = {mt::p_m, mt::p_p, mt::p_l };
      for(int j=0;j<dof;j++){
	type[j]   = tf[j];
	factor[j] = 1.;
      }
    }

    ~Pol_Info(){if(type) delete[] type;if(factor)delete[] factor;}
    //Polarisation type: 'c'=circular 'l'=linear 'h'=helicity
    char p_type;
    //direction for linear polarisation
    double angle;
    int      num;
    int    * type;
    double * factor;

    void init(int i)  { num=i; type=new int[num]; factor=new double[num]; }
  };
}
#endif
