#ifndef Pol_Info_h
#define Pol_Info_h

#include "Flavour.H"
#include "MathTools.H"
#include "Message.H"
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
		 cmprop=12,  //center of momentum propagator
		 p_p  = 1,   //polarisation +1 (circular)
		 p_m  =-1,   //             -1 (circular)
		 p_l  = 0,   //             longitudinal
		 p_s  = 2,   //             scalar
		 p_si = 7,   //             massless scalar
		 p_l0 = 3,   //             linear with direction
		 p_l1 = 4,   //             linear, perpendicular to p_l0
	  	 p_lh = 5,   //             horizontaly (linear)
		 p_lp = 6,    //             perpendicular (linear)
		 p_none = 8,

		 p_t1 = 101,  //dummies for graviton polarization tensors
		 p_t2 = 102,
		 p_t3 = 103,
		 p_t4 = 104,
		 p_t5 = 105
  };
};

  class Pol_Info{
  public:
    Pol_Info(const Pol_Info & p) {
      num=p.num;
      p_type=p.p_type;
      angle=p.angle;
      type=new int[num];
      factor=new double[num];
      for(int i=0;i<num;i++){
	type[i]=p.type[i];
	factor[i]=p.factor[i];
      }
    }
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
      if(fl.IsFermion())                                 { dof = 2;p_type='h';};
      if(fl.IsVector() &&  AMATOOLS::IsZero(fl.Mass()))  { dof = 2;p_type='c';}
      if(fl.IsVector() && !AMATOOLS::IsZero(fl.Mass()))  {

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
    //Polarisation type: 'c'=circular 'l'=linear 'h'=helicity 't'=tensor
    char p_type;
    //direction for linear polarisation
    double angle;
    int      num;
    int    * type;
    double * factor;
    void init(int i){num=i;type=new int[num];factor=new double[num];}
  };

  class tensor_struc{                //defines the spin2 tensors
  public:
    void Get_Pol_Combos(int num, vector<vector<int> >* pol, vector<int>* sign){
      pol->clear();
      sign->clear();
      vector<int> cc;cc.push_back(8);cc.push_back(8); 
      sign->push_back(1);
      switch(num){
      case mt::p_t1:
	cc[0]=mt::p_p;cc[1]=mt::p_p;
	pol->push_back(cc);
	break;
      case mt::p_t2:
	cc[0]=mt::p_p;cc[1]=mt::p_l;
	pol->push_back(cc);
	break;
      case mt::p_t3:
	cc[0]=mt::p_p;cc[1]=mt::p_m;
	pol->push_back(cc);
	cc[0]=mt::p_l;cc[1]=mt::p_l;
	pol->push_back(cc);
	sign->push_back(-1);
	break;
      case mt::p_t4:
	cc[0]=mt::p_m;cc[1]=mt::p_l;
	pol->push_back(cc);
	break;
      case mt::p_t5:
	cc[0]=mt::p_m;cc[1]=mt::p_m;
	pol->push_back(cc);
	break;
      default: AORGTOOLS::msg.Error()<<"Invalid tensor type: "<<num<<endl;abort();
      }
    }
    double Get_Tfactor(int num){
      if(num<mt::p_t1) return 1.;
      switch(num){
      case mt::p_t1:return 1.;
      case mt::p_t2:return 2.;
      case mt::p_t3:return 2./3.;
      case mt::p_t4:return 2.;
      case mt::p_t5:return 1.;
      }
    }
  };
}
#endif
