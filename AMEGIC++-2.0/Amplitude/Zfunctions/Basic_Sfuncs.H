
#ifndef Basic_Sfuncs_h
#define Basic_Sfuncs_h

#include "Flavour.H"
#include "Vector.H"
#include "MyComplex.H"
#include "Pfunc.H"
#include "MathTools.H"
#include "Pol_Info.H"
#include <list>
 
namespace AMEGIC {

  class Momfunc {
  public:
    int*      arg;
    int       argnum;
    AMATOOLS::Vec4D     mom,mom_img;
    double    mass,angle;
    Complex   cplxmass2;
    mt::momtype type;
    int       on;
    Momfunc() {arg = 0;argnum = 0;}
    Momfunc(const Momfunc& m) {
      arg = 0;
      argnum = 0;
      *this = m;
    }

    ~Momfunc() {if (arg) delete[] arg;}
      

    Momfunc& operator=(const Momfunc& m) {
      if (this!=&m) {
	argnum = m.argnum;
	if (arg) delete[] arg;
	if (argnum>0) {
	  arg = new int[argnum];
	  for (short int i=0;i<argnum;i++) arg[i] = m.arg[i];
	}
	else std::cout<<"Number of arguments is zero in Momfunc operator=()"<<std::endl;
	mom     = m.mom;
	mom_img = m.mom_img;
	mass    = m.mass;
	cplxmass2    = m.cplxmass2;
	angle   = m.angle;
	type    = m.type;
	on      = m.on;
      }
    }
};

class Basic_Sfuncs {
  std::vector<Momfunc> Momlist;
  APHYTOOLS::Flavour* fl;
  int Nmom,nvec;
  int* b;
  double** delta4x4;
  Complex* _eta;
  Complex* _mu;
  AMATOOLS::Vec4D k0,k1;
  Complex** _S0;
  Complex** _S1;
  int** calc_st;
  int  momcount;
  AMATOOLS::Vec4D* p;
  int k0_n;
  int  Initialize_Momlist();
  void Calc_Momlist();
  void calcS(int i, int j);
  AMATOOLS::Vec4D Mom_dir(int i);
public:
  Basic_Sfuncs(int,int, APHYTOOLS::Flavour*,int*);
  ~Basic_Sfuncs();
  double N(int,int);
  int epsilon(int,int,int,int);
  int setS(AMATOOLS::Vec4D*);
  void ResetS_GT(double);
  void set_k0(int);
  inline Complex S0(int i,int j) {i=AMATOOLS::iabs(i);j=AMATOOLS::iabs(j);if(!calc_st[i][j])calcS(i,j); return _S0[i][j];}
  inline Complex S1(int i,int j) {i=AMATOOLS::iabs(i);j=AMATOOLS::iabs(j);if(!calc_st[i][j])calcS(i,j); return _S1[i][j];}
  inline Complex mu(int i) {return (i>0) ? _mu[i] : -_mu[AMATOOLS::iabs(i)];}
  inline Complex eta(int i) {return _eta[AMATOOLS::iabs(i)];}
  inline int    GetNmomenta() {return Nmom;}
  inline int    Sign(int i)   {return b[AMATOOLS::iabs(i)];}
  inline APHYTOOLS::Flavour Flav(int i)   {return fl[AMATOOLS::iabs(i)];}
  inline AMATOOLS::Vec4D  Momentum(int i){return Momlist[i].mom;}
  inline AMATOOLS::Vec4D  Momentum_img(int i){return Momlist[i].mom_img;}
  int    Build_Momlist(std::list<Pfunc*>&); 
  int    Build_TensorPolarisations(int); 
  int    Build_Polarisations(int,char,double angle=0.);
  int    Build_Polarisations(int,APHYTOOLS::Flavour);
  void   Print_Momlist();
  int    Get_Mom_Number(Pfunc*);
  void   PropPolarisation(int,std::list<Pfunc*>&,std::vector<int>&);
  int    Get_Pol_Number(int,int,double,int check=0);
  void   Initialize();
  inline bool iscplx(int i) {return (Momlist[i].type==mt::p_p||Momlist[i].type==mt::p_m || Momlist[i].type==mt::p_l || Momlist[i].type==mt::p_s|| Momlist[i].type==mt::p_si);} 
};
}
#endif
