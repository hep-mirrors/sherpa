#ifndef String_Generator_h
#define String_Generator_h

#include "Kabbala.H"
#include "MyComplex.H"
#include "Pfunc.H"
#include "String_Tree.H"
#include "Basic_Func.H"
#include <vector>

namespace AMEGIC {
  class Basic_Sfuncs;
  class Values;

  class ZXlist {
  public:
    int     narg;
    int*    arg;
    Kabbala value;
    int     zlist;
    sknot*  sk;
    int     on;
    ZXlist()  {arg = 0;sk = 0;on = 1;narg = 0;zlist = -1;}
    ZXlist(const ZXlist& copy) {arg = 0;*this = copy;}
    ~ZXlist() {if (arg!=0) delete[] arg;} 

    ZXlist& operator=(const ZXlist& copy) {
      if (this!=&copy) {
	narg  = copy.narg;
	value = copy.value;
	zlist = copy.zlist;
	on    = copy.on;
	sk    = copy.sk;
	if (arg!=0) delete[] arg;
	arg = new int[narg];
	for (short int i=0;i<narg;i++) arg[i] = copy.arg[i];
      }
      return *this;
    }
  };

  class Virtual_String_Generator {
  public:
    virtual ~Virtual_String_Generator() {}
    virtual void Reset() {}
    virtual int Get_Cnumber(Complex) {
      std::cerr<<"Virtual Method called !"<<std::endl;return 0;}
#ifdef Kabbala_on
    virtual Kabbala Get_CZnumber(Complex,std::string) {
      std::cerr<<"Virtual Method called !"<<std::endl;return Kabbala(std::string(""),Complex(0.,0.));}
    virtual Kabbala Get_Znumber(int*,Complex*,Complex) { 
      std::cerr<<"Virtual Method called !"<<std::endl;return Kabbala(std::string(""),Complex(0.,0.));}
    virtual Kabbala Get_Xnumber(int*,Complex*,Complex) {
      std::cerr<<"Virtual Method called !"<<std::endl;return Kabbala(std::string(""),Complex(0.,0.));} 
    virtual Kabbala Get_Ynumber(int*,Complex*,Complex) {
      std::cerr<<"Virtual Method called !"<<std::endl;return Kabbala(std::string(""),Complex(0.,0.));} 
    virtual Kabbala Get_Pnumber(Pfunc*,int) {   
      std::cerr<<"Virtual Method called !"<<std::endl;return Kabbala(std::string(""),Complex(0.,0.));} 
    virtual Kabbala Get_Massnumber(int,APHYTOOLS::Flavour,Complex) {
      std::cerr<<"Virtual Method called !"<<std::endl;return Kabbala(std::string(""),Complex(0.,0.));}  
    virtual Kabbala Get_Enumber(Complex) {
      std::cerr<<"Virtual Method called !"<<std::endl;return Kabbala(std::string(""),Complex(0.,0.));}
    virtual Kabbala Get_Snumber(const int,const int,Complex) {     
      std::cerr<<"Virtual Method called !"<<std::endl;return Kabbala(std::string(""),Complex(0.,0.));}
  virtual Kabbala Get_Scplxnumber(const int,const int,Complex)     
    {std::cerr<<"Virtual Method called !"<<std::endl;return Kabbala(std::string(""),Complex(0.,0.));}
#else
    virtual Kabbala Get_CZnumber(Complex,std::string) {
      std::cerr<<"Virtual Method called !"<<std::endl;return Complex(0.,0.);}
    virtual Kabbala Get_Znumber(int*,Complex*,Complex) { 
      std::cerr<<"Virtual Method called !"<<std::endl;return Complex(0.,0.);}
    virtual Kabbala Get_Xnumber(int*,Complex*,Complex) {
      std::cerr<<"Virtual Method called !"<<std::endl;return Complex(0.,0.);} 
    virtual Kabbala Get_Ynumber(int*,Complex*,Complex) {
      std::cerr<<"Virtual Method called !"<<std::endl;return Complex(0.,0.);} 
    virtual Kabbala Get_Pnumber(Pfunc*,int) {   
      std::cerr<<"Virtual Method called !"<<std::endl;return Complex(0.,0.);} 
    virtual Kabbala Get_Massnumber(int,APHYTOOLS::Flavour,Complex) {
      std::cerr<<"Virtual Method called !"<<std::endl;return Complex(0.,0.);}  
    virtual Kabbala Get_Enumber(Complex) {
      std::cerr<<"Virtual Method called !"<<std::endl;return Complex(0.,0.);}
    virtual Kabbala Get_Snumber(const int,const int,Complex) {     
      std::cerr<<"Virtual Method called !"<<std::endl;return Complex(0.,0.);}
  virtual Kabbala Get_Scplxnumber(const int,const int,Complex)     
#endif    
    virtual void     Calculate(Values*) {
      std::cerr<<"Virtual Method called !"<<std::endl;}
    virtual int      Massless(int) {
      std::cerr<<"Virtual Method called !"<<std::endl;
      return 0;
    }
    virtual Kabbala* Get_Kabbala(const std::string&) {
      std::cerr<<"Virtual Method called !"<<std::endl;return 0;}
    virtual inline int ZX_Max_Number() {
      std::cerr<<"Virtual Method called !"<<std::endl;return 0;}
    virtual inline std::vector<int>* GetFlavours() {
      std::cerr<<"Virtual Method called !"<<std::endl;return 0;}
    virtual inline ZXlist* Get_ZXl(const int&) {
      std::cerr<<"Virtual Method called !"<<std::endl;return 0;}
    virtual inline int Coupl_Max_Number() {
      std::cerr<<"Virtual Method called !"<<std::endl;return 0;}
    virtual inline Complex Get_Coupl(const int&) {
      std::cerr<<"Virtual Method called !"<<std::endl;return Complex(0.,0.);}
    virtual inline void SetOn(const int&) {
      std::cerr<<"Virtual Method called !"<<std::endl;}
    virtual inline void SetOn(const std::string&) {
      std::cerr<<"Virtual Method called !"<<std::endl;}
    virtual inline void SetOff(const int&) {
      std::cerr<<"Virtual Method called !"<<std::endl;}
    virtual inline int couplingsLength() {
      std::cerr<<"Virtual Method called !"<<std::endl;
      return 0;
    }
    virtual inline int zxlLength() {
      std::cerr<<"Virtual Method called !"<<std::endl;
      return 0;
    }
  };
  


  class String_Generator : public Virtual_String_Generator,
			   public Basic_Yfunc, 
			   public Basic_Zfunc,
			   public Basic_Xfunc,
			   public Basic_Mfunc,
			   public Basic_Vfunc,
			   public Basic_Pfunc,
			   public Basic_MassTermfunc {
    std::vector<ZXlist>  zxl;
    std::vector<Complex> couplings;
    std::vector<int>     flavours;
    String_Tree stree;
    int ZX_Count(); 
    int Z_Number(int*,int*);
    int X_Number(int*,int*);
    int Y_Number(int*,int*);
    int H_Number(int*,int*);
    int ZXY_Number(int,int,int*,int,int*);
    int GetNumber(int,Complex);
    void Print();
  public:
    String_Generator(Basic_Sfuncs* _BS) : 
      Basic_Func(this,_BS), 
      Basic_Yfunc(this,_BS), 
      Basic_Zfunc(this,_BS), 
      Basic_Xfunc(this,_BS), 
      Basic_Mfunc(this,_BS), 
      Basic_Vfunc(this,_BS), 
      Basic_Pfunc(this,_BS),
      Basic_MassTermfunc(this,_BS) {Reset();}
    ~String_Generator();
    void Reset();
    int Z_Count(); 
    int X_Count(); 
    int E_Count(); 
    
    int Get_Cnumber(Complex);
    int Get_Fnumber(int);
    
    Kabbala Number(int,Complex);
    Kabbala Get_CZnumber(Complex,std::string);
    Kabbala Get_Znumber(int*,Complex*,Complex);
    Kabbala Get_Xnumber(int*,Complex*,Complex); 
    Kabbala Get_Ynumber(int*,Complex*,Complex); 
    Kabbala Get_Pnumber(Pfunc*,int);
    Kabbala Get_Massnumber(int,APHYTOOLS::Flavour,Complex);
    Kabbala Get_Enumber(Complex); 
    Kabbala Get_Snumber(int,int,Complex); 
    Kabbala Get_Scplxnumber(const int,const int,Complex); 
    void Calculate(Values*); 
    int  Massless(int); 
    Kabbala* Get_Kabbala(const std::string&);
    inline int ZX_Max_Number()             { return zxl.size(); }
    inline std::vector<int>* GetFlavours() {return &flavours;}
    inline ZXlist* Get_ZXl(const int& i)   { return &zxl[i]; }
    inline int Coupl_Max_Number()          { return couplings.size(); }
    inline Complex Get_Coupl(const int& i) { return couplings[i]; }
    inline void SetOff(const int& i)       { zxl[i].on = 0; }
    inline void SetOn(const int& i)        { zxl[i].on = 1; }
    void   SetOn(const std::string&);

    inline int zxlLength()                 { return zxl.size(); }
    inline ZXlist zxlist(int i)            { return zxl[i]; }
    inline int couplingsLength()           { return couplings.size(); }
    inline Complex cplist(int i)           { return couplings[i]; } 
    inline int FlavoursLength()            { return flavours.size(); }
    inline int Fllist(int i)               { return flavours[i]; }  
    inline String_Tree Stree()             { return stree; }
  };
  
  class No_String_Generator : public Virtual_String_Generator {
  public:
    void Reset() {};
    int Get_Cnumber(Complex)                         {return 0;}
#ifdef Kabbala_on
    Kabbala Get_CZnumber(Complex value,std::string)       {return Kabbala(std::string(""),value);}
    Kabbala Get_Znumber(int*,Complex*,Complex value) {return Kabbala(std::string(""),value);}
    Kabbala Get_Xnumber(int*,Complex*,Complex value) {return Kabbala(std::string(""),value);} 
    Kabbala Get_Ynumber(int*,Complex*,Complex value) {return Kabbala(std::string(""),value);} 
    Kabbala Get_Pnumber(Pfunc* pl,int)               {return Kabbala(std::string(""),pl->value);} 
    Kabbala Get_Massnumber(int,APHYTOOLS::Flavour,Complex value){return Kabbala(std::string(""),value);}
    Kabbala Get_Enumber(Complex value)               {return Kabbala(std::string(""),value);}
    Kabbala Get_Snumber(int,int,Complex value)     {return Kabbala(std::string(""),value);}
    Kabbala Get_Scplxnumber(const int,const int,Complex value)     {return Kabbala(std::string(""),value);}
#else
    Kabbala Get_CZnumber(Complex value,std::string)       {return value;}
    Kabbala Get_Znumber(int*,Complex*,Complex value) {return value;}
    Kabbala Get_Xnumber(int*,Complex*,Complex value) {return value;} 
    Kabbala Get_Ynumber(int*,Complex*,Complex value) {return value;} 
    Kabbala Get_Pnumber(Pfunc* pl,int)               {return pl->value;} 
    Kabbala Get_Massnumber(int,APHYTOOLS::Flavour,Complex value){return value;}
    Kabbala Get_Enumber(Complex value)               {return value;}
    Kabbala Get_Scplxnumber(const int&,const int&,Complex value)     {return value;}
#endif    
    void     Calculate(Values*)                       {}  
    int      Massless(int)                            { return 0; }  
    Kabbala* Get_Kabbala(const std::string&)          { return 0;}
    inline   int zxlLength()                          { return 0; }
    inline   int couplingsLength()                    { return 0; }
  };
}
#include "Values.H"

#endif









