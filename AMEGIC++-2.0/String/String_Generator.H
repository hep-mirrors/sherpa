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
    int       narg, * arg;
    Kabbala   value;
    int       zlist;
    sknot   * sk;
    int       on;

    ZXlist()                    { arg = 0;sk = 0;on = 1;narg = 0;zlist = -1; }
    ZXlist(const ZXlist & copy) { arg = 0;*this = copy; }
    ~ZXlist()                   { if (arg!=0) delete[] arg; } 

    ZXlist& operator=(const ZXlist& copy);
  };

  class Virtual_String_Generator {
  public:
    virtual ~Virtual_String_Generator() {}
    virtual void      Reset(int = 0)    {}
    virtual void      ReplaceZXlist(Virtual_String_Generator*) {}   
    virtual void      ReStore()         {}
    virtual int       GetCnumber(Complex)                          = 0;
    virtual Kabbala   GetCZnumber(Complex,std::string)             = 0;
    virtual Kabbala   GetZnumber(int*,Complex*,Complex)            = 0; 
    virtual Kabbala   GetXnumber(int*,Complex*,Complex)            = 0;
    virtual Kabbala   GetYnumber(int*,Complex*,Complex)            = 0;
    virtual Kabbala   GetPnumber(Pfunc *,int)                      = 0;   
    virtual Kabbala   GetMassnumber(int,ATOOLS::Flavour,Complex)   = 0;
    virtual Kabbala   GetEnumber(Complex)                          = 0;
    virtual Kabbala   GetSnumber(const int,const int,Complex)      = 0;     
    virtual Kabbala   GetScplxnumber(const int,const int,Complex)  = 0;      
    virtual void      Calculate(Values*)                           = 0;
    virtual int       Massless(int)                                = 0;
    virtual Kabbala * GetKabbala(const std::string&)               = 0;
    virtual int       ZXMaxNumber()                                = 0;
    virtual std::vector<int>* GetFlavours()                        = 0;
    virtual ZXlist*   GetZXl(const int&)                           = 0;
    virtual Complex   GetCoupling(const int&)                      = 0;
    virtual std::vector<Complex>* GetCouplList()                   = 0;
    virtual std::vector<ZXlist>*  GetZXlist()                      = 0;
    virtual void      SetOn(const int&)                            = 0;
    virtual void      SetOn(const std::string&)                    = 0;
    virtual void      SetOff(const int&)                           = 0;
    virtual int       NumberOfCouplings()                          = 0;
    virtual int       NumberOfzxl()                                = 0;
  };
  


  class String_Generator : public Virtual_String_Generator,
			   public Basic_Yfunc, 
			   public Basic_Zfunc,
			   public Basic_Xfunc,
			   public Basic_Mfunc,
			   public Basic_Vfunc,
			   public Basic_Pfunc,
			   public Basic_MassTermfunc {
    std::vector<ZXlist>  * p_zxl,       * p_zxlsave;
    std::vector<Complex> * p_couplings, * p_couplingssave;
    std::vector<int>     * p_flavours,  * p_flavourssave;
    String_Tree            m_stree,m_sthelp;
    bool                   m_copied;
    int  ZXCount(); 
    int  ZNumber(int *,int *);
    int  XNumber(int *,int *);
    int  YNumber(int *,int *);
    int  HNumber(int *,int *);
    int  ZXYNumber(int,int,int *,int,int *);
    int  GetNumber(int,Complex);
    void Print();
  public:
    String_Generator(Basic_Sfuncs *);
    ~String_Generator();
    void Reset(int = 0);
    void ReplaceZXlist(Virtual_String_Generator*);
    void ReStore();
    int ZCount(); 
    int XCount(); 
    int ECount(); 
    
    int GetCnumber(Complex);
    int GetFnumber(int);
    
    Kabbala     Number(int,Complex);
    Kabbala     GetCZnumber(Complex,std::string);
    Kabbala     GetZnumber(int *,Complex *,Complex);
    Kabbala     GetXnumber(int *,Complex *,Complex); 
    Kabbala     GetYnumber(int *,Complex *,Complex); 
    Kabbala     GetPnumber(Pfunc *,int);
    Kabbala     GetMassnumber(int,ATOOLS::Flavour,Complex);
    Kabbala     GetEnumber(Complex); 
    Kabbala     GetSnumber(int,int,Complex); 
    Kabbala     GetScplxnumber(const int,const int,Complex); 
    void        Calculate(Values*); 
    int         Massless(int); 
    Kabbala   * GetKabbala(const std::string&);
    int         ZXMaxNumber()               { return p_zxl->size(); }
    std::vector<int>* GetFlavours()         { return p_flavours;}
    ZXlist    * GetZXl(const int & i)       { return &(*p_zxl)[i]; }
    Complex     GetCoupling(const int & i)  { return (*p_couplings)[i]; }
    std::vector<Complex>* GetCouplList()    { return p_couplings; }
    std::vector<ZXlist>*  GetZXlist()       { return p_zxl; }
    void        SetOff(const int & i)       { (*p_zxl)[i].on = 0; }
    void        SetOn(const int & i)        { (*p_zxl)[i].on = 1; }
    void        SetOn(const std::string &);

    int         NumberOfzxl()               { return p_zxl->size(); }
    ZXlist      zxlist(int i)               { return (*p_zxl)[i]; }
    int         NumberOfCouplings()         { return p_couplings->size(); }
    Complex     cplist(int i)               { return (*p_couplings)[i]; } 
    int         NumberOfFlavours()          { return p_flavours->size(); }
    int         GetFlavour(int i)           { return (*p_flavours)[i]; }  
    String_Tree GetStree()                  { return m_stree; }
  };
  
  class No_String_Generator : public Virtual_String_Generator {
  public:
    void Reset()                                                {};
    int       GetCnumber(Complex)                               { return 0;}
    Kabbala   GetCZnumber(Complex value,std::string)            { return Kabbala(std::string(""),value);}
    Kabbala   GetZnumber(int*,Complex*,Complex value)           { return Kabbala(std::string(""),value);}
    Kabbala   GetXnumber(int*,Complex*,Complex value)           { return Kabbala(std::string(""),value);} 
    Kabbala   GetYnumber(int*,Complex*,Complex value)           { return Kabbala(std::string(""),value);} 
    Kabbala   GetPnumber(Pfunc * pl,int)                        { return Kabbala(std::string(""),pl->value);} 
    Kabbala   GetMassnumber(int,ATOOLS::Flavour,Complex value)  { return Kabbala(std::string(""),value);}
    Kabbala   GetEnumber(Complex value)                         { return Kabbala(std::string(""),value);}
    Kabbala   GetSnumber(int,int,Complex value)                 { return Kabbala(std::string(""),value);}
    Kabbala   GetScplxnumber(const int,const int,Complex value) { return Kabbala(std::string(""),value);}
    int       Massless(int)                                     { return 0; }  
    Kabbala * GetKabbala(const std::string&)                    { return 0;}
    int       NumberOfzxl()                                     { return 0; }
    int       NumberOfCouplings()                               { return 0; }
    int       ZXMaxNumber()                                     { return 0; }   
    std::vector<int>* GetFlavours()                             { return NULL; }
    ZXlist*   GetZXl(const int&)                                { return NULL; }
    Complex   GetCoupling(const int&)                           { return Complex(0.,0.); }
    std::vector<Complex>* GetCouplList()                        { return NULL; }
    std::vector<ZXlist>*  GetZXlist()                           { return NULL; }
    void      Calculate(Values*)                                {}  
    void      SetOn(const int&)                                 {}  
    void      SetOn(const std::string&)                         {}  
    void      SetOff(const int&)                                {}  
  };

}
#include "Values.H"

#endif









