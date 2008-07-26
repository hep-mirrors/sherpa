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
    ATOOLS::Kabbala   value;
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
    virtual ATOOLS::Kabbala   GetCZnumber(Complex,std::string)             = 0;
    virtual ATOOLS::Kabbala   GetZnumber(int*,Complex*,Complex)            = 0; 
    virtual ATOOLS::Kabbala   GetXnumber(int*,Complex*,Complex)            = 0;
    virtual ATOOLS::Kabbala   GetYnumber(int*,Complex*,Complex)            = 0;
    virtual ATOOLS::Kabbala   GetPnumber(Pfunc *,int)                      = 0;   
    virtual ATOOLS::Kabbala   GetMassnumber(int,ATOOLS::Flavour,Complex)   = 0;
    virtual ATOOLS::Kabbala   GetMnumber(ATOOLS::Flavour,Complex)          = 0;
    virtual ATOOLS::Kabbala   GetEnumber(Complex)                          = 0;
    virtual ATOOLS::Kabbala   GetSnumber(const int,const int,Complex)      = 0;     
    virtual ATOOLS::Kabbala   GetScplxnumber(const int,const int,Complex)  = 0;      
    virtual ATOOLS::Kabbala   GetEpsnumber(int*,int,Complex)               = 0;     
    virtual ATOOLS::Kabbala   GetSFnumber(Complex,int)                     = 0;
    virtual void      Calculate(Values*)                           = 0;
    virtual int       Massless(int)                                = 0;
    virtual ATOOLS::Kabbala * GetKabbala(const std::string&)               = 0;
    virtual int       ZXMaxNumber()                                = 0;
    virtual std::vector<long int>* GetFlavours()                   = 0;
    virtual ZXlist*   GetZXl(const int&)                           = 0;
    virtual Complex   GetCoupling(const int&)                      = 0;
    virtual std::vector<Complex>* GetCouplList()                   = 0;
    virtual std::vector<ZXlist>*  GetZXlist()                      = 0;
    virtual void      SetOn(const int&)                            = 0;
    virtual void      SetOn(const std::string&)                    = 0;
    virtual void      SetOff(const int&)                           = 0;
    virtual int       NumberOfCouplings()                          = 0;
    virtual int       NumberOfzxl()                                = 0;
    virtual int       UsesFunction(size_t i) {return 0;}
    virtual void      SetCouplings(std::map<std::string,Complex>*) = 0;
    virtual std::map<std::string,Complex>* GetCouplings()          = 0;
    virtual void      WriteCouplings(std::ofstream&)               = 0;
    virtual int       ReadCouplings(std::ifstream&)                = 0;
    virtual void      UpdateCouplings(std::map<std::string,Complex> & cmap) = 0;
  };
  


  class String_Generator : public Virtual_String_Generator,
			   public Basic_Yfunc, 
			   public Basic_Zfunc,
			   public Basic_Xfunc,
			   public Basic_Mfunc,
			   public Basic_Vfunc,
			   public Basic_Pfunc,
			   public Basic_MassTermfunc,
			   public Basic_Epsilonfunc,
			   public Unitarityfunc {
    std::vector<ZXlist>  * p_zxl,       * p_zxlsave;
    std::vector<Complex> * p_couplings, * p_couplingssave;
    std::vector<long int>* p_flavours,  * p_flavourssave;
    std::vector<int>       m_zuse;
    String_Tree            m_stree,m_sthelp;
    bool                   m_copied;
    std::map<std::string,Complex>* p_couplmap;
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
    int GetFnumber(long int);
    
    ATOOLS::Kabbala     Number(int,Complex);
    ATOOLS::Kabbala     GetCZnumber(Complex,std::string);
    ATOOLS::Kabbala     GetZnumber(int *,Complex *,Complex);
    ATOOLS::Kabbala     GetXnumber(int *,Complex *,Complex); 
    ATOOLS::Kabbala     GetYnumber(int *,Complex *,Complex); 
    ATOOLS::Kabbala     GetPnumber(Pfunc *,int);
    ATOOLS::Kabbala     GetMassnumber(int,ATOOLS::Flavour,Complex);
    ATOOLS::Kabbala     GetMnumber(ATOOLS::Flavour,Complex);
    ATOOLS::Kabbala     GetEnumber(Complex); 
    ATOOLS::Kabbala     GetSnumber(int,int,Complex); 
    ATOOLS::Kabbala     GetScplxnumber(const int,const int,Complex); 
    ATOOLS::Kabbala     GetEpsnumber(int*,int,Complex);
    ATOOLS::Kabbala     GetSFnumber(Complex,int);
    void        Calculate(Values*); 
    int         Massless(int); 
    ATOOLS::Kabbala   * GetKabbala(const std::string&);
    int         ZXMaxNumber()               { return p_zxl->size(); }
    std::vector<long int>* GetFlavours()    { return p_flavours;}
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
    int         UsesFunction(size_t i)      { return m_zuse[i]; }
    void        SetCouplings(std::map<std::string,Complex>* cplm) { p_couplmap=cplm; }
    std::map<std::string,Complex>* GetCouplings()                 { return p_couplmap; }
    void        WriteCouplings(std::ofstream&);
    int         ReadCouplings(std::ifstream&);
    void        UpdateCouplings(std::map<std::string,Complex> & cmap);
  };
  
  class No_String_Generator : public Virtual_String_Generator {
  public:
    void Reset()                                                {};
    int       GetCnumber(Complex)                               { return 0;}
    ATOOLS::Kabbala   GetCZnumber(Complex value,std::string)            { return ATOOLS::Kabbala(std::string(""),value);}
    ATOOLS::Kabbala   GetZnumber(int*,Complex*,Complex value)           { return ATOOLS::Kabbala(std::string(""),value);}
    ATOOLS::Kabbala   GetXnumber(int*,Complex*,Complex value)           { return ATOOLS::Kabbala(std::string(""),value);} 
    ATOOLS::Kabbala   GetYnumber(int*,Complex*,Complex value)           { return ATOOLS::Kabbala(std::string(""),value);} 
    ATOOLS::Kabbala   GetPnumber(Pfunc * pl,int)                        { return ATOOLS::Kabbala(std::string(""),pl->value);} 
    ATOOLS::Kabbala   GetMassnumber(int,ATOOLS::Flavour,Complex value)  { return ATOOLS::Kabbala(std::string(""),value);}
    ATOOLS::Kabbala   GetMnumber(ATOOLS::Flavour,Complex value)         { return ATOOLS::Kabbala(std::string(""),value);}
    ATOOLS::Kabbala   GetEnumber(Complex value)                         { return ATOOLS::Kabbala(std::string(""),value);}
    ATOOLS::Kabbala   GetSnumber(int,int,Complex value)                 { return ATOOLS::Kabbala(std::string(""),value);}
    ATOOLS::Kabbala   GetScplxnumber(const int,const int,Complex value) { return ATOOLS::Kabbala(std::string(""),value);}
    ATOOLS::Kabbala   GetEpsnumber(int*,int,Complex value)              { return ATOOLS::Kabbala(std::string(""),value);}
    ATOOLS::Kabbala   GetSFnumber(Complex value ,int)                   { return ATOOLS::Kabbala(std::string(""),value);}
    int       Massless(int)                                     { return 0; }  
    ATOOLS::Kabbala * GetKabbala(const std::string&)                    { return 0;}
    int       NumberOfzxl()                                     { return 0; }
    int       NumberOfCouplings()                               { return 0; }
    int       ZXMaxNumber()                                     { return 0; }   
    std::vector<long int>* GetFlavours()                        { return NULL; }
    ZXlist*   GetZXl(const int&)                                { return NULL; }
    Complex   GetCoupling(const int&)                           { return Complex(0.,0.); }
    std::vector<Complex>* GetCouplList()                        { return NULL; }
    std::vector<ZXlist>*  GetZXlist()                           { return NULL; }
    void      Calculate(Values*)                                {}  
    void      SetOn(const int&)                                 {}  
    void      SetOn(const std::string&)                         {}  
    void      SetOff(const int&)                                {}  
    void      SetCouplings(std::map<std::string,Complex>*)      {}
    std::map<std::string,Complex>* GetCouplings()               { return NULL; }
    void      WriteCouplings(std::ofstream&)                    {}
    int       ReadCouplings(std::ifstream&)                     { return 0; }
    void      UpdateCouplings(std::map<std::string,Complex> & cmap) {}
  };

}
#include "Values.H"

#endif









