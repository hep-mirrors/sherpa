#ifndef Vertex_h
#define Vertex_h

#include "MyComplex.H"
#include "Flavour.H"
#include "Color_Function.H"
#include "Lorentz_Function.H"
//#include "Kabbala.H"


namespace AMEGIC {

  class Single_Vertex {
  public:
    APHYTOOLS::Flavour   in[4];
    //number of legs
    int       nleg;
    short int on;
    Complex   cpl[4];
    std::string    Str;
    Color_Function*   Color;
    //number of Color_Functions
    int               ncf;
    Lorentz_Function* Lorentz;
    //number of Lorentz_Functions
    int               nlf;
    
    Single_Vertex() {
      nleg    = 3;
      Color   = new Color_Function;
      Lorentz = new Lorentz_Function;
    }

    ~Single_Vertex() {
      delete Color;
      delete Lorentz;
    }

    Single_Vertex& operator=(const Single_Vertex& v) {
      if (this!=&v) {
	for (short int i=0;i<4;i++) in[i]  = v.in[i];
	for (short int j=0;j<4;j++) cpl[j] = v.cpl[j];
	
	nleg = v.nleg;
	Str  = v.Str;
	on   = v.on;
	ncf  = v.ncf;
	nlf  = v.nlf;
	*Color   = *v.Color;
	*Lorentz = *v.Lorentz;
      }
    }
  };
  
  class Vertex {
    Single_Vertex* v;
    Single_Vertex* v4;
    void           CheckEqual(APHYTOOLS::Flavour**,short int&);
    int            nvertex;
    int            n4vertex;
    void           GenerateVertex();
    int            CheckExistence(Single_Vertex&);
    int            FermionRule(Single_Vertex& probe);
    int            SetVertex(Single_Vertex&,Single_Vertex&,int,int,int,int a=99);

    void           ColorExchange(Color_Function*,int,int,int);
    void           LorentzExchange(Lorentz_Function*,int,int,int,int);
    void           Conjugate(Color_Function*);
  public:
    int on;
    Vertex();
    void Print();
    ~Vertex();
    inline Single_Vertex* GetVertex(int n)   {return &v[n];}
    // use [] for 3 leg vertices
    inline Single_Vertex* operator[] (int n) {return &v[n];}  
    inline const Single_Vertex* operator[] (int n) const {return &v[n];}
    // use () for 4 leg vertices
    inline Single_Vertex* operator() (int n) {return &v4[n];}  
    inline const Single_Vertex* operator() (int n) const {return &v4[n];}
    inline const int Max_Number() const {return nvertex;}
    inline const int Max_Number4() const {return n4vertex;}
    void TexOutput();
    void AddVertex(Single_Vertex*);
    int  FindVertex(Single_Vertex*);
  };
  
  std::ostream& operator<<(std::ostream&, const Single_Vertex&);
  std::ostream& operator<<(std::ostream&, const Vertex&);
}
#endif








