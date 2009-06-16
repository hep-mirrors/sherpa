#include "MODEL/Interaction_Models/Lorentz_Function.H"

using namespace MODEL;
using namespace ATOOLS;

class LF_None: public Lorentz_Function {
public:
  LF_None(): Lorentz_Function("None") {}
  int NofIndex() const { return 0; }
  std::string String(int shortversion) const 
  { return "0"; }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_None());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_None,LFNone_Getter,"None","")

class LF_Gamma: public Lorentz_Function {
public:  
  LF_Gamma(): Lorentz_Function("Gamma") {}
  int NofIndex() const { return 1; }
  std::string String(int shortversion) const 
  { 
    // Gam[0]
    return "Gam["+Str(0)+"]"; 
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Gamma());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Gamma,LFGamma_Getter,"Gamma","")
class LF_Gab: public Lorentz_Function {
public:  
  LF_Gab(): Lorentz_Function("Gab") {}
  int NofIndex() const { return 2; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation(1,0,1);
    AddPermutation(1,1,0);  
  }
  std::string String(int shortversion) const 
  { 
    // G[0,1]
    return "G["+Str(0)+","+Str(1)+"]"; 
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Gab());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Gab,LFGab_Getter,"Gab","")
class LF_Gauge3: public Lorentz_Function {
public:  
  LF_Gauge3(): Lorentz_Function("Gauge3") {}
  int NofIndex() const { return 3; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2);
    AddPermutation(-1,0,2,1);  
    AddPermutation(-1,1,0,2);
    AddPermutation(-1,2,1,0);  
    AddPermutation( 1,1,2,0);
    AddPermutation( 1,2,0,1);  
  }
  std::string String(int shortversion) const 
  {
    return "1";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Gauge3());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Gauge3,LFGauge3_Getter,"Gauge3","")
class LF_Gauge4: public Lorentz_Function {
public:  
  LF_Gauge4(): Lorentz_Function("Gauge4") {}
  int NofIndex() const { return 4; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2,3);
    AddPermutation( 1,1,0,2,3);
    AddPermutation( 1,0,1,3,2);
    AddPermutation( 1,1,0,3,2);
    AddPermutation( 1,2,3,1,0);    
    AddPermutation( 1,3,2,1,0);
    AddPermutation( 1,2,3,0,1);
    AddPermutation( 1,3,2,0,1);
  }
  std::string String(int shortversion) const 
  {
    std::string help;
    //(2G(0,1)*G(2,3)-G(0,2)*G(1,3)-G(0,3)*G(1,2))
    help  = std::string("(2*G[")  + Str(0) + std::string(",") + Str(1) + std::string("]*");
    help += std::string("G[")  + Str(2) + std::string(",") + Str(3) + std::string("]-");
    help += std::string("G[")  + Str(0) + std::string(",") + Str(2) + std::string("]*");
    help += std::string("G[")  + Str(1) + std::string(",") + Str(3) + std::string("]-");
    help += std::string("G[")  + Str(0) + std::string(",") + Str(3) + std::string("]*");
    help += std::string("G[")  + Str(1) + std::string(",") + Str(2) + std::string("])");
    return help;
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Gauge4());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Gauge4,LFGauge4_Getter,"Gauge4","")
class LF_Gluon4: public Lorentz_Function {
public:  
  LF_Gluon4(): Lorentz_Function("Gluon4") {}
  int NofIndex() const { return 4; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2,3);
    AddPermutation(-1,2,1,0,3);
    AddPermutation(-1,0,3,2,1);
    AddPermutation( 1,2,3,0,1);
    AddPermutation( 1,1,0,3,2);
    AddPermutation(-1,3,0,1,2);
    AddPermutation(-1,1,2,3,0);
    AddPermutation( 1,3,2,1,0);
  }
  std::string String(int shortversion) const 
  {
    std::string help;
    //G(0,1)*G(2,3)-G(0,3)*G(2,1)
    if (shortversion) {
      help += std::string("G4[") + Str(0) + std::string(",") + 
	Str(1) + std::string(",") + 
	Str(2) + std::string(",") + 
	Str(3) + std::string("]");
    }
    else {
      help  = std::string("(G[")  + Str(0) + std::string(",") + Str(1) + std::string("]*");
      help += std::string("G[")  + Str(2) + std::string(",") + Str(3) + std::string("]-");
      help += std::string("G[")  + Str(0) + std::string(",") + Str(3) + std::string("]*");
      help += std::string("G[")  + Str(2) + std::string(",") + Str(1) + std::string("])");
    }
    return help;
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Gluon4());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Gluon4,LFGluon4_Getter,"Gluon4","")
class LF_SSV: public Lorentz_Function {
public:  
  LF_SSV(): Lorentz_Function("SSV") {}
  int NofIndex() const { return 3; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2);
    AddPermutation(-1,1,0,2);
  }
  std::string String(int shortversion) const 
  {
    //P[0,2]-P[1,2]
    std::string help = std::string("P[") + Str(0) + std::string(",") + Str(2) +std::string("]-"); 
    return help + std::string("P[") + Str(1) + std::string(",") + Str(2) +std::string("]");
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_SSV());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_SSV,LFSSV_Getter,"SSV","")
class LF_SSS: public Lorentz_Function {
public:  
  LF_SSS(): Lorentz_Function("SSS") {}
  int NofIndex() const { return 0; }
  std::string String(int shortversion) const 
  {
    return "1";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_SSS());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_SSS,LFSSS_Getter,"SSS","")
class LF_FFS: public Lorentz_Function {
public:  
  LF_FFS(): Lorentz_Function("FFS") {}
  int NofIndex() const { return 0; }
  std::string String(int shortversion) const 
  {
    return "1";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_FFS());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_FFS,LFFFS_Getter,"FFS","")
class LF_Pol: public Lorentz_Function {
public:  
  LF_Pol(): Lorentz_Function("Pol") {}
  int NofIndex() const { return 1; }
  std::string String(int shortversion) const 
  {
    // Eps[0]
    return "Eps["+Str(0)+"]";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Pol());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Pol,LFPol_Getter,"Pol","")
class LF_VVSS: public Lorentz_Function {
public:  
  LF_VVSS(): Lorentz_Function("VVSS") {}
  int NofIndex() const { return 2; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation(1,0,1);
    AddPermutation(1,1,0);  
  }
  std::string String(int shortversion) const 
  { 
    // G(2V2S)[0,1]
    return "G(2V2S)["+Str(0)+","+Str(1)+"]"; 
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_VVSS());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_VVSS,LFVVSS_Getter,"VVSS","")
class LF_SSSS: public Lorentz_Function {
public:  
  LF_SSSS(): Lorentz_Function("SSSS") {}
  int NofIndex() const { return 0; }
  std::string String(int shortversion) const 
  {
    return "1";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_SSSS());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_SSSS,LFSSSS_Getter,"SSSS","")
class LF_AGauge4: public Lorentz_Function {
public: 
  LF_AGauge4(): Lorentz_Function("AGauge4") {}
  int NofIndex() const { return 4; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2,3);
    AddPermutation( 1,1,0,2,3);
    AddPermutation( 1,0,1,3,2);
    AddPermutation( 1,1,0,3,2);
    AddPermutation( 1,2,3,1,0);    
    AddPermutation( 1,3,2,1,0);
    AddPermutation( 1,2,3,0,1);
    AddPermutation( 1,3,2,0,1);
  }
  std::string String(int shortversion) const 
  {
    return "1";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_AGauge4());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_AGauge4,LFAGauge4_Getter,"AGauge4","")
class LF_AGauge3: public Lorentz_Function {
public:  
  LF_AGauge3(): Lorentz_Function("AGauge3") {}
  int NofIndex() const { return 3; }
  bool CutVectors() { return true; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2);
    AddPermutation(-1,0,2,1);  
    AddPermutation(-1,1,0,2);
    AddPermutation(-1,2,1,0);  
    AddPermutation( 1,1,2,0);
    AddPermutation( 1,2,0,1);  
  }
  std::string String(int shortversion) const 
  {
    // (P[0,2]-P[1,2])*G(0,1)+(P[1,0]-P[2,0])*G(1,2)+(P[2,1]-P[0,1])*G(2,0)
    std::string help;
    if (shortversion) {
      help += std::string("V3[") + Str(0) + std::string(",") + 
	Str(1) + std::string(",") + 
	Str(2) + std::string("]");
    }
    else {
      help  = std::string("(P[") + Str(0) + std::string(",") + Str(2) + std::string("]-");
      help += std::string("P[")  + Str(1) + std::string(",") + Str(2) + std::string("])*");
      help += std::string("G[")  + Str(0) + std::string(",") + Str(1) + std::string("]");
	  
      help += std::string("+");
	  
      help += std::string("(P[") + Str(1) + std::string(",") + Str(0) + std::string("]-");
      help += std::string("P[")  + Str(2) + std::string(",") + Str(0) + std::string("])*");
      help += std::string("G[")  + Str(1) + std::string(",") + Str(2) + std::string("]");
	  
      help += std::string("+");
	  
      help += std::string("(P[") + Str(2) + std::string(",") + Str(1) + std::string("]-");
      help += std::string("P[")  + Str(0) + std::string(",") + Str(1) + std::string("])*");
      help += std::string("G[")  + Str(2) + std::string(",") + Str(0) + std::string("]");
    }
    return help;
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_AGauge3());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_AGauge3,LFAGauge3_Getter,"AGauge3","")

class LF_AZZZ: public Lorentz_Function {
public:  
  LF_AZZZ(): Lorentz_Function("AZZZ") {}
  int NofIndex() const { return 3; }
  bool CutVectors() { return true; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation(1,0,1,2);
    AddPermutation(1,0,2,1);  
    AddPermutation(1,1,0,2);
    AddPermutation(1,2,1,0);  
    AddPermutation(1,1,2,0);
    AddPermutation(1,2,0,1);  
  }
  std::string String(int shortversion) const 
  {
    std::string help;
    help += std::string("AZZZ[") + Str(0) + std::string(",") + 
      Str(1) + std::string(",") + 
      Str(2) + std::string("]");
    return help;
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_AZZZ());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_AZZZ,LFAZZZ_Getter,"AZZZ","")


class LF_AZZG: public Lorentz_Function {
public:  
  LF_AZZG(): Lorentz_Function("AZZG") {}
  int NofIndex() const { return 3; }
  bool CutVectors() { return true; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation(1,0,1,2);
    AddPermutation(1,0,2,1);  
  }
  std::string String(int shortversion) const 
  {
    std::string help;
    help += std::string("AZZG[") + Str(0) + std::string(",") + 
      Str(1) + std::string(",") + 
      Str(2) + std::string("]");
    return help;
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_AZZG());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_AZZG,LFAZZG_Getter,"AZZG","")


class LF_AZGG: public Lorentz_Function {
public:  
  LF_AZGG(): Lorentz_Function("AZGG") {}
  int NofIndex() const { return 3; }
  bool CutVectors() { return true; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation(1,0,1,2);
    AddPermutation(1,2,1,0);  
  }
  std::string String(int shortversion) const 
  {
    std::string help;
    help += std::string("AZGG[") + Str(0) + std::string(",") + 
      Str(1) + std::string(",") + 
      Str(2) + std::string("]");
    return help;
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_AZGG());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_AZGG,LFAZGG_Getter,"AZGG","")

class LF_FFT: public Lorentz_Function {
public: 
  LF_FFT(): Lorentz_Function("FFT") {}
  int NofIndex() const { return 1; }
  std::string String(int shortversion) const 
  {
    return "FFT["+Str(0)+","+Str(1)+"]";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_FFT());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_FFT,LFFFT_Getter,"FFT","")
class LF_VVT: public Lorentz_Function {
public:  
  LF_VVT(): Lorentz_Function("VVT") {}
  int NofIndex() const { return 3; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2);
    AddPermutation( 1,1,0,2);
  }
  std::string String(int shortversion) const 
  {
    return "VVT["+Str(0)+","+Str(1)+","+Str(2)+"]";    
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_VVT());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_VVT,LFVVT_Getter,"VVT","")
class LF_SST: public Lorentz_Function {
public: 
  LF_SST(): Lorentz_Function("SST") {}
  int NofIndex() const { return 3; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2);
    AddPermutation( 1,1,0,2);
  }
  std::string String(int shortversion) const 
  {
    return "1";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_SST());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_SST,LFSST_Getter,"SST","")
class LF_FFVT: public Lorentz_Function {
public:  
  LF_FFVT(): Lorentz_Function("FFVT") {}
  int NofIndex() const { return 2; }
  std::string String(int shortversion) const 
  {
    return "FFVT["+Str(0)+","+Str(1)+"]";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_FFVT());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_FFVT,LFFFVT_Getter,"FFVT","")
class LF_VVVT: public Lorentz_Function {
public: 
  LF_VVVT(): Lorentz_Function("VVVT") {}
  int NofIndex() const { return 4; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2,3);
    AddPermutation(-1,0,2,1,3);  
    AddPermutation(-1,1,0,2,3);
    AddPermutation(-1,2,1,0,3);  
    AddPermutation( 1,1,2,0,3);
    AddPermutation( 1,2,0,1,3);  
  }
  std::string String(int shortversion) const 
  {
    return "1";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_VVVT());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_VVVT,LFVVVT_Getter,"VVVT","")
class LF_SSST: public Lorentz_Function {
public:
  LF_SSST(): Lorentz_Function("SSST") {}
  int NofIndex() const { return 1; }
  std::string String(int shortversion) const 
  {
    return "1";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_SSST());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_SSST,LFSSST_Getter,"SSST","")
class LF_FFGS: public Lorentz_Function {
public:
  LF_FFGS(): Lorentz_Function("FFGS") {}
  int NofIndex() const { return 0; }
  std::string String(int shortversion) const 
  {
    return "1";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_FFGS());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_FFGS,LFFFGS_Getter,"FFGS","")
class LF_VVGS: public Lorentz_Function {
public:
  LF_VVGS(): Lorentz_Function("VVGS") {}
  int NofIndex() const { return 3; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2);
    AddPermutation( 1,1,0,2);
  }
  std::string String(int shortversion) const 
  {
    return "VVGS["+Str(0)+","+Str(1)+","+Str(2)+"]";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_VVGS());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_VVGS,LFVVGS_Getter,"VVGS","")
class LF_SSGS: public Lorentz_Function {
public: 
  LF_SSGS(): Lorentz_Function("SSGS") {}
  int NofIndex() const { return 2; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1);
    AddPermutation( 1,1,0);
  }
  std::string String(int shortversion) const 
  {
    return "1";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_SSGS());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_SSGS,LFSSGS_Getter,"SSGS","")
class LF_FFVGS: public Lorentz_Function {
public:
  LF_FFVGS(): Lorentz_Function("FFVGS") {}
  int NofIndex() const { return 1; }
  std::string String(int shortversion) const 
  {
    return "FFVGS["+Str(0)+"]";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_FFVGS());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_FFVGS,LFFFVGS_Getter,"FFVGS","")
class LF_Triangle: public Lorentz_Function {
public:
  LF_Triangle(): Lorentz_Function("Triangle") {}
  int NofIndex() const { return 2; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation(1,0,1);
    AddPermutation(1,1,0);
  }
  std::string String(int shortversion) const 
  {
    // G[0,1]
    return "T["+Str(0)+","+Str(1)+"]";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Triangle());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Triangle,LFTriangle_Getter,"Triangle","")
class LF_Box: public Lorentz_Function {
public:
  LF_Box(): Lorentz_Function("Box") {}
  int NofIndex() const { return 3; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2);
    AddPermutation(-1,0,2,1);  
    AddPermutation(-1,1,0,2);
    AddPermutation(-1,2,1,0);  
    AddPermutation( 1,1,2,0);
    AddPermutation( 1,2,0,1);  
  }
  std::string String(int shortversion) const 
  {
    // G[0,1]
    return "B["+Str(0)+","+Str(1)+","+Str(2)+"]";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Box());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Box,LFBox_Getter,"Box","")
class LF_C4GS: public Lorentz_Function {
public:
  LF_C4GS(): Lorentz_Function("C4GS") {}
  int NofIndex() const { return 2; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation(1,0,1);
    AddPermutation(1,1,0);  
  }
  std::string String(int shortversion) const 
  {
    // G[0,1]
    return "AddOn5Vertex["+Str(0)+","+Str(1)+"]";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_C4GS());
    *copy=*this;
    return copy;
  }
}; 
DEFINE_LF_GETTER(LF_C4GS,LFC4GS_Getter,"C4GS","")

class LF_PseudoTriangle: public Lorentz_Function {
public:
  LF_PseudoTriangle(): Lorentz_Function("PseudoTriangle") {}
  int NofIndex() const { return 2; }
  bool CutVectors() { return true; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation(1,0,1);
    AddPermutation(1,1,0);
  }
  std::string String(int shortversion) const 
  {
    // G[0,1]
    return "PsT["+Str(0)+","+Str(1)+"]";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_PseudoTriangle());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_PseudoTriangle,LFPseudoTriangle_Getter,"PseudoTriangle","")
class LF_PseudoBox: public Lorentz_Function {
public:
  LF_PseudoBox(): Lorentz_Function("PseudoBox") {}
  int NofIndex() const { return 3; }
  bool CutVectors() { return true; }
  void InitPermutation() 
  {
    Lorentz_Function::InitPermutation(); 
    AddPermutation( 1,0,1,2);
    AddPermutation(-1,0,2,1);  
    AddPermutation(-1,1,0,2);
    AddPermutation(-1,2,1,0);  
    AddPermutation( 1,1,2,0);
    AddPermutation( 1,2,0,1);  
  }
  std::string String(int shortversion) const 
  {
    // G[0,1]
    return "PsB["+Str(0)+","+Str(1)+","+Str(2)+"]";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_PseudoBox());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_PseudoBox,LFPseudoBox_Getter,"PseudoBox","")

