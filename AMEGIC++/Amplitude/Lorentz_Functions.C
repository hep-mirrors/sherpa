#include "AMEGIC++/Amplitude/Lorentz_Function.H"

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
DEFINE_LF_GETTER(LF_None,"None","")

class LF_Gamma: public Lorentz_Function {
public:  
  LF_Gamma(): Lorentz_Function("FFV")
  { SetParticleArg(2,1,0); }
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
DEFINE_LF_GETTER(LF_Gamma,"FFV","")
class LF_Gab: public Lorentz_Function {
public:  
  LF_Gab(): Lorentz_Function("VVS")
  { SetParticleArg(0,1); }
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
DEFINE_LF_GETTER(LF_Gab,"VVS","")
class LF_Gauge3: public Lorentz_Function {
public:  
  LF_Gauge3(): Lorentz_Function("VVV")
  { SetParticleArg(0,1,2); }
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
    return "VVV["+Str(0)+","+Str(1)+","+Str(2)+"]";
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Gauge3());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Gauge3,"VVV","")
class LF_Gauge4: public Lorentz_Function {
public:  
  LF_Gauge4(): Lorentz_Function("VVVV")
  { SetParticleArg(0,1,2,3); }
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
    help  = std::string("(G[")  + Str(0) + std::string(",") + Str(2) + std::string("]*");
    help += std::string("G[")  + Str(1) + std::string(",") + Str(3) + std::string("]+");
    help += std::string("G[")  + Str(0) + std::string(",") + Str(3) + std::string("]*");
    help += std::string("G[")  + Str(1) + std::string(",") + Str(2) + std::string("]-");
    help += std::string("2*G[")  + Str(0) + std::string(",") + Str(1) + std::string("]*");
    help += std::string("G[")  + Str(2) + std::string(",") + Str(3) + std::string("])");
    return help;
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Gauge4());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Gauge4,"VVVV","")
class LF_Gluon4A: public Lorentz_Function {
public:  
  LF_Gluon4A(): Lorentz_Function("VVVVA")
  { SetParticleArg(0,1,2,3); }
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
      help += std::string("G4A[") + Str(0) + std::string(",") + 
	Str(1) + std::string(",") + 
	Str(2) + std::string(",") + 
	Str(3) + std::string("]");
    }
    else {
      help  = std::string("(G[")  + Str(0) + std::string(",") + Str(3) + std::string("]*");
      help += std::string("G[")  + Str(2) + std::string(",") + Str(1) + std::string("]-");
      help += std::string("G[")  + Str(0) + std::string(",") + Str(2) + std::string("]*");
      help += std::string("G[")  + Str(3) + std::string(",") + Str(1) + std::string("])");
    }
    return help;
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Gluon4A());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Gluon4A,"VVVVA","")
class LF_Gluon4B: public Lorentz_Function {
public:  
  LF_Gluon4B(): Lorentz_Function("VVVVB")
  { SetParticleArg(0,1,2,3); }
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
      help += std::string("G4B[") + Str(0) + std::string(",") + 
	Str(1) + std::string(",") + 
	Str(2) + std::string(",") + 
	Str(3) + std::string("]");
    }
    else {
      help  = std::string("(G[")  + Str(0) + std::string(",") + Str(3) + std::string("]*");
      help += std::string("G[")  + Str(1) + std::string(",") + Str(2) + std::string("]-");
      help += std::string("G[")  + Str(0) + std::string(",") + Str(1) + std::string("]*");
      help += std::string("G[")  + Str(3) + std::string(",") + Str(2) + std::string("])");
    }
    return help;
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Gluon4B());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Gluon4B,"VVVVB","")
class LF_Gluon4C: public Lorentz_Function {
public:  
  LF_Gluon4C(): Lorentz_Function("VVVVC")
  { SetParticleArg(0,1,2,3); }
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
      help += std::string("G4C[") + Str(0) + std::string(",") + 
	Str(1) + std::string(",") + 
	Str(2) + std::string(",") + 
	Str(3) + std::string("]");
    }
    else {
      help  = std::string("(G[")  + Str(0) + std::string(",") + Str(2) + std::string("]*");
      help += std::string("G[")  + Str(1) + std::string(",") + Str(3) + std::string("]-");
      help += std::string("G[")  + Str(0) + std::string(",") + Str(1) + std::string("]*");
      help += std::string("G[")  + Str(2) + std::string(",") + Str(3) + std::string("])");
    }
    return help;
  }
  Lorentz_Function *GetCopy() const 
  {
    Lorentz_Function *copy(new LF_Gluon4C());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_Gluon4C,"VVVVC","")
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
DEFINE_LF_GETTER(LF_SSV,"SSV","")

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
DEFINE_LF_GETTER(LF_SSS,"SSS","")
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
DEFINE_LF_GETTER(LF_FFS,"FFS","")
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
DEFINE_LF_GETTER(LF_Pol,"Pol","")
class LF_VVSS: public Lorentz_Function {
public:  
  LF_VVSS(): Lorentz_Function("VVSS")
  { SetParticleArg(0,1); }
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
DEFINE_LF_GETTER(LF_VVSS,"VVSS","")
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
DEFINE_LF_GETTER(LF_SSSS,"SSSS","")

class LF_HVV: public Lorentz_Function {
public:
  LF_HVV(): Lorentz_Function("HVV")
  { SetParticleArg(1,2); }
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
    Lorentz_Function *copy(new LF_HVV());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_HVV,"HVV","")
class LF_HVVV: public Lorentz_Function {
public:
  LF_HVVV(): Lorentz_Function("HVVV")
  { SetParticleArg(1,2,3); }
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
    Lorentz_Function *copy(new LF_HVVV());
    *copy=*this;
    return copy;
  }
};
DEFINE_LF_GETTER(LF_HVVV,"HVVV","")
class LF_C4GS: public Lorentz_Function {
public:
  LF_C4GS(): Lorentz_Function("C4GS")
  { SetParticleArg(1,2); }
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
DEFINE_LF_GETTER(LF_C4GS,"C4GS","")
