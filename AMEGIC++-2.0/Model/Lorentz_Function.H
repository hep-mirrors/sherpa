#ifndef Lorentz_Function_h
#define Lorentz_Function_h

#include <strstream>
//#include <sstream>
#include <vector>
#include "MyString.H"

namespace AMEGIC {

  class lf {
  public:
    enum code {None   = 0, 
	       Gamma  = 1, 
	       Gab    = 2, 
	       Gauge3 = 3, 
	       Gauge4 = 4, 
	       Gluon4 = 5, 
	       SSV    = 6, 
	       SSS    = 7, 
	       FFS    = 8, 
	       Pol    =10,
	       VVSS   =11,//the generic interaction of 2 vectors and 2 scalars (Gab)
	       SSSS   =12,//the generic interaction of 4 scalars
	       Gamma_NC   = 51,
	       Photon3_NC = 53,
	       Photon4_NC = 55,
	       FFT    =60,
	       VVT    =61,
	       SST    =62,
	       FFVT   =63,
	       VVVT   =64,
	       SSST   =65,
               FFGS   =70,
	       VVGS   =71,
	       SSGS   =72,
	       FFVGS  =73,
	       Unknown=99}; 
  };
  
  class Lorentz_Function {
    int permcount;
    std::vector<int*> permlist;
    std::vector<int>  signlist;

    std::string Str(int a) const
    {
      std::strstream sstream;
      sstream<<partarg[a];
      std::string help;
      sstream>>help;
      return help;
    } 
  public:
    Lorentz_Function() {
      type = lf::Unknown;
      Next = 0;
      for (short int i=0;i<4;i++) partarg[i] = -1;
    }

    Lorentz_Function(lf::code _type) :type(_type) {
      Next = 0;
      for (short int i=0;i<4;i++) partarg[i] = -1;
    }
    
    lf::code  type; 
    int  partarg[4];
    Lorentz_Function* Next;
    int NofIndex() {
      switch (type) {
      case lf::Gamma : return 1;
      case lf::Gab   : return 2;
      case lf::Gauge3: return 3;
      case lf::Gauge4: return 4;
      case lf::Gluon4: return 4;
      case lf::SSV   : return 3;	
      case lf::Pol   : return 1;	
      case lf::VVSS  : return 2;	
      case lf::Gamma_NC   : return 1;
      case lf::Photon3_NC : return 3;
      case lf::Photon4_NC : return 4;
      case lf::FFT   : return 1;
      case lf::VVT   : return 3;
      case lf::SST   : return 3;
      case lf::FFVT  : return 2;
      case lf::VVVT  : return 4;
      case lf::SSST  : return 1;
      case lf::VVGS  : return 3;
      case lf::FFVGS : return 1;
      case lf::SSGS  : return 2;
      default: return 0;
      }
      //default (includes FFS, SSS and SSSS )
      return 0;
    }

    void SetParticleArg(int a=-1, int b=-1,int c=-1,int d=-1) {
      partarg[0] = a;partarg[1] = b;
      partarg[2] = c;partarg[3] = d;
    }

    int  GetSign();
    void AddPermutation(int,int,int,int,int);
    void InitPermutation();
    int  ResetPermutation();
    int  NextPermutation();

    std::string String(int shortversion=0) const {
      if (type==lf::SSS)  return std::string("1");
      if (type==lf::FFS)  return std::string("1");
      if (type==lf::SSSS) return std::string("1");
      std::string help;
      switch (type) {
      case lf::Gamma:  
	// Gam[0]
	help = std::string("Gam[") + Str(0) + std::string("]");break;
      case lf::Pol:  
	// Eps[0]
	help = std::string("Eps[") + Str(0) + std::string("]");break;
      case lf::Gab:   
	// G[0,1]
	help = std::string("G[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
      case lf::VVSS:   
	// G(2V2S)[0,1]
	help = std::string("G(2V2S)[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
      case lf::Gauge3: 
	// (P[0,2]-P[1,2])*G(0,1)+(P[1,0]-P[2,0])*G(1,2)+(P[2,1]-P[0,1])*G(2,0)
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
	break;
      case lf::SSV:
	//P[0,2]-P[1,2]
	help = std::string("P[") + Str(0) + std::string(",") + Str(2) +std::string("]-"); 
	help += std::string("P[") + Str(1) + std::string(",") + Str(2) +std::string("]");
	break;
      case lf::Gauge4: 
	//(2G(0,1)*G(2,3)-G(0,2)*G(1,3)-G(0,3)*G(1,2))
	help  = std::string("(2*G[")  + Str(0) + std::string(",") + Str(1) + std::string("]*");
	help += std::string("G[")  + Str(2) + std::string(",") + Str(3) + std::string("]-");
	help += std::string("G[")  + Str(0) + std::string(",") + Str(2) + std::string("]*");
	help += std::string("G[")  + Str(1) + std::string(",") + Str(3) + std::string("]-");
	help += std::string("G[")  + Str(0) + std::string(",") + Str(3) + std::string("]*");
	help += std::string("G[")  + Str(1) + std::string(",") + Str(2) + std::string("])");
	break;
      case lf::Gluon4:
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
	break;
      case lf::Gamma_NC:  
	// Gam_NC[0]
	help  = std::string("Exp[") + Str(0) + std::string("]*");
	help += std::string("Gam[") + Str(0) + std::string("]");break;
      case lf::Photon3_NC: 
	// (P[0,2]-P[1,2])*G(0,1)+(P[1,0]-P[2,0])*G(1,2)+(P[2,1]-P[0,1])*G(2,0)
	if (shortversion) {
	  help += std::string("P3NC[") + Str(0) + std::string(",") + 
	                            Str(1) + std::string(",") + 
                                    Str(2) + std::string("]");
	}
	else {
	  help  = std::string("Sin[p") + Str(0) + std::string("^") + Str(1) +std::string("]*(");

	  help += std::string("(P[") + Str(0) + std::string(",") + Str(2) + std::string("]-");
	  help += std::string("P[")  + Str(1) + std::string(",") + Str(2) + std::string("])*");
	  help += std::string("G[")  + Str(0) + std::string(",") + Str(1) + std::string("]");
	  
	  help += std::string("+");
	  
	  help += std::string("(P[") + Str(1) + std::string(",") + Str(0) + std::string("]-");
	  help += std::string("P[")  + Str(2) + std::string(",") + Str(0) + std::string("])*");
	  help += std::string("G[")  + Str(1) + std::string(",") + Str(2) + std::string("]");
	  
	  help += std::string("+");
	  
	  help += std::string("(P[") + Str(2) + std::string(",") + Str(1) + std::string("]-");
	  help += std::string("P[")  + Str(0) + std::string(",") + Str(1) + std::string("])*");
	  help += std::string("G[")  + Str(2) + std::string(",") + Str(0) + std::string("])");
	}
	break;
      case lf::Photon4_NC:
	//G(0,1)*G(2,3)-G(0,3)*G(2,1)
	if (shortversion) {
	  help += std::string("P4NC[") + Str(0) + std::string(",") + 
	                            Str(1) + std::string(",") + 
	                            Str(2) + std::string(",") + 
                                    Str(3) + std::string("]");
	}
	else {
	  //still to be done
	  help  = std::string("(G[")  + Str(0) + std::string(",") + Str(1) + std::string("]*");
	  help += std::string("G[")  + Str(2) + std::string(",") + Str(3) + std::string("]-");
	  help += std::string("G[")  + Str(0) + std::string(",") + Str(3) + std::string("]*");
	  help += std::string("G[")  + Str(2) + std::string(",") + Str(1) + std::string("])");
	}
	break;
      case lf::FFT:
	help = std::string("FFT[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
      case lf::FFVT:
	help = std::string("FFVT[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
      case lf::FFVGS:
	help = std::string("FFVGS[") + Str(0) + std::string("]");break;
      case lf::VVT:
	help = std::string("VVT[") + Str(0) + std::string(",") + Str(1) 
	     + std::string(",") + Str(2) + std::string("]");break;    
      case lf::VVGS:
	help = std::string("VVGS[") + Str(0) + std::string(",") + Str(1) 
	     + std::string(",") + Str(2) + std::string("]");break;
      }     
      return help;
    }
  };

  // Format of a Lorentz_Function to transport via MPI:
  struct MPI_Lorentz_Function {
    
    int m_type;
    int    m_partarg[4];         

   };

  //for MPI usage
  
  void Lorentz_Function2MPI(const Lorentz_Function * , MPI_Lorentz_Function &);
  
  Lorentz_Function * MPI2Lorentz_Function(const MPI_Lorentz_Function &);
  
  std::ostream & operator<<(std::ostream &, const MPI_Lorentz_Function &);

}
#endif















