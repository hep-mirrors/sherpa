#ifndef Zfunc_Calc_h
#define Zfunc_Calc_h

#include <strstream>
#include <vector>
#include "Basic_Func.H"

#define Scalar_Args

namespace AMEGIC {

  class zl {
  public:
    enum code {Y      = -1, 
	       Z      = 0, 
	       VVV    = 1,
	       VVS    = 2,
	       SSV    = 3,
	       SSS    = 4,
	       VVVV   = 5,
	       V4     = 6,
	       G4     = 7,
	       V5     = 9,
               V4V3   = 10, 
               VVSS   = 11,
	       VVSS4  = 12,//the generic interaction of 2 vectors and 2 scalars
	       SSSS   = 13,//the generic interaction of 4 scalars
	       ZNC    = 50,
	       VVVNC  = 51,
	       FFT    = 60,
	       VVT    = 61,
	       SST    = 62,
	       FFVT   = 63,
	       VVVT   = 64,
	       SSST   = 65,
	       FFGS   = 70,
	       VVGS   = 71,
	       SSGS   = 72,
	       FFVGS  = 73,
	       Unknown= 99}; 
  };
  
  class Zfunc_Calc : public virtual Basic_Func {
  public:
    int                      ncoupl,narg,pn;
    std::vector<Lorentz_Function> lorentzlist;
    zl::code                 type; 

    Zfunc_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {} 
    virtual ~Zfunc_Calc() {};
    virtual Kabbala Do() {std::cerr<<"Error: Virtual method Zfunc_Calc::Do() called!"<<std::endl;}
    virtual int GetScalarNumb() {return 0;}
  };
  
  class Y_Calc : public Zfunc_Calc, public Basic_Yfunc {
  public:
    Y_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) 
      : Basic_Func(_sgen,_BS), Zfunc_Calc(_sgen,_BS), Basic_Yfunc(_sgen,_BS) { 
      type     = zl::Y;
      ncoupl=2;narg=2;pn=1;
#ifdef Scalar_Args
      narg=3;
#endif
      lorentzlist.push_back(Lorentz_Function(lf::FFS));
    }
#ifdef Scalar_Args
    int GetScalarNumb() {return 1;}
#endif
    Kabbala Do();
  };

  class Z_Calc : public Zfunc_Calc, 
		 public Basic_Yfunc, 
		 public Basic_Zfunc,
                 public Basic_Xfunc,
                 public Basic_Vfunc,
                 public Basic_Mfunc,
                 public Basic_Pfunc {
  public:
    Z_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::Z;
      ncoupl=4;narg=4;pn=1;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist[0].SetParticleArg(0);
      lorentzlist[1].SetParticleArg(0);
    }
    Kabbala Do();
  };

  class VVS_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Mfunc,
		   public Basic_Vfunc,
		   public Basic_Pfunc {
  public:
    VVS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::VVS;
      ncoupl=5;narg=4;pn=2;
#ifdef Scalar_Args
      narg=5;
#endif
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gab));
      lorentzlist[0].SetParticleArg(0);
      lorentzlist[1].SetParticleArg(1);
      lorentzlist[2].SetParticleArg(0,1);
    }
#ifdef Scalar_Args
    int GetScalarNumb() {return 1;}
#endif
    Kabbala Do();
  };

class VVSS4_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Mfunc,
		   public Basic_Vfunc,
		   public Basic_Pfunc {
  public:
    VVSS4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::VVSS4;
      ncoupl=5;narg=4;pn=2;
#ifdef Scalar_Args
      narg=6;
#endif
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::VVSS));
      lorentzlist[0].SetParticleArg(0);
      lorentzlist[1].SetParticleArg(1);
      lorentzlist[2].SetParticleArg(0,1);
    }
#ifdef Scalar_Args
    int GetScalarNumb() {return 2;}
#endif
    Kabbala Do();
  };

  class SSV_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Mfunc,
		   public Basic_Vfunc,
		   public Basic_Pfunc {
  public:
    SSV_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::SSV;
      ncoupl=7;narg=6;pn=3;
#ifdef Scalar_Args
      //      narg=8;
#endif
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::SSV));
      lorentzlist[0].SetParticleArg(2);
      lorentzlist[1].SetParticleArg(0,1,2);
    }
#ifdef Scalar_Args
    //    int GetScalarNumb() {return 2;}
#endif
    Kabbala Do();

  };

  class SSS_Calc : public Zfunc_Calc {
  public:
    SSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS),
      Zfunc_Calc(_sgen,_BS)
    { 
      type     = zl::SSS;
      ncoupl=1;narg=0;pn=0;
#ifdef Scalar_Args
      narg=3;
#endif
      lorentzlist.push_back(Lorentz_Function(lf::SSS));
    }
#ifdef Scalar_Args
    int GetScalarNumb() {return 3;}
#endif
    Kabbala Do();
  };
  
  class SSSS_Calc : public Zfunc_Calc {
  public:
    SSSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS),
      Zfunc_Calc(_sgen,_BS)
    { 
      type     = zl::SSSS;
      ncoupl=1;narg=0;pn=0;
#ifdef Scalar_Args
      narg=4;
#endif
      lorentzlist.push_back(Lorentz_Function(lf::SSSS));
    }
#ifdef Scalar_Args
    int GetScalarNumb() {return 4;}
#endif
    Kabbala Do();
  };

  class VVSS_Calc : public Zfunc_Calc, 
		    public Basic_Zfunc,
		    public Basic_Xfunc,
		    public Basic_Mfunc,
		    public Basic_Vfunc,
		    public Basic_Pfunc {
  public:
    VVSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::VVSS;
      ncoupl=6;narg=4;pn=3;
#ifdef Scalar_Args
      narg=6;
#endif
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gab));
      lorentzlist.push_back(Lorentz_Function(lf::Gab));
      lorentzlist[0].SetParticleArg(0);
      lorentzlist[1].SetParticleArg(1);
      lorentzlist[2].SetParticleArg(0,2);
      lorentzlist[3].SetParticleArg(2,1);
    }
#ifdef Scalar_Args
    int GetScalarNumb() {return 2;}
#endif
    Kabbala Do();
  };


  class ZNC_Calc : public Zfunc_Calc, 
		   public Basic_Yfunc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Mfunc,
		   public Basic_VNCfunc,
		   public Basic_Pfunc {
  public:
    ZNC_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS), 
      Basic_VNCfunc(_sgen,_BS) 
    { 
      type     = zl::ZNC;
      ncoupl=4;narg=4;pn=1;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma_NC));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma_NC));
      lorentzlist[0].SetParticleArg(0);
      lorentzlist[1].SetParticleArg(0);
    }
    Kabbala Do();
  };

  class VVV_Calc : public Zfunc_Calc, 
		   public Basic_Yfunc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Mfunc,
		   public Basic_Vfunc,
		   public Basic_Pfunc {
    Kabbala GGG();
  public:
    VVV_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::VVV;
      ncoupl=10;narg=6;pn=3;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gauge3));
      for (short int i=0;i<3;i++) lorentzlist[i].SetParticleArg(i);
      lorentzlist[3].SetParticleArg(0,1,2);     
    }
    Kabbala Do();
  };

  class VVVNC_Calc : public Zfunc_Calc, 
		   public Basic_Yfunc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Mfunc,
		   public Basic_Vfunc,
		   public Basic_VNCfunc,
		   public Basic_Pfunc {
  public:
    VVVNC_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_VNCfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::VVVNC;
      ncoupl=10;narg=6;pn=3;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma_NC));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma_NC));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma_NC));
      lorentzlist.push_back(Lorentz_Function(lf::Photon3_NC));
      for (short int i=0;i<3;i++) lorentzlist[i].SetParticleArg(i);
      lorentzlist[3].SetParticleArg(0,1,2);     
    }
    Kabbala Do();
  };

  class V4_Calc : public Zfunc_Calc, 
		  public Basic_Yfunc, 
		  public Basic_Zfunc,
		  public Basic_Xfunc,
		  public Basic_Mfunc,
		  public Basic_Vfunc,
		  public Basic_Pfunc {
    Kabbala Massless();
  public:
    V4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::V4;
      ncoupl=9;narg=8;pn=4;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gauge4));
      for (short int i=0;i<4;i++) lorentzlist[i].SetParticleArg(i);
      lorentzlist[4].SetParticleArg(0,1,2,3);     
    }
    Kabbala Do();
  };

  class G4_Calc : public Zfunc_Calc, 
		  public Basic_Zfunc,
		  public Basic_Xfunc,
		  public Basic_Vfunc,
		  public Basic_Pfunc {
  public:
    G4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::G4;
      ncoupl=9;narg=8;pn=4;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gluon4));
      for (short int i=0;i<4;i++) lorentzlist[i].SetParticleArg(i);
      lorentzlist[4].SetParticleArg(0,1,2,3);     
    }
    Kabbala Do(); 
  };

  class VVVV_Calc : public Zfunc_Calc, 
		    public Basic_Yfunc, 
		    public Basic_Zfunc,
		    public Basic_Xfunc,
		    public Basic_Mfunc,
		    public Basic_Vfunc,
		    public Basic_Pfunc {
    Kabbala GGGG();
    Kabbala Triple_M1(const int &,const int &,const int &,const int &,const int &);
    Kabbala Triple_M2(const int &,const int &,const int &,const int &,const int &);
    Kabbala Double_M1(const int &,const int &,const int &,const int &,const int &);
    Kabbala Double_M2(const int &,const int &,const int &,const int &,const int &);
    Kabbala Double_M3(const int &,const int &,const int &,const int &,const int &);
    Kabbala Single_M (const int &,const int &,const int &,const int &,const int &);
  public:
    VVVV_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::VVVV;
      ncoupl=16;narg=8;pn=5;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gauge3));
      lorentzlist.push_back(Lorentz_Function(lf::Gauge3));
      for (short int i=0;i<4;i++) lorentzlist[i].SetParticleArg(i);
      lorentzlist[4].SetParticleArg(1,0,4);           
      lorentzlist[5].SetParticleArg(-4,2,3);     
    }
    Kabbala Do();
  };

  class V5_Calc : public Zfunc_Calc, 
		  public Basic_Yfunc, 
		  public Basic_Zfunc,
		  public Basic_Xfunc,
		  public Basic_Mfunc,
		  public Basic_Vfunc,
		  public Basic_Pfunc {
    Kabbala G5();
    Kabbala ZXXX1(const int&,const int&,const int&,
		  const int&,const int&,const int&,
		  const int&);
    Kabbala ZXXX2(const int&,const int&,const int&,
		  const int&,const int&,const int&,
		  const int&);
    Kabbala ZXXX3(const int&,const int&,const int&,
		  const int&,const int&,const int&,
		  const int&);
  public:
    V5_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::V5;
      ncoupl=19;narg=10;pn=7;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gauge3));
      lorentzlist.push_back(Lorentz_Function(lf::Gauge3));
      lorentzlist.push_back(Lorentz_Function(lf::Gauge3));
      for (short int i=0;i<5;i++) lorentzlist[i].SetParticleArg(i);
      lorentzlist[5].SetParticleArg(0,5,6);           
      lorentzlist[6].SetParticleArg(-5,1,2);     
      lorentzlist[7].SetParticleArg(-6,3,4);     
    }
    Kabbala Do();
  };

  class V4V3_Calc : public Zfunc_Calc, 
		    public Basic_Yfunc, 
		    public Basic_Zfunc,
		    public Basic_Xfunc,
		    public Basic_Mfunc,
		    public Basic_Vfunc,
		    public Basic_Pfunc {
    Kabbala ZXX(const int&,const int&,const int&,const int&);
    Kabbala G4G3();
  public:
    V4V3_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Pfunc(_sgen,_BS) 
    { 
      type     = zl::V4V3;
      ncoupl=18;narg=10;pn=6;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gluon4));
      lorentzlist.push_back(Lorentz_Function(lf::Gauge3));
      for (short int i=0;i<5;i++) lorentzlist[i].SetParticleArg(i);
      lorentzlist[5].SetParticleArg(0,1,5,2);           
      lorentzlist[6].SetParticleArg(-5,3,4);     
    }
    Kabbala Do();
  };
class FFT_Calc : public Zfunc_Calc, 
		   public Basic_Yfunc,
		   public Basic_Xfunc,
		   public Basic_Vfunc {
  public:
    FFT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS)
    { 
      type     = zl::FFT;
      ncoupl=4;narg=4;pn=3;
      lorentzlist.push_back(Lorentz_Function(lf::FFT));
      lorentzlist[0].SetParticleArg(0);
    }
    Kabbala Do();
  };

  class VVT_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Vfunc {
  public:
    VVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS) 
    { 
      type     = zl::VVT;
      ncoupl=7;narg=6;pn=3;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::VVT));
      lorentzlist[0].SetParticleArg(0);
      lorentzlist[1].SetParticleArg(1);
      lorentzlist[2].SetParticleArg(0,1,2);
    }
    Kabbala Do();
  };

  class SST_Calc : public Zfunc_Calc, 
		   public Basic_Vfunc {
  public:
    SST_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Vfunc(_sgen,_BS) 
    { 
      type     = zl::SST;
      ncoupl=7;narg=6;pn=3;
      lorentzlist.push_back(Lorentz_Function(lf::SST));
      lorentzlist[0].SetParticleArg(0,1,2);
    }
    Kabbala Do();
  };

  class FFVT_Calc : public Zfunc_Calc, 
		    public Basic_Yfunc,
		    public Basic_Xfunc,
		    public Basic_Zfunc,
		    public Basic_Mfunc,
		    public Basic_Vfunc {
  public:
    FFVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Zfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS)
      { 
	type     = zl::FFVT;
	ncoupl=7;narg=6;pn=2;
	lorentzlist.push_back(Lorentz_Function(lf::FFVT));
	lorentzlist.push_back(Lorentz_Function(lf::Gamma));
	lorentzlist[0].SetParticleArg(0,1);
	lorentzlist[1].SetParticleArg(0);
      }
    Kabbala Do();
  };

  class VVVT_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Vfunc {
  public:
    VVVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS) 
    { 
      type     = zl::VVVT;
      ncoupl=9;narg=8;pn=4;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::VVVT));
      for (short int i=0;i<3;i++) lorentzlist[i].SetParticleArg(i);
      lorentzlist[3].SetParticleArg(0,1,2,3);     
    }
    Kabbala Do();
  };

  class SSST_Calc : public Zfunc_Calc, 
		   public Basic_Vfunc {
  public:
    SSST_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Vfunc(_sgen,_BS) 
    { 
      type     = zl::SSST;
      ncoupl=2;narg=2;pn=1;
#ifdef Scalar_Args
      narg=5;
#endif
      lorentzlist.push_back(Lorentz_Function(lf::SSST));
      lorentzlist[0].SetParticleArg(0);     
    }
#ifdef Scalar_Args
    int GetScalarNumb() {return 3;}
#endif
    Kabbala Do();
  };

  class FFGS_Calc : public Zfunc_Calc, 
		   public Basic_Yfunc,
		   public Basic_Xfunc {
  public:
    FFGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS) 
    { 
      type     = zl::FFGS;
      ncoupl=3;narg=2;pn=3;
#ifdef Scalar_Args
      narg=3;
#endif
      lorentzlist.push_back(Lorentz_Function(lf::FFGS));
    }
#ifdef Scalar_Args
    int GetScalarNumb() {return 1;}
#endif
    Kabbala Do();
  };

  class VVGS_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Vfunc {
  public:
    VVGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS) 
    { 
      type     = zl::VVGS;
      ncoupl=8;narg=6;pn=3;
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::VVGS));
      lorentzlist[0].SetParticleArg(0);
      lorentzlist[1].SetParticleArg(1);
      lorentzlist[2].SetParticleArg(0,1,2);
    }
    Kabbala Do();
  };

  class SSGS_Calc : public Zfunc_Calc, 
		   public Basic_Vfunc {
  public:
    SSGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Vfunc(_sgen,_BS) 
    { 
      type     = zl::SSGS;
      ncoupl=6;narg=4;pn=2;
#ifdef Scalar_Args
      narg=5;
#endif
      lorentzlist.push_back(Lorentz_Function(lf::SSGS));
      lorentzlist[0].SetParticleArg(0,1);
    }
#ifdef Scalar_Args
    int GetScalarNumb() {return 1;} 
#endif
    Kabbala Do();
  };
  
  class FFVGS_Calc : public Zfunc_Calc, 
		     public Basic_Yfunc, 
		     public Basic_Zfunc,
		     public Basic_Xfunc,
		     public Basic_Vfunc,
		     public Basic_Mfunc {
  public:
    FFVGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
      Basic_Func(_sgen,_BS), 
      Zfunc_Calc(_sgen,_BS),
      Basic_Yfunc(_sgen,_BS), 
      Basic_Zfunc(_sgen,_BS), 
      Basic_Xfunc(_sgen,_BS), 
      Basic_Vfunc(_sgen,_BS), 
      Basic_Mfunc(_sgen,_BS)      { 
	type     = zl::FFVGS;
	ncoupl=4;narg=4;pn=1;
#ifdef Scalar_Args
      narg=5;
#endif
      lorentzlist.push_back(Lorentz_Function(lf::Gamma));
      lorentzlist.push_back(Lorentz_Function(lf::FFVGS));
      lorentzlist[0].SetParticleArg(0);
      lorentzlist[1].SetParticleArg(0);
    }
#ifdef Scalar_Args
    int GetScalarNumb() {return 1;}
#endif
    Kabbala Do();
  };

}
#endif



























