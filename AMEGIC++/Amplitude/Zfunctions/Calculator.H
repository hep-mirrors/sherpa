#ifndef AMEGIC_Amplitude_Zfunctions_Calculator_H
#define AMEGIC_Amplitude_Zfunctions_Calculator_H

#include "AMEGIC++/Amplitude/Zfunctions/Zfunc_Calc.H"

namespace AMEGIC {
  
  class Y_Calc : public Zfunc_Calc, public Basic_Yfunc {
  public:
    Y_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS); 
    ~Y_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 1; }
#endif
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };

  class Z_Calc : public Zfunc_Calc, 
		 public Basic_Zfunc,
                 public Basic_Xfunc,
                 public Basic_Mfunc,
		 public Basic_Vfunc {
  public:
    Z_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~Z_Calc() {}
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };

  class VVS_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Mfunc,
		   public Basic_Vfunc {
  public:
    VVS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~VVS_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 1; }
#endif
    ATOOLS::Kabbala Do();
  };

  class SSV_Calc : public Zfunc_Calc, 
		   public Basic_Xfunc,
		   public Basic_Mfunc,
		   public Basic_Vfunc {
  public:
    SSV_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~SSV_Calc() {}
    ATOOLS::Kabbala Do();
  };

  class SSS_Calc : public Zfunc_Calc {
  public:
    SSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~SSS_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 3; }
#endif
    ATOOLS::Kabbala Do();
  };
  
  class SSSS_Calc : public Zfunc_Calc {
  public:
    SSSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~SSSS_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 4; }
#endif
    ATOOLS::Kabbala Do();
  };

  class VVSS_Calc : public Zfunc_Calc, 
		    public Basic_Zfunc,
		    public Basic_Xfunc,
		    public Basic_Mfunc,
		    public Basic_Vfunc {
  public:
    VVSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~VVSS_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 2; }
#endif
    ATOOLS::Kabbala Do();
  };

  class VVV_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Mfunc,
		   public Basic_Vfunc {
    ATOOLS::Kabbala GGG();
  public:
    VVV_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~VVV_Calc() {}
    ATOOLS::Kabbala Do();
  };

  class V4_Calc : public Zfunc_Calc, 
		  public Basic_Zfunc,
		  public Basic_Xfunc,
		  public Basic_Mfunc,
		  public Basic_Vfunc {
    ATOOLS::Kabbala Massless();
  public:
    V4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~V4_Calc() {}
    ATOOLS::Kabbala Do();
  };

  class G4A_Calc : public Zfunc_Calc, 
		  public Basic_Zfunc,
		  public Basic_Xfunc,
		  public Basic_Vfunc {
  public:
    G4A_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~G4A_Calc() {}
    ATOOLS::Kabbala Do(); 
  };

  class G4B_Calc : public Zfunc_Calc, 
		  public Basic_Zfunc,
		  public Basic_Xfunc,
		  public Basic_Vfunc {
  public:
    G4B_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~G4B_Calc() {}
    ATOOLS::Kabbala Do(); 
  };

  class G4C_Calc : public Zfunc_Calc, 
		  public Basic_Zfunc,
		  public Basic_Xfunc,
		  public Basic_Vfunc {
  public:
    G4C_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~G4C_Calc() {}
    ATOOLS::Kabbala Do(); 
  };

  class Triangle_Calc : public Zfunc_Calc, 
			public Basic_Zfunc,
			public Basic_Xfunc,
			public Basic_Vfunc {
  public:
    Triangle_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~Triangle_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 1; }
#endif
    ATOOLS::Kabbala Do();
  };

  class Box_Calc : public Zfunc_Calc, 
			public Basic_Zfunc,
			public Basic_Xfunc,
			public Basic_Vfunc {
  public:
    Box_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~Box_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 1; }
#endif
    ATOOLS::Kabbala Do();
  };

  class PentagonA_Calc : public Zfunc_Calc, 
			public Basic_Zfunc,
			public Basic_Xfunc,
			public Basic_Vfunc {
  public:
    PentagonA_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~PentagonA_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 1; }
#endif
    ATOOLS::Kabbala Do();
  };

  class PentagonB_Calc : public Zfunc_Calc, 
			public Basic_Zfunc,
			public Basic_Xfunc,
			public Basic_Vfunc {
  public:
    PentagonB_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~PentagonB_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 1; }
#endif
    ATOOLS::Kabbala Do();
  };

  class PentagonC_Calc : public Zfunc_Calc, 
			public Basic_Zfunc,
			public Basic_Xfunc,
			public Basic_Vfunc {
  public:
    PentagonC_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~PentagonC_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 1; }
#endif
    ATOOLS::Kabbala Do();
  };

}
#endif
