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

class VVSS4_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Mfunc,
		   public Basic_Vfunc {
  public:
    VVSS4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS); 
    ~VVSS4_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 2; }
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

  class G4_Calc : public Zfunc_Calc, 
		  public Basic_Zfunc,
		  public Basic_Xfunc,
		  public Basic_Vfunc {
  public:
    G4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~G4_Calc() {}
    ATOOLS::Kabbala Do(); 
  };

  /*class VVVV_Calc : public Zfunc_Calc, 
		    public Basic_Zfunc,
		    public Basic_Xfunc,
		    public Basic_Mfunc,
		    public Basic_Vfunc {
    ATOOLS::Kabbala GGGG();
    ATOOLS::Kabbala Triple_M1(const int &,const int &,const int &,const int &,const int &);
    ATOOLS::Kabbala Triple_M2(const int &,const int &,const int &,const int &,const int &);
    ATOOLS::Kabbala Double_M1(const int &,const int &,const int &,const int &,const int &);
    ATOOLS::Kabbala Double_M2(const int &,const int &,const int &,const int &,const int &);
    ATOOLS::Kabbala Double_M3(const int &,const int &,const int &,const int &,const int &);
    ATOOLS::Kabbala Single_M (const int &,const int &,const int &,const int &,const int &);
  public:
    VVVV_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~VVVV_Calc() {}
    ATOOLS::Kabbala Do();
    };*/

  /*class V5_Calc : public Zfunc_Calc, 
		  public Basic_Zfunc,
		  public Basic_Xfunc,
		  public Basic_Mfunc,
		  public Basic_Vfunc {
    ATOOLS::Kabbala G5();
    ATOOLS::Kabbala ZXXX1(const int&,const int&,const int&,
		  const int&,const int&,const int&,
		  const int&);
    ATOOLS::Kabbala ZXXX2(const int&,const int&,const int&,
		  const int&,const int&,const int&,
		  const int&);
    ATOOLS::Kabbala ZXXX3(const int&,const int&,const int&,
		  const int&,const int&,const int&,
		  const int&);
  public:
    V5_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~V5_Calc() {}
    ATOOLS::Kabbala Do();
    };*/

  /*class V4V3_Calc : public Zfunc_Calc, 
		    public Basic_Zfunc,
		    public Basic_Xfunc,
		    public Basic_Vfunc,
		    public Basic_Mfunc {
    ATOOLS::Kabbala ZXX(const int&,const int&,const int&,const int&);
    ATOOLS::Kabbala G4G3();
  public:
    V4V3_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~V4V3_Calc() {}
    ATOOLS::Kabbala Do();
    };*/


  class FFT_Calc : public Zfunc_Calc, 
		   public Basic_Yfunc,
		   public Basic_Xfunc,
		   public Basic_Vfunc {
  public:
    FFT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS); 
    ~FFT_Calc() {}
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };

  class VVT_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Vfunc {
  public:
    VVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS); 
    ~VVT_Calc() {}
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };

  class SST_Calc : public Zfunc_Calc, 
		   public Basic_Vfunc {
  public:
    SST_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS); 
    ~SST_Calc() {}
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };

  class FFVT_Calc : public Zfunc_Calc, 
		    public Basic_Xfunc,
		    public Basic_Zfunc,
		    public Basic_Mfunc,
		    public Basic_Vfunc {
  public:
    FFVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~FFVT_Calc() {}
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };

  class VVVT_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Vfunc {
  public:
    VVVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~VVVT_Calc() {}
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };

  class SSST_Calc : public Zfunc_Calc, 
		   public Basic_Vfunc {
  public:
    SSST_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~SSST_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 3; }
#endif
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };

  class FFGS_Calc : public Zfunc_Calc, 
		   public Basic_Yfunc,
		   public Basic_Xfunc {
  public:
    FFGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~FFGS_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() {return 1;}
#endif
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };

  class VVGS_Calc : public Zfunc_Calc, 
		   public Basic_Zfunc,
		   public Basic_Xfunc,
		   public Basic_Vfunc {
  public:
    VVGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~VVGS_Calc() {}
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };

  class SSGS_Calc : public Zfunc_Calc, 
		   public Basic_Vfunc {
  public:
    SSGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~SSGS_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() {return 1;} 
#endif
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };
  
  class FFVGS_Calc : public Zfunc_Calc, 
		     public Basic_Zfunc,
		     public Basic_Xfunc,
		     public Basic_Mfunc,
		     public Basic_Vfunc {
  public:
    FFVGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS); 
    ~FFVGS_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() {return 1;}
#endif
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
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

  class Pentagon_Calc : public Zfunc_Calc, 
			public Basic_Zfunc,
			public Basic_Xfunc,
			public Basic_Vfunc {
  public:
    Pentagon_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~Pentagon_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 1; }
#endif
    ATOOLS::Kabbala Do();
  };

  class PseudoTriangle_Calc : public Zfunc_Calc, 
			      public Basic_Zfunc,
			      public Basic_Xfunc,
			      public Basic_Vfunc,
			      public Basic_Epsilonfunc {
  public:
    PseudoTriangle_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~PseudoTriangle_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 1; }
#endif
    ATOOLS::Kabbala Do();
  };

  class PseudoBox_Calc : public Zfunc_Calc, 
			 public Basic_Zfunc,
			 public Basic_Xfunc,
			 public Basic_Vfunc,
			 public Basic_Epsilonfunc {
  public:
    PseudoBox_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~PseudoBox_Calc() {}
#ifdef Scalar_Args
    int GetScalarNumb() { return 1; }
#endif
    ATOOLS::Kabbala Do();
  };

  class SSVgen_Calc : public Zfunc_Calc,
		      public Basic_Xfunc,
		      public Basic_Mfunc,
		      public Basic_Vfunc {
  public:
    SSVgen_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    ~SSVgen_Calc() {}
    ATOOLS::Kabbala Do();
    void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		 Point *const p,Point *const pf,Point *&pb,
		 int *lfnumb,int *canumb);
  };
}
#endif
