#ifndef Flavour_H
#define Flavour_H

/*  Flavour.H: 
 *      Declarations for particle Flavour classes (for Flavour.C)
 *      The use of an explicit Init() function is ugly. Changes are needed.
 */

#include <fstream>
#include <string> 
#include <vector>
#include <set>
#include <iostream>

namespace ATOOLS {

  /////////////////////////////////////////////////////////////////////////////////////

  const long int MAX_PARTICLES=50000;

  void ParticleInit(std::string);
  void SetMassless();

  /////////////////////////////////////////////////////////////////////////////////////

  class kf {
  public:
    enum code {
      start                  = 0,

      /* ------------------------------------------
	Standard Model elementaries
      ------------------------------------------ */
      d                      = 1,
      u                      = 2,
      s                      = 3,
      c                      = 4,
      b                      = 5,
      t                      = 6,
      
      e                      = 11,
      nue                    = 12,
      mu                     = 13,
      numu                   = 14,
      tau                    = 15,
      nutau                  = 16,
      
      gluon                  = 21,
      photon                 = 22,
      W                      = 23,
      Z                      = 24,
      h                      = 25,
      
      /* ------------------------------------------
	LED - has to get other numbers to comply 
	       with the PDG numbering scheme
      ------------------------------------------ */
      graviton               = 27,
      gscalar                = 28,

      /* ------------------------------------------
	THDM - has to get other numbers to comply 
	       with the PDG numbering scheme
      ------------------------------------------ */
      h0                     = 31,
      H0                     = 32,
      A0                     = 33,
      Hmin                   = 34,
      
      /* ------------------------------------------
	SUSY - has to get other numbers to comply 
	       with the PDG numbering scheme
      ------------------------------------------ */
      Chargino1              = 41,
      Chargino2              = 42,
      Neutralino1            = 43,
      Neutralino2            = 44,
      Neutralino3            = 45,
      Neutralino4            = 46,
      Gluino                 = 47,
      
      sUpL                   = 51,
      sCharmL                = 52,
      sTop2                  = 53,
      sUpR                   = 54,
      sCharmR                = 55,
      sTop1                  = 56,
      
      sDownL                 = 61,
      sStrangeL              = 62,
      sBottom2               = 63,
      sDownR                 = 64,
      sStrangeR              = 65,
      sBottom1               = 66,
      
      sElectronL             = 71,
      sMuL                   = 72,
      sTau2                  = 73,
      sElectronR             = 74,
      sMuR                   = 75,
      sTau1                  = 76,
      
      sNu1                   = 81,
      sNu2                   = 82,
      sNu3                   = 83,
      
      /* ------------------------------------------
	helpers
      ------------------------------------------ */
      isrphoton              = 26,
      pol                    = 70,
      /* ------------------------------------------
	containers
      ------------------------------------------ */
      lepton                 = 90,
      neutrino               = 91,
      fermion                = 92,
      jet                    = 93,
      quark                  = 94,
      string                 = 95,
      cluster                = 96,
      bjet                   = 97,
      /* ------------------------------------------
	unspecified hadrons
      ------------------------------------------ */
      m_charged              = 100,
      m_neutral              = 101,
      b_charged              = 102,
      b_neutral              = 103,
      l_charged              = 104,
      l_neutral              = 105,
      
      /* ------------------------------------------
        diquarks
      ------------------------------------------ */
      dd_1                   = 1103,
      ud_0                   = 2101,
      ud_1                   = 2103,
      uu_1                   = 2203,
      sd_0                   = 3101,
      sd_1                   = 3103,
      su_0                   = 3201,
      su_1                   = 3203,
      ss_1                   = 3303,
      cd_0                   = 4101,
      cd_1                   = 4103,
      cu_0                   = 4201,
      cu_1                   = 4203,
      cs_0                   = 4301,
      cs_1                   = 4303,
      cc_1                   = 4403,
      bd_0                   = 5101,
      bd_1                   = 5103,
      bu_0                   = 5201,
      bu_1                   = 5203,
      bs_0                   = 5301,
      bs_1                   = 5303,
      bc_0                   = 5401,
      bc_1                   = 5403,
      bb_1                   = 5503,   

      /* ------------------------------------------
	specific hadrons
      ------------------------------------------ */
      pi                     = 111,
      rho_770                = 113,
      a_2_1320               = 115,
      rho_3_1690             = 117,
      K_L                    = 130,
      pi_plus                = 211,
      rho_770_plus           = 213,
      a_2_1320_plus          = 215,
      rho_3_1690_plus        = 217,
      eta                    = 221,
      omega_782              = 223,
      f_2_1270               = 225,
      omega_3_1670           = 227,
      f_4_2050               = 229,
      K_S                    = 310,
      K                      = 311,
      K_star_892             = 313,
      K_2_star_1430          = 315,
      K_3_star_1780          = 317,
      K_4_star_2045          = 319,
      K_plus                 = 321,
      K_star_892_plus        = 323,
      K_2_star_1430_plus     = 325,
      K_3_star_1780_plus     = 327,
      K_4_star_2045_plus     = 329,
      eta_prime_958          = 331,
      phi_1020               = 333,
      f_2_prime_1525         = 335,
      phi_3_1850             = 337,
      D_plus                 = 411,
      D_star_2010_plus       = 413,
      D_2_star_2460_plus     = 415,
      D                      = 421,
      D_star_2007            = 423,
      D_2_star_2460          = 425,
      D_s_plus               = 431,
      D_s_star_plus          = 433,
      eta_c_1S               = 441,
      J_psi_1S               = 443,
      chi_c2_1P              = 445,
      B                      = 511,
      B_star                 = 513,
      B_plus                 = 521,
      B_star_plus            = 523,
      B_s                    = 531,
      Upsilon_1S             = 553,
      chi_b2_1P              = 555,
      Delta_1620_minus       = 1112,
      Delta_1232_minus       = 1114,
      Delta_1905_minus       = 1116,
      Delta_1950_minus       = 1118,
      Delta_1620             = 1212,
      N_1520                 = 1214,
      Delta_1905             = 1216,
      N_2190                 = 1218,
      n                      = 2112,
      Delta_1232             = 2114,
      N_1675                 = 2116,
      Delta_1950             = 2118,
      Delta_1620_plus        = 2122,
      N_1520_plus            = 2124,
      Delta_1905_plus        = 2126,
      N_2190_plus            = 2128,
      p_plus                 = 2212,
      Delta_1232_plus        = 2214,
      N_1675_plus            = 2216,
      Delta_1950_plus        = 2218,
      Delta_1620_plus_plus   = 2222,
      Delta_1232_plus_plus   = 2224,
      Delta_1905_plus_plus   = 2226,
      Delta_1950_plus_plus   = 2228,
      Sigma_minus            = 3112,
      Sigma_1385_minus       = 3114,
      Sigma_1775_minus       = 3116,
      Sigma_2030_minus       = 3118,
      Lambda                 = 3122,
      Lambda_1520            = 3124,
      Lambda_1820            = 3126,
      Lambda_2100            = 3128,
      Sigma                  = 3212,
      Sigma_1385             = 3214,
      Sigma_1775             = 3216,
      Sigma_2030             = 3218,
      Sigma_plus             = 3222,
      Sigma_1385_plus        = 3224,
      Sigma_1775_plus        = 3226,
      Sigma_2030_plus        = 3228,
      Xi_minus               = 3312,
      Xi_1530_minus          = 3314,
      Xi                     = 3322,
      Xi_1530                = 3324,
      Omega_minus            = 3334,
      Sigma_c_2455           = 4112,
      Lambda_c_plus          = 4122,
      Xi_c                   = 4132,
      Sigma_c_2455_plus      = 4212,
      Sigma_c_2455_plus_plus = 4222,
      Xi_c_plus              = 4232,
      Omega_c                = 4332,
      Lambda_b               = 5122,
      b_1_1235               = 10113,
      pi_2_1670              = 10115,
      b_1_1235_plus          = 10213,
      pi_2_1670_plus         = 10215,
      h_1_1170               = 10223,
      K_0_star_1430          = 10311,
      K_1_1270               = 10313,
      K_2_1770               = 10315,
      K_0_star_1430_plus     = 10321,
      K_1_1270_plus          = 10323,
      K_2_1770_plus          = 10325,
      D_1_2420               = 10423,
      D_s1_2536_plus         = 10433,
      chi_c0_1P              = 10441,
      chi_b0_1P              = 10551,
      Delta_1700_minus       = 11114,
      Delta_1930_minus       = 11116,
      Delta_1930             = 11216,
      N_1440                 = 12112,
      Delta_1700             = 12114,
      N_1680                 = 12116,
      Delta_1930_plus        = 12126,
      N_1440_plus            = 12212,
      Delta_1700_plus        = 12214,
      N_1680_plus            = 12216,
      Delta_1700_plus_plus   = 12224,
      Delta_1930_plus_plus   = 12226,
      Sigma_1660_minus       = 13112,
      Sigma_1670_minus       = 13114,
      Sigma_1915_minus       = 13116,
      Lambda_1405            = 13122,
      Lambda_1690            = 13124,
      Lambda_1830            = 13126,
      Sigma_1660             = 13212,
      Sigma_1670             = 13214,
      Sigma_1915             = 13216,
      Sigma_1660_plus        = 13222,
      Sigma_1670_plus        = 13224,
      Sigma_1915_plus        = 13226,
      Xi_1820_minus          = 13314,
      Xi_1820                = 13324,
      Lambda_c_2593_plus     = 14122,
      a_1_1260               = 20113,
      a_1_1260_plus          = 20213,
      f_1_1285               = 20223,
      K_1_1400               = 20313,
      K_2_1820               = 20315,
      K_1_1400_plus          = 20323,
      K_2_1820_plus          = 20325,
      f_1_1420               = 20333,
      chi_c1_1P              = 20443,
      chi_b1_1P              = 20553,
      Delta_1910_minus       = 21112,
      Delta_1920_minus       = 21114,
      Delta_1910             = 21212,
      N_1700                 = 21214,
      N_1535                 = 22112,
      Delta_1920             = 22114,
      Delta_1910_plus        = 22122,
      N_1700_plus            = 22124,
      N_1535_plus            = 22212,
      Delta_1920_plus        = 22214,
      Delta_1910_plus_plus   = 22222,
      Delta_1920_plus_plus   = 22224,
      Sigma_1750_minus       = 23112,
      Sigma_1940_minus       = 23114,
      Lambda_1600            = 23122,
      Lambda_1890            = 23124,
      Lambda_2110            = 23126,
      Sigma_1750             = 23212,
      Sigma_1940             = 23214,
      Sigma_1750_plus        = 23222,
      Sigma_1940_plus        = 23224,
      rho_1700               = 30113,
      rho_1700_plus          = 30213,
      omega_1600             = 30223,
      K_star_1680            = 30313,
      K_star_1680_plus       = 30323,
      psi_3770               = 30443,
      Delta_1600_minus       = 31114,
      N_1720                 = 31214,
      N_1650                 = 32112,
      Delta_1600             = 32114,
      N_1720_plus            = 32124,
      N_1650_plus            = 32212,
      Delta_1600_plus        = 32214,
      Delta_1600_plus_plus   = 32224,
      Lambda_1670            = 33122,
      N_1710                 = 42112,
      N_1710_plus            = 42212,
      Lambda_1800            = 43122,
      Lambda_1810            = 53122,
      pi_1300                = 100111,
      rho_1450               = 100113,
      pi_1300_plus           = 100211,
      rho_1450_plus          = 100213,
      eta_1295               = 100221,
      omega_1420             = 100223,
      K_star_1410            = 100313,
      K_star_1410_plus       = 100323,
      phi_1680               = 100333,
      f_2_2010               = 100335,
      psi_2S                 = 100443,
      Upsilon_2S             = 100553,
      chi_b2_2P              = 100555,
      chi_b0_2P              = 110551,
      chi_b1_2P              = 120553,
      Upsilon_3S             = 200553,
      Upsilon_4S             = 300553,
      a_0_980                = 9000111,
      a_0_980_plus           = 9000211,
      psi_4040               = 9000443,
      Upsilon_10860          = 9000553,
      f_0_980                = 9010221,
      psi_4160               = 9010443,
      Upsilon_11020          = 9010553,
      f_0_1500               = 9020221,
      psi_4415               = 9020443,
      f_J_1710               = 9030221,
      f_2_2300               = 9050225,
      f_2_2340               = 9060225,

      none                   = 99999
    };
  };

  /////////////////////////////////////////////////////////////////////////////////////


  class Kf_To_Int {
    static int is_initialised;
    kf::code   kf_tab[MAX_PARTICLES];
    int        anz;
  public:
    void       Init();
    int        ToInt(kf::code);
    kf::code   FromInt(int);
    kf::code   FromString(std::string st);
    int IsInitialised() { return is_initialised;}
  };

  /////////////////////////////////////////////////////////////////////////////////////



  class Flavour; 

  struct Part_Info {
    kf::code kfc;
    double m,w,yuk;
    int    iq,isow,sp,masssign;
    bool   str,Maj,on,stbl,msv,hadron;
    char * n;

    std::vector<Flavour *> flavs;
    bool   group;

    Part_Info() { kfc=kf::none;m=0;iq=0;n=0;group=0; }
    Part_Info(kf::code kfc_, double mass, double width, int icharge, int isoweak, 
	      bool strong, int spin, bool Majorana, bool Take, bool stable,bool massive,
	      char* name,int _masssign);
    Part_Info(kf::code kfc_, double mass, double width, int icharge, int isoweak,
	      int spin, bool Take ,bool stable , char* name);
    void       Add(const Flavour);
    int  const Size();
    bool const Includes(Flavour);
    Flavour    operator[](const int);
    void       SetGroup();

    ~Part_Info() { delete[] n; };
    void operator = (const Part_Info& pi) {
      char* np = new char[strlen(pi.n)+1];
      strcpy(np,pi.n);
      delete[] n;
      n        = np;
      kfc      = pi.kfc;
      m        = pi.m;
      w        = pi.w;
      yuk      = pi.yuk;
      iq       = pi.iq;
      isow     = pi.isow;
      str      = pi.str;
      sp       = pi.sp;
      Maj      = pi.Maj;
      on       = pi.on;
      stbl     = pi.stbl;
      msv      = pi.msv;
      masssign = pi.masssign;
      flavs    = pi.flavs;
      group    = pi.group;
    }
  };

  extern Part_Info particles[MAX_PARTICLES];
  extern Kf_To_Int kf_table;

  /////////////////////////////////////////////////////////////////////////////////////


  class Flavour {
  protected:
    kf::code kfc;
    int    anti;
    inline int Index() const { return kf_table.ToInt(kfc); } 
    //    friend bool operator==(const Flavour&, const Flavour&);
  public:
    Flavour()                          { kfc=kf::none; anti=0;}
    Flavour(kf::code kfc_) : kfc(kfc_) { anti=0; }
    Flavour(int kfint);
    const Flavour Bar() const;
    kf::code Kfcode()     const { return kfc; }
    inline operator int() const { return (anti)? -int(kfc) : int(kfc) ; }

    int  Ctq();
    void FromCtq(int code);
    std::string TexName() const;

    const int  Size() const {
      return particles[ Index() ].Size();
    }
    const bool Includes(Flavour _fl) const {
      if ( (kfc==_fl.kfc) && (anti==_fl.anti))  return 1;
      if ( particles[ Index() ].Size()==1)      return 0;
      if ( particles[ Index() ].Includes(_fl) ) return 1;
      return 0;
    }

    Flavour operator[](int) const;

    Flavour& operator=(const Flavour& fl) {
      if (this!=&fl) {
	kfc   = fl.kfc;
	anti  = fl.anti;
      }
      return *this;
    }


    inline double Charge() const {
      double c=particles[ Index() ].iq/3.;
      return (anti) ? -c : c;
    }
    inline int IntCharge() const {
      int iq=particles[Index()].iq;
      return (anti) ? -iq : iq;
    }
    inline double IsoWeak() const {
      double c=particles[ Index() ].isow/2.;
      return (anti) ? -c : c;
    }
    inline bool Strong() const {
      return (particles[ Index() ].str !=0);
    }
    inline double Spin() const {
      return particles[ Index() ].sp/2.;
    }
    inline int IntSpin() const {
      return particles[ Index() ].sp;
    }
    inline bool Majorana() const {
      return particles[ Index() ].Maj;
    }
    inline bool IsOn() const {
      return particles[ Index() ].on;
    }
    inline bool IsStable() const {
      return particles[ Index() ].stbl;
    }
    inline void SetOn(bool on) const {
      particles[ Index() ].on = on;
    }
    inline void SetMassSign(int ms) const {
      particles[ Index() ].masssign = ms;
    }
    inline int MassSign() const {
      return particles[ Index() ].masssign;
    }
    inline bool IsMassive() const {
      if (particles[ Index() ].m > 0) return particles[ Index() ].msv;
      return 0;
    }
    inline double Mass() const {
      if (particles[ Index() ].msv) return particles[Index()].m;
      return 0.;
    }
    inline double PSMass() const {
      return particles[Index()].m;
    }
    inline void SetMass(double mass) const {
      particles[Index()].m = mass; 
    }
    inline double Yuk() const {
      return particles[Index()].yuk; 
    }
    inline void SetYuk(double _yuk) const {
      particles[Index()].yuk = _yuk; 
    }
    inline double Width() const {
      return particles[Index()].w; 
    }
    inline void SetWidth(double width) const {
      particles[Index()].w = width; 
    }
    inline char * Name() const {
      return particles[Index()].n; 
    }
    inline bool IsFermion() const {
      return (IntSpin()==1);
    }
    inline bool IsBoson() const {
      return ((IntSpin()==4) || (IntSpin()==2) || (IntSpin()==0)); 
    }
    inline bool IsScalar() const {
      return (IntSpin()==0); 
    }
    inline bool IsVector() const {
      return (IntSpin()==2); 
    }
    inline bool IsTensor() const {
      return (IntSpin()==4); 
    }
    inline int LeptonFamily() const {
      if (IsLepton()) return ((abs(kfc)-9)/2);
      else return 0;
    }
    inline int QuarkFamily() const {
      if ((IsQuark()) && (abs(kfc) < 10)) return ((abs(kfc)+1)/2);
      else return 0;
    }
    inline double BaryonNumber() {
      if (IsQuark() || IsSquark()) { return anti ? -1./3. : 1./3.; }
      return 0.;
    }
    inline int LeptonNumber() {
      if (IsLepton() || IsSlepton() || IsSneutrino()) { return anti ? -1 : 1; }
      return 0;
    }
    inline int StrongCharge() {
      if (IsQuark()) { return anti ? -3 : 3; }
      if (IsGluon()) { return -9; }
      return 0;
    }
    inline bool IsQuark() const {
      return ((Strong()>0) && (IntSpin()==1) && !(Majorana())); 
    }
    inline bool IsGluon() const {
      return (kfc==kf::gluon); 
    }
    inline bool IsJet() const {
      return (kfc==kf::jet); 
    }
    inline bool IsPhoton() const {
      return ((kfc==kf::photon) || (kfc==kf::isrphoton));
    }
    inline bool IsSquark() const  {
      return ((Strong()>0) && (IntSpin()==0) && !(Majorana())); 
    }
    inline bool IsGluino()  const {
      return (kfc==kf::Gluino); 
    }
    inline bool IsIno() {
      return (IsGluino() || IsNeutralino() || IsChargino());
    }                              
    inline bool IsKK() const {
      return (kfc==kf::graviton || kfc==kf::gscalar);
    }
    inline bool IsLepton() const {
      if ((abs(kfc)>10) && (abs(kfc)<20)) return 1;
      return 0;      
    }

    inline bool IsSlepton() const {
      if ((abs(kfc)>70) && (abs(kfc)<80)) return 1;
      return 0;      
    }
    inline bool IsSneutrino() const {
      if ((abs(kfc)>80) && (abs(kfc)<84)) return 1;
      return 0;      
    }
    inline bool IsChargino() const {
      if ((abs(kfc)>40) && (abs(kfc)<43)) return 1;
      return 0;      
    }
    inline bool IsNeutralino() const {
      if ((abs(kfc)>42) && (abs(kfc)<47)) return 1;
      return 0;      
    }
    inline bool IsUptype() const {
      switch (abs(kfc)) {
      case kf::u:     return 1;
      case kf::d:     return 0;
      case kf::s:     return 0;
      case kf::c:     return 1;
      case kf::b:     return 0;
      case kf::t:     return 1;
      case kf::e:     return 0;
      case kf::nue:   return 1;
      case kf::mu:    return 0;
      case kf::numu:  return 1;
      case kf::tau:   return 0;
      case kf::nutau: return 1;

      case kf::sDownL:    return 0;
      case kf::sDownR:    return 0;
      case kf::sStrangeL: return 0;
      case kf::sStrangeR: return 0;
      case kf::sBottom2:  return 0;
      case kf::sBottom1:  return 0; 
      case kf::sUpL:      return 1;
      case kf::sUpR:      return 1;
      case kf::sCharmL:   return 1;
      case kf::sCharmR:   return 1;
      case kf::sTop2:     return 1;
      case kf::sTop1:     return 1;
      case kf::sElectronL:return 0;
      case kf::sElectronR:return 0;
      case kf::sMuL:      return 0;
      case kf::sMuR:      return 0;
      case kf::sTau2:     return 0;
      case kf::sTau1:     return 0;
      case kf::sNu1:      return 1;
      case kf::sNu2:      return 1;
      case kf::sNu3:      return 1;
      default:break;
      }
      return 0;
    } 
    inline bool IsDowntype() const {
      switch (abs(kfc)) {
      case kf::u:     return 0;
      case kf::d:     return 1;
      case kf::s:     return 1;
      case kf::c:     return 0;
      case kf::b:     return 1;
      case kf::t:     return 0;
      case kf::e:     return 1;
      case kf::nue:   return 0;
      case kf::mu:    return 1;
      case kf::numu:  return 0;
      case kf::tau:   return 1;
      case kf::nutau: return 0;

      case kf::sDownL:    return 1;
      case kf::sDownR:    return 1;
      case kf::sStrangeL: return 1;
      case kf::sStrangeR: return 1;
      case kf::sBottom2:  return 1;
      case kf::sBottom1:  return 1; 
      case kf::sUpL:      return 0;
      case kf::sUpR:      return 0;
      case kf::sCharmL:   return 0;
      case kf::sCharmR:   return 0;
      case kf::sTop2:     return 0;
      case kf::sTop1:     return 0;
      case kf::sElectronL:return 1;
      case kf::sElectronR:return 1;
      case kf::sMuL:      return 1;
      case kf::sMuR:      return 1;
      case kf::sTau2:     return 1;
      case kf::sTau1:     return 1;
      case kf::sNu1:      return 0;
      case kf::sNu2:      return 0;
      case kf::sNu3:      return 0;
      default:break;
      }
      return 0;
    } 
    inline const bool IsLightQuark() const { return (abs(kfc)>=1 && abs(kfc)<=3);}
    inline const bool IsDiQuark() const {
      if(abs(kfc)>=1103 && abs(kfc)<=5505) {
	// better?:  if (abs(kfc)%100 <=3) return true;
	double help=abs(kfc)/100.0-int(abs(kfc)/100.0); if(help<0.031) return true;}
      return false;
    }
    inline const bool IsLightDiQuark() const {
      if(abs(kfc)>=1103 && abs(kfc)<=3303) {
	double help=abs(kfc)/100.0-int(abs(kfc)/100.0); if(help<0.031) return true;}
      return false;
    }
    inline const bool IsAnti() const {     
      return anti;
    }
    inline const bool IsHadron() const {     
      if (abs(kfc)<100) return 0;
      if (kfc-100*int(kfc/100)<10) return 0;
      return 1;
    }

    inline const bool IsSusy() const {     
      if ((30<kfc)&&(kfc<100)) return 1;
      return 0;
    }

    friend std::ostream& operator<<(std::ostream& os, const Flavour& f);
    int HepEvt();
    void FromHepEvt(int code);
    void FromHadron(int code);

    static int GetIntID(double value);
    static int Rot(int,int,int);
    int FlagID();
    int PropertiesID();
    static int ID_SM();
    static int ID_MSSM();
    static int ID_Had();
    static int WriteOut();
  };

  class Fl_Iter {
    int current;
  public:
    Fl_Iter() { current=0; }
    Flavour first() {
      current=kf_table.ToInt(kf::start)+1;
      return Flavour(current);
    }
    Flavour next() {
      if(Flavour(current)!=Flavour(kf::none))
	++current;
      return Flavour(current);
    }
  };

  struct Flavour_Pair {
    Flavour fl1;
    Flavour fl2;
  };

  class Flavour_Sorting_Criterion {
  public :
    bool operator() (const Flavour & _fl1,const Flavour & _fl2) {
      int fl1 = _fl1.Kfcode(), fl2 = _fl2.Kfcode();
      if (fl1<fl2) return 0;
      if (fl1==fl2 && _fl1.IsAnti() && !(_fl2.IsAnti())) return 0;
      return 1;
    }
  };

  typedef std::set<Flavour,Flavour_Sorting_Criterion> FlavourSet;
  typedef FlavourSet::iterator FlSetIter;
  typedef const FlavourSet::iterator FlSetConstIter;

  /*!
    \file 
    \brief contains Flavour
  */

  /*!
    \class Flavour
    \brief main class for accessing properties of particles
  */

}


#endif  // Flavour.H
