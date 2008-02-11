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
#include <map>
#include <iostream>

namespace ATOOLS {

  /////////////////////////////////////////////////////////////////////////////

  const long int MAX_PARTICLES=50000;

  void ParticleInit(std::string);
  void SetMassless();

  /////////////////////////////////////////////////////////////////////////////

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
      b_prime                = 7,
      t_prime                = 8,

      e                      = 11,
      nue                    = 12,
      mu                     = 13,
      numu                   = 14,
      tau                    = 15,
      nutau                  = 16,
      tau_prime              = 17,
      nutau_prime            = 18,
      
      gluon                  = 21,
      photon                 = 22,
      Z                      = 23,
      Wplus                  = 24,
      h0                     = 25,
      
      Z0_2                   = 32,
      Z0_3                   = 33,
      Wplus_2                = 34,
      H0                     = 35,
      A0                     = 36,
      Hplus                  = 37,
      Hplusplus              = 38,

      //Graviton               = 39,
      //GScalar                = 40,

      /* ------------------------------------------
	helpers
      ------------------------------------------ */
      isrphoton              = 26,
      pol                    = 70,
      gluonqgc               = 99,

      // dummy for shrinked gluon in 
      // 4gluon-higgs vertex

      shgluon                = 29,

      /* ------------------------------------------
	LED - has to get other numbers to comply 
	       with the PDG numbering scheme
      ------------------------------------------ */
      graviton               = 27,
      gscalar                = 28,

      /* ------------------------------------------
	SUSY - has to get other numbers to comply 
	       with the PDG numbering scheme
      ------------------------------------------ */
      Chargino1              = 1000024,
      Chargino2              = 1000037,
      Neutralino1            = 1000022,
      Neutralino2            = 1000023,
      Neutralino3            = 1000025,
      Neutralino4            = 1000035,
      Gluino                 = 1000021,
      
      sUpL                   = 1000002,
      sCharmL                = 1000004,
      sTop1                  = 1000006,
      sUpR                   = 2000002,
      sCharmR                = 2000004,
      sTop2                  = 2000006,
      
      sDownL                 = 1000001,
      sStrangeL              = 1000003,
      sBottom1               = 1000005,
      sDownR                 = 2000001,
      sStrangeR              = 2000003,
      sBottom2               = 2000005,

      
      sElectronL             = 1000011,
      sMuL                   = 1000013,
      sTau1                  = 1000015,
      sElectronR             = 2000011,
      sMuR                   = 2000013,
      sTau2                  = 2000015,
      
      sNu1                   = 1000012,
      sNu2                   = 1000014,
      sNu3                   = 1000016,
      
      /* ------------------------------------------
	KK modes: numbers 5k000xx, 6k000xx for k-th resonance 
                  (left-, right-handed)
      ------------------------------------------ */

      KK_dL_1                = 5100001,
      KK_uL_1                = 5100002,
      KK_sL_1                = 5100003,
      KK_cL_1                = 5100004,
      KK_b2_1                = 5100005,
      KK_t2_1                = 5100006,
      
      KK_eL_1                = 5100011,
      KK_nueL_1              = 5100012,
      KK_muL_1               = 5100013,
      KK_numuL_1             = 5100014,
      KK_tau2_1              = 5100015,
      KK_nutauL_1            = 5100016,
      
      KK_gluon_1             = 5100021,
      KK_B1_1                = 5100022,
      KK_W_1                 = 5100023,
      KK_W3_1                = 5100024,

      KK_H_1                 = 5100035,
      KK_A_1                 = 5100036,
      KK_Hplus_1             = 5100037,

      KK_dR_1                = 6100001,
      KK_uR_1                = 6100002,
      KK_sR_1                = 6100003,
      KK_cR_1                = 6100004,
      KK_b1_1                = 6100005,
      KK_t1_1                = 6100006,
      
      KK_eR_1                = 6100011,
      KK_muR_1               = 6100013,
      KK_tau1_1              = 6100015,

      KK_dL_2                = 5200001,
      KK_uL_2                = 5200002,
      KK_sL_2                = 5200003,
      KK_cL_2                = 5200004,
      KK_b2_2                = 5200005,
      KK_t2_2                = 5200006,
      
      KK_eL_2                = 5200011,
      KK_nueL_2              = 5200012,
      KK_muL_2               = 5200013,
      KK_numuL_2             = 5200014,
      KK_tau2_2              = 5200015,
      KK_nutauL_2            = 5200016,
      
      KK_gluon_2             = 5200021,
      KK_B1_2                = 5200022,
      KK_W_2                 = 5200023,
      KK_W3_2                = 5200024,

      KK_H_2                 = 5200035,
      KK_A_2                 = 5200036,
      KK_Hplus_2             = 5200037,

      KK_dR_2                = 6200001,
      KK_uR_2                = 6200002,
      KK_sR_2                = 6200003,
      KK_cR_2                = 6200004,
      KK_b1_2                = 6200005,
      KK_t1_2                = 6200006,
      
      KK_eR_2                = 6200011,
      KK_muR_2               = 6200013,
      KK_tau1_2              = 6200015,
      


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
      seaquark               = 98,
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
	 ------------------------------------------
	 ------------------------------------------

 	 specific hadrons

	 ------------------------------------------
	 ------------------------------------------
         ------------------------------------------ */


      /* ------------------------------------------

         mesons

         ------------------------------------------ */
      pi                     = 111,
      K_L                    = 130,
      pi_plus                = 211,
      eta                    = 221,
      K_S                    = 310,
      K                      = 311,
      K_plus                 = 321,
      eta_prime_958          = 331,
      D_plus                 = 411,
      D                      = 421,
      D_s_plus               = 431,
      eta_c_1S               = 441,
      B                      = 511,
      B_plus                 = 521,
      B_s                    = 531,
      B_c                    = 541,
      eta_b                  = 551,

      rho_770                = 113,
      rho_770_plus           = 213,
      omega_782              = 223,
      K_star_892             = 313,
      K_star_892_plus        = 323,
      phi_1020               = 333,
      D_star_2010_plus       = 413,
      D_star_2007            = 423,
      D_s_star_plus          = 433,
      J_psi_1S               = 443,
      B_star                 = 513,
      B_star_plus            = 523,
      B_s_star               = 533,
      B_c_star               = 543,
      Upsilon_1S             = 553,

      a_2_1320               = 115,
      a_2_1320_plus          = 215,
      f_2_1270               = 225,
      K_2_star_1430          = 315,
      K_2_star_1430_plus     = 325,
      f_2_prime_1525         = 335,
      D_2_star_2460_plus     = 415,
      D_2_star_2460          = 425,
      D_s2_star_2573         = 435,
      chi_c2_1P              = 445,
      B_2_star               = 515,
      B_2_star_plus          = 525,
      B_s2_star              = 535,
      B_c2_star              = 545,
      chi_b2_1P              = 555,

      rho_3_1690             = 117,
      rho_3_1690_plus        = 217,
      omega_3_1670           = 227,
      K_3_star_1780          = 317,
      K_3_star_1780_plus     = 327,
      phi_3_1850             = 337,

      a_4_2040               = 119,
      a_4_2040_plus          = 219,
      f_4_2050               = 229,
      K_4_star_2045          = 319,
      K_4_star_2045_plus     = 329,
      f_J_2220               = 339,

      a_0_980                = 10111,
      a_0_980_plus           = 10211,
      f_0_980                = 10221,
      K_0_star_1430          = 10311,
      K_0_star_1430_plus     = 10321,
      f_0_1370               = 10331,
      D_0_star_plus          = 10411,
      D_0_star               = 10421,
      D_s0_star              = 10431,
      chi_c0_1P              = 10441,
      B_0_star               = 10511,
      B_0_star_plus          = 10521,
      B_s0_star              = 10531,
      B_c0_star              = 10541,
      chi_b0_1P              = 10551,

      b_1_1235               = 10113,
      b_1_1235_plus          = 10213,
      h_1_1170               = 10223,
      K_1_1270               = 10313,
      K_1_1270_plus          = 10323,
      h_1_1380               = 10333,
      D_1_2420_plus          = 10413,
      D_1_2420               = 10423,
      D_s1_2536_plus         = 10433,
      h_c1                   = 10443,
      B_1                    = 10513,
      B_1_plus               = 10523,
      B_s1                   = 10533,
      B_c1                   = 10543,
      h_b1                   = 10553,

      pi_2_1670              = 10115,
      pi_2_1670_plus         = 10215,
      eta_2_1645             = 10225,
      K_2_1770               = 10315,
      K_2_1770_plus          = 10325,
      eta_2_1870             = 10335,

      a_1_1260               = 20113,
      a_1_1260_plus          = 20213,
      f_1_1285               = 20223,
      K_1_1400               = 20313,
      K_1_1400_plus          = 20323,
      f_1_1420               = 20333,
      D_1_H_plus             = 20413,
      D_1_H                  = 20423,
      D_s1_H                 = 20433,
      B_1_H                  = 20513,
      B_1_H_plus             = 20523,
      B_s1_H                 = 20533,
      B_c1_H                 = 20533,
      chi_c1_1P              = 20543,
      chi_b1_1P              = 20553,

      K_2_1820               = 20315,
      K_2_1820_plus          = 20325,

      rho_1700               = 30113,
      rho_1700_plus          = 30213,
      omega_1600             = 30223,
      K_star_1680            = 30313,
      K_star_1680_plus       = 30323,
      f1_1900_fict           = 30333,
      psi_3770               = 30443,

      pi_1300                = 100111,
      pi_1300_plus           = 100211,
      eta_1295               = 100221,
      K_1460                 = 100311,
      K_1460_plus            = 100321,
      eta_1475               = 100331,

      rho_1450               = 100113,
      rho_1450_plus          = 100213,
      omega_1420             = 100223,
      K_star_1410            = 100313,
      K_star_1410_plus       = 100323,
      phi_1680               = 100333,
      psi_2S                 = 100443,
      Upsilon_2S             = 100553,

      f_2_2010               = 100335,
      chi_b2_2P              = 100555,

      chi_b0_2P              = 110551,
      chi_b1_2P              = 120553,
      Upsilon_3S             = 200553,
      Upsilon_4S             = 300553,

      a_0_1450               = 9000111,
      a_0_1450_plus          = 9000211,
      f_0_600                = 9000221,

      psi_4040               = 9000443,
      Upsilon_10860          = 9000553,

      f_0_1710               = 9010221,
      psi_4160               = 9010443,
      Upsilon_11020          = 9010553,
      f_0_1500               = 9020221,
      psi_4415               = 9020443,
      f_J_1710               = 9030221,
      f_2_2300               = 9050225,
      f_2_2340               = 9060225,

      /* ------------------------------------------

         baryons

         ------------------------------------------ */

      n                      = 2112,
      p_plus                 = 2212,
      Sigma_minus            = 3112,
      Lambda                 = 3122,
      Sigma                  = 3212,
      Sigma_plus             = 3222,
      Xi_minus               = 3312,
      Xi                     = 3322,
      Sigma_c_2455           = 4112,
      Lambda_c_plus          = 4122,
      Sigma_c_2455_plus      = 4212,
      Sigma_c_2455_plus_plus = 4222,
      Xi_c_2466              = 4132,
      Xi_c_2466_plus         = 4232,
      Xi_c_2574              = 4312,
      Xi_c_2574_plus         = 4322,
      Omega_c_0              = 4332,
      Sigma_b_5820_minus     = 5212,
      Lambda_b               = 5122,
      Sigma_b_5820           = 5212,
      Sigma_b_5820_plus      = 5222,
      Xi_b_5840              = 5132,
      Xi_b_5840_minus        = 5232,
      Xi_b_5960              = 5312,
      Xi_b_5960_minus        = 5322,
      Omega_b_0              = 5332,

      Delta_1232_minus       = 1114,
      Delta_1232             = 2114,
      Delta_1232_plus        = 2214,
      Delta_1232_plus_plus   = 2224,
      Sigma_1385_minus       = 3114,
      Sigma_1385             = 3214,
      Sigma_1385_plus        = 3224,
      Xi_1530_minus          = 3314,
      Xi_1530                = 3324,
      Omega_minus            = 3334,
      Sigma_c_2520           = 4114,
      Lambda_c_2625_plus     = 4124,
      Sigma_c_2520_plus      = 4214,
      Sigma_c_2520_plus_plus = 4224,
      Xi_c_2645              = 4314,
      Xi_c_2645_plus         = 4324,
      Omega_c_star           = 4334,
      Sigma_b_5840_minus     = 5114,
      Lambda_b_5940          = 5124,
      Sigma_b_5840           = 5214,
      Sigma_b_5840_plus      = 5224,
      Xi_b_5940              = 5314,
      Xi_b_5940_minus        = 5324,
      Omega_b_star           = 5334,

      N_1440                 = 12112,
      N_1440_plus            = 12212,
      Sigma_1660_minus       = 13112,
      Lambda_1600            = 13122,
      Sigma_1660             = 13212,
      Sigma_1660_plus        = 13222,
      Xi_1800_minus_fict     = 13312,
      Xi_1800_fict           = 13322,

      N_1535                 = 22112,
      N_1535_plus            = 22212,
      Sigma_1620_minus       = 23112,
      Lambda_1670            = 23122,
      Sigma_1620             = 23212,
      Sigma_1620_plus        = 23222,
      Xi_1720_minus_fict     = 23312,
      Xi_1720_fict           = 23322,
      Lambda_1405            = 22132,

      N_1520                 = 12114,
      N_1520_plus            = 12214,
      Sigma_1670_minus       = 13114,
      Lambda_1690            = 13124,
      Sigma_1670             = 13214,
      Sigma_1670_plus        = 13224,
      Xi_1820_minus          = 13314,
      Xi_1820                = 13324,
      Lambda_1520            = 12134,

      Delta_1700_minus       = 9011114,
      Delta_1700             = 9012114,
      Delta_1700_plus        = 9012214,
      Delta_1700_plus_plus   = 9012224,
      Sigma_1940_minus_fict  = 9013114,
      Sigma_1940_fict        = 9013214,
      Sigma_1940_plus_fict   = 9013224,
      Xi_2100_minus_fict     = 9013314,
      Xi_2100_fict           = 9013324,
      Omega_2300_minus_fict  = 9013334,

      none                   = 99999
    };
  };

  /////////////////////////////////////////////////////////////////////////////////////


  class Kf_To_Int {
    static int is_initialised;
    kf::code   kf_tab[MAX_PARTICLES];
    int        anz;
    std::map<kf::code,int> m_kftoint;
  public:
    void       Init();
    int        ToInt(kf::code);

    kf::code   FromString(std::string st);
    int IsInitialised() { return is_initialised;}
    int Size()          { return anz; }

    inline kf::code FromInt(int code) { return kf_tab[code]; }

  };

  /////////////////////////////////////////////////////////////////////////////



  class Flavour; 

  struct Part_Info {
    kf::code kfc;
    double m,w,yuk;
    int    iq,isow,sp,stbl,masssign;
    bool   str,Maj,on,msv,hadron;
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

    ~Part_Info();
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

  /////////////////////////////////////////////////////////////////////////////


  class Flavour {
  protected:
    kf::code kfc;
    int      anti;
    inline int Index() const { return kf_table.ToInt(kfc); } 
    //    friend bool operator==(const Flavour&, const Flavour&);
  public:
    Flavour()                          { kfc=kf::none; anti=0;}
    //NEW for ADICIC!
    Flavour(const Flavour& fin) : kfc(fin.kfc), anti(fin.anti) {}
    //***************
    Flavour(kf::code kfc_) : kfc(kfc_) { anti=0; }
    Flavour(kf::code, bool);    //NEW for ADICIC!
    inline Flavour(int kf): kfc(kf_table.FromInt(abs(kf))), anti(kf<0?1:0) {}
    const Flavour Bar() const;
    kf::code Kfcode()     const { return kfc; }
    inline operator int() const { return (anti)? -int(kfc) : int(kfc) ; }

    int  Ctq();
    void FromCtq(int code);
    std::string TexName() const;
    std::string IDName() const;
    std::string ShellName() const;

    const int  Size() const {
      return particles[ Index() ].Size();
    }
    const bool Includes(Flavour _fl) const {
      if ( (kfc==_fl.kfc) && (anti==_fl.anti))  return 1;
      int idx(Index());
      if ( particles[idx].Size()==1)      return 0;
      if ( particles[idx].Includes(_fl) ) return 1;
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
    inline void SetIntCharge(const int iq) 
    {
      particles[Index()].iq = iq;
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
      return particles[ Index() ].stbl&1;
    }
    inline int Stable() const {
      return particles[ Index() ].stbl;
    }
    inline void SetStable(const int stable) const {
      particles[ Index() ].stbl = stable;
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
      int idx(Index());
      return particles[idx].m ? particles[idx].msv : 0;
    }
    inline double Mass() const {
      int idx(Index());
      return particles[idx].msv ? particles[idx].m : 0.0;
    }
    inline double SelMass() const {
      int idx(Index());
      return particles[idx].msv && !IsKK() ? particles[idx].m : 0.0;
    }
    inline double PSMass() const {
      return particles[Index()].m;
    }
    inline void SetMassOn(bool on) const {
      particles[ Index() ].msv = on;
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
      return (kfc==kf::gluon || kfc==kf::shgluon); 
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
    inline bool IsLepton() const {
      if ((abs(kfc)>10) && (abs(kfc)<20)) return 1;
      return 0;      
    }

    inline bool IsSlepton() const {
      if (((abs(kfc)>1000010) && (abs(kfc)<1000016)) ||
	  ((abs(kfc)>2000010) && (abs(kfc)<2000016))) return 1;
      return 0;      
    }
    inline bool IsSneutrino() const {
      if ((abs(kfc)>1000011) && (abs(kfc)<1000017)) return 1;
      return 0;      
    }
    inline bool IsChargino() const {
      if ((abs(kfc)==1000024) || (abs(kfc)==1000037)) return 1;
      return 0;      
    }
    inline bool IsNeutralino() const {
      if ((abs(kfc)==1000022) || (abs(kfc)==1000023) ||
	  (abs(kfc)==1000025) || (abs(kfc)==1000035)) return 1;
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
      case kf::tau_prime:   return 0;
      case kf::nutau_prime: return 1;

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
      case kf::tau_prime:   return 1;
      case kf::nutau_prime: return 0;

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
      int exc(abs(kfc)/1000000);
      if (exc==1 || exc==2 || exc==3 || exc==4 || exc==5 || exc==6) return 0;
      return 1;
    }
    inline const bool IsBaryon() const {
      if (abs(kfc)<1000)             return 0;
      if (kfc-1000*int(kfc/1000)==0) return 0;
      return !IsDiQuark();
    }
    inline const bool IsB_Hadron() const {
      if (abs(kfc)<100)                            return 0;
      if (kfc-100*int(kfc/100)<10)                 return 0;
      if (abs((kfc-100*int(kfc/100))/10)==5)       return 1;
      if (abs((kfc-1000*int(kfc/1000))/100)==5)    return 1;
      if (abs((kfc-10000*int(kfc/10000))/1000)==5) return 1;
      return 0;
    }
    inline const bool IsC_Hadron() const {
      if (abs(kfc)<100)                            return 0;
      if (kfc-100*int(kfc/100)<10)                 return 0;
      if (abs((kfc-100*int(kfc/100))/10)==4)       return 1;
      if (abs((kfc-1000*int(kfc/1000))/100)==4)    return 1;
      if (abs((kfc-10000*int(kfc/10000))/1000)==4) return 1;
      return 0;
    }

    inline const bool IsSusy() const {     
      if ((35!=kfc)&&(30<kfc)&&(kfc<100)) return 1;
      return 0;
    }

    inline const bool IsKK() const {  
      if (kfc==kf::graviton || kfc==kf::gscalar) return 1;
      if (5000000<kfc && kfc<7000000) return 1;
      return 0;
    }

    inline const int KKGeneration() const {
      if (!IsKK()) return 0;
      return ((kfc-1000000*(kfc/1000000))/100000);
    }

     inline bool Is5VDummy() const {
      return (kfc==kf::shgluon); 
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

  class Flavour_Sorting_Criterion {
  public :
    bool operator() (const Flavour & _fl1,const Flavour & _fl2) {
      int fl1 = _fl1.Kfcode(), fl2 = _fl2.Kfcode();
      if (fl1<fl2) return 0;
      if (fl1==fl2 && _fl1.IsAnti() && !(_fl2.IsAnti())) return false;
      return true;
    }
  };

  typedef std::set<Flavour,Flavour_Sorting_Criterion> FlavourSet;
  typedef FlavourSet::iterator FlSetIter;
  typedef FlavourSet::const_iterator FlSetConstIter;

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




