#ifndef YFS_Main_YFS_H
#define YFS_Main_YFS_H



#include "MODEL/Main/Running_AlphaQED.H"


namespace ATOOLS {
  class Blob;
  class Data_Reader;
}

namespace YFS {
  class Semi_YFS {
  private:
    std::string                 m_name;
    bool                        m_photons_added;
    bool                        m_success;
    //Yfs(ATOOLS::Data_Reader* reader=NULL);
    
    public:
      void   Initial();
        double Gen_V();
      double Beta();
      double Jac(double v);
      double Emass(double v,double alp, double s);
      double nbar();   //m_RN = -\sum{\log{v_i}}_i^N
      int    Photon_Mult();
      void   initial();
      double Gen_K();
      double Gen_Y();
      double Gen_Gamma(double alp);
      double Gen_exp(double alp);
      double Photon_Ang();
      double Total();
      // bool   YfsMode();
      double Weight(double x, double alp, const int &semiyfs);
      double bornxs(double s,double theta);
      double sigma_crude(double s,double v);
      double max_xs(double s,double vmin, double vmax);
      double m_rho;
      double m_intial_S;
      double m_BeamEnergy1;
      double m_BeamEnergy2;
      double mass_f;
      double m_Jac;
      double m_alqed;
      double m_Emass;
      double m_g;
      double m_gp;
      double m_epsilon;
      double m_Kmax;  
      double m_Kmin;
      double m_alpha_A;
      double m_yi;
      double m_ki;
      double m_var;
      double s;
      double v;
      double vmin;
      double vmax;
      double xsmax;
      double theta;
      int    m_Photon_Number; 
      int    m_yfsmode;
      double s_ir_cutoff;

      Semi_YFS();
      ~Semi_YFS();
        
      


    bool AddRadiation(ATOOLS::Blob *);
    bool CheckStateBeforeTreatment(ATOOLS::Blob *);

    inline bool        Done_Successfully() { return m_success; }
    inline bool        Added_Anything()    { return m_photons_added; }
    inline std::string Name()             { return m_name; }

    static double               s_alpha;
    static bool                 s_use_running_parameters;
    static inline void Reset_Alpha_QED()
    { s_alpha=MODEL::aqed->AqedThomson(); }
    static inline void Set_Alpha_QED(const double& scale)
    { if (s_use_running_parameters) s_alpha=(*MODEL::aqed)(scale); }

    static int                  s_YFS_mode;
    static bool                 s_use_me;
    static int                  s_ir_cutoffframe;
    static double               s_accu;
    static int                  s_n_max;
    static int                  s_n_min;
    static double               s_dr_cut;
    static bool                 s_strict;
    static double               s_increase_max_weight;
    static double               s_reduce_max_en;
    static bool                 s_check_first;
    static int                  s_ff_recscheme;
    static int                  s_fi_recscheme;
    int                         s_YFS;
    int                         s_seed;
    double                      s_var;

 

  };

}

#endif


