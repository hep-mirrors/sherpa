#ifndef Dipole_Jet_Finder_h
#define Dipole_Jet_Finder_h

#include "Particle_List.H"
#include "Selector.H"
#include "Poincare.H"

namespace ATOOLS {
  struct dipjet_type {
    enum code {
      ariadne  = 0,
      csshower = 1,
      unknown  = -1
    };
  };

  struct dipjet_mode {
    enum code {
      dec     = 0,
      ee      = 1,
      DIS     = 2,
      hh      = 4,
      unknown = -1
    };
  };

  struct Reconstructed_Dipole {
  public:
    int i,j,k;
  };

  

  class Dipole_Jet_Finder : public Selector_Base {
    double            m_ycut, m_stot, m_kt2cut;
    dipjet_type::code m_type;
    dipjet_mode::code m_mode;
    bool              m_coloursorted;

    int                  m_nin,m_nout;
    double               m_shat, m_actual_y;
    double               m_kt2minII, m_kt2minIF, m_kt2minFI, m_kt2minFF, m_kt2min;
    Reconstructed_Dipole m_recII, m_recIF, m_recFI, m_recFF, m_rec;

    std::vector<Flavour> m_flavs;
    std::vector<Vec4D>   m_moms;

    void   II_Winner();
    void   FI_Winner();
    void   IF_Winner();
    void   FF_Winner();
    double KT_min_1I();
    double KT_min_II(const int,const int,const int);
    double KT_min_IF(const int,const int,const int,int &);
    double KT_min_FI(const int,const int,const int,int &);
    double KT_min_FF(const int,const int,const int,int &);

 public:
    Dipole_Jet_Finder(const int,Flavour *,const double,
		      const dipjet_type::code=dipjet_type::ariadne,
		      const dipjet_mode::code=dipjet_mode::hh,
		      const bool colsorted=false);
    Dipole_Jet_Finder(const double,
		      const dipjet_type::code=dipjet_type::ariadne,
		      const dipjet_mode::code=dipjet_mode::hh,
		      const bool=false);
    ~Dipole_Jet_Finder();

    void   Init(const Vec4D *);
    bool   Trigger(const Vec4D *);
    void   BuildCuts(Cut_Data *);
    void   UpdateCuts(double,double,Cut_Data *);
    double ActualValue() const;
  };

  inline double Dipole_Jet_Finder::ActualValue() const { return m_actual_y; }

}

#endif
