#ifndef ATOOLS_Phys_Durham_Algorithm_H
#define ATOOLS_Phys_Durham_Algorithm_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include <vector>

namespace ANALYSIS {

  class Durham_Algorithm : public Jet_Algorithm_Base {
    int    m_mode;
    double m_ycut;//, m_r2min;
    double m_sprime;

    int    m_matrixsize;
    double ** p_yij;
    int    *  p_imap;
    int    m_vectorsize;
    Vec4D  *  p_moms;
    int    *  p_bflag;

    Particle_List       * p_jets;
    std::vector<double> * p_kts;

    double Y12(const Vec4D &,const Vec4D &) const;
    double DCos12(const Vec4D &,const Vec4D &) const;

    void AddToKtlist(double );
    void AddToJetlist(const Vec4D &, int);
  public:
    static double Kt2(const Vec4D & p);    

    Durham_Algorithm(ATOOLS::Particle_Qualifier_Base * const qualifier);
    ~Durham_Algorithm();

    void   Init(int);
    void   InitMoms(int);
    bool   ConstructJets(const Particle_List * ,Particle_List * ,std::vector<double> * ,double);

    void   Ymin(Vec4D *,int *,int);
  };

  inline double Durham_Algorithm::Kt2(const Vec4D & p)
  {
    return sqr(p[1])+sqr(p[2]);
  }



}

#endif








