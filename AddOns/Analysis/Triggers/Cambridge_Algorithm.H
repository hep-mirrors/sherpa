#ifndef ATOOLS_Phys_Cambridge_Algorithm_H
#define ATOOLS_Phys_Cambridge_Algorithm_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include <vector>

using namespace ATOOLS;

namespace ANALYSIS {

  class Cambridge_Algorithm : public Jet_Algorithm_Base {
    int    m_mode, m_njets;
    double m_ycut, m_ymax;
    double m_sprime;

    int    m_matrixsize;
    double ** p_yij;
    int    *  p_imap;
    int    m_vectorsize;
    Vec4D  *  p_moms;
    bool   *  p_bflag;

    Particle_List       * p_jets;
    std::vector<double> * p_kts;

    inline double DCos12(const Vec4D &,const Vec4D &) const;
    inline double R12(const Vec4D &,const Vec4D &) const;

    void AddToKtlist(double );
    void AddToJetlist(const Vec4D &, bool);
  public:
    static double Kt2(const Vec4D & p);    

    Cambridge_Algorithm(ATOOLS::Particle_Qualifier_Base * const qualifier, int mode=0);
    ~Cambridge_Algorithm();

    void   Init(int);
    void   InitMoms(int);
    bool   ConstructJets(const Particle_List * ,Particle_List * ,std::vector<double> * ,double);

    void   Ymin(Vec4D *,bool *,int);
  };

  inline double Cambridge_Algorithm::Kt2(const Vec4D & p)
  {
    return sqr(p[1])+sqr(p[2]);
  }

  inline double Cambridge_Algorithm::DCos12(const Vec4D & p1,const Vec4D & p2) const
  {
    double s  = p1[1]*p2[1]+p1[2]*p2[2]+p1[3]*p2[3];
    double b1 = p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3];
    double b2 = p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3];
    return s/sqrt(b1*b2);
  }

  inline double Cambridge_Algorithm::R12(const Vec4D & p1, const Vec4D & p2) const
  {
    return 2.*(1.-DCos12(p1,p2))/m_sprime;
  }
}

#endif








