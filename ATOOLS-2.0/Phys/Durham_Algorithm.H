#ifndef Durham_Algorithm_H
#define Durham_Algorithm_H

#include "Vector.H"
#include "Particle_List.H"
#include <vector>

namespace ATOOLS {

  class Durham_Algorithm  {
    int    m_mode;
    double m_ycut;//, m_r2min;
    double m_sprime;

    int    m_matrixsize;
    //double ** p_ktij;
    int    *  p_imap;
    //double *  p_kis;
    bool   m_bflag;

    Particle_List       * p_jets;
    std::vector<double> * p_kts;

    double Y12(const Vec4D &,const Vec4D &) const;
    double DCos12(const Vec4D &,const Vec4D &) const;

    void AddToKtlist(double );
    void AddToJetlist(const Vec4D &, bool);

    void SortE();
    void SortPT();
  public:
    static double Kt2(const Vec4D & p);    

    Durham_Algorithm(double =1.0);
    ~Durham_Algorithm();

    void   Init(int);
    bool   ConstructJets(const Particle_List * ,Particle_List * ,std::vector<double> * ,double);

    void   Ymin(Vec4D *,bool *,int);
    void   Setbflag(bool bf) { m_bflag=bf; }
  };

  inline double Durham_Algorithm::Kt2(const Vec4D & p)
  {
    return sqr(p[1])+sqr(p[2]);
  }



}

#endif








