#ifndef Kt_Algorithm_H
#define Kt_Algorithm_H

#include "Vector.H"
#include "Particle_List.H"
#include "Particle_Qualifier.H"
#include <vector>

namespace ATOOLS {
  class Order_PT {
  public:
    int operator()(ATOOLS::Particle * a, ATOOLS::Particle * b) {
      if (a->Momentum().PPerp2() > b->Momentum().PPerp2()) return 1;
      return 0;
    }
    int operator()(ATOOLS::Vec4D & a, ATOOLS::Vec4D & b) {
      if (a.PPerp2() > b.PPerp2()) return 1;
      return 0;
    }
  };
  
  class Order_E {
  public:
    int operator()(ATOOLS::Particle * a, ATOOLS::Particle * b) {
      if (a->Momentum()[0] > b->Momentum()[0]) return 1;
      return 0;
    }
  };

  class Jet_Algorithm_Base {
  protected:
    bool   m_bflag;
    ATOOLS::Particle_Qualifier_Base * p_qualifier;
  public:
    Jet_Algorithm_Base(ATOOLS::Particle_Qualifier_Base * const qualifier) : 
      p_qualifier(qualifier) {}

    virtual bool ConstructJets(const Particle_List * ,Particle_List * ,std::vector<double> * ,double)=0;

    void Setbflag(bool bf) { m_bflag=bf; }
    void SortE(Particle_List *);
    void SortPT(Particle_List *);
    
    virtual ~Jet_Algorithm_Base();
  };

  class Kt_Algorithm : public Jet_Algorithm_Base  {
    int    m_mode;
    double m_ymin, m_r2min;

    int    m_matrixsize;
    double ** p_ktij;
    int    *  p_imap;
    double *  p_kis;

    Particle_List       * p_jets;
    std::vector<double> * p_kts;

    double DRap12(const Vec4D &,const Vec4D &) const;
    double DEta12(const Vec4D &,const Vec4D &) const;
    double DPhi12(const Vec4D &,const Vec4D &) const;
    double CosDPhi12(const Vec4D &,const Vec4D &) const;
    double DCos12(const Vec4D &,const Vec4D &) const;

    double R2(const Vec4D &p1, const Vec4D &p2) const;
      
    void AddToKtlist(double );
    void AddToJetlist(const Vec4D &, bool);
  public:
    static double Kt2(const Vec4D & p);    

    Kt_Algorithm(ATOOLS::Particle_Qualifier_Base * const qualifier);
    ~Kt_Algorithm();

    void   Init(int);
    bool   ConstructJets(const Particle_List *,Particle_List * ,std::vector<double> * ,double);

    double Ktmin(Vec4D *,bool *,int);
  };


  inline double Kt_Algorithm::Kt2(const Vec4D & p)
  {
    return sqr(p[1])+sqr(p[2]);
  }

  inline double Kt_Algorithm::R2(const Vec4D &p1, const Vec4D &p2) const
  {
    return (sqr(DEta12(p1,p2)) + sqr(DPhi12(p1,p2)));
  }


}

#endif








