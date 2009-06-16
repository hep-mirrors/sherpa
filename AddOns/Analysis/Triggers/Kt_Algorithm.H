#ifndef ATOOLS_Phys_Kt_Algorithm_H
#define ATOOLS_Phys_Kt_Algorithm_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "AddOns/Analysis/Tools/Particle_Qualifier.H"
#include <vector>

using namespace ATOOLS;

namespace ANALYSIS {
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
  
  class Order_ET {
  public:
    int operator()(ATOOLS::Particle * a, ATOOLS::Particle * b) {
      if (a->Momentum().EPerp() > b->Momentum().EPerp()) return 1;
      return 0;
    }
    int operator()(ATOOLS::Vec4D & a, ATOOLS::Vec4D & b) {
      if (a.EPerp() > b.EPerp()) return 1;
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
    ATOOLS::Particle_Qualifier_Base * p_qualifier;
    std::vector<double> m_jetnumbers;
    int   m_bflag;
    const ATOOLS::Blob_List *p_bl;
  public:
    Jet_Algorithm_Base(ATOOLS::Particle_Qualifier_Base * const qualifier) : 
      p_qualifier(qualifier), m_bflag(true) {}

    virtual bool ConstructJets(const Particle_List * ,Particle_List * ,std::vector<double> * ,double)=0;

    void Setbflag(int bf) { m_bflag=bf; }
    void SetBlobList(const ATOOLS::Blob_List *bl) { p_bl=bl; }
    void SortE(Particle_List *);
    void SortPT(Particle_List *);
    
    virtual ~Jet_Algorithm_Base();
    std::vector<double> & JetNumbers() { return m_jetnumbers; }
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
    void AddToJetlist(const Vec4D &, int);
  public:
    static double Kt2(const Vec4D & p);    

    Kt_Algorithm(ATOOLS::Particle_Qualifier_Base * const qualifier);
    ~Kt_Algorithm();

    void   Init(int);
    bool   ConstructJets(const Particle_List *,Particle_List * ,std::vector<double> * ,double);

    double Ktmin(Vec4D *,int *,int);
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








