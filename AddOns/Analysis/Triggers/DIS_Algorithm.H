#ifndef ATOOLS_Phys_DIS_Algorithm_H
#define ATOOLS_Phys_DIS_Algorithm_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "AddOns/Analysis/Tools/Particle_Qualifier.H"
#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include <vector>

using namespace ATOOLS;

namespace ANALYSIS {
  
  class DIS_Algorithm : public Jet_Algorithm_Base  {
    int    m_mode;
    double m_ymin, m_r2min;

    int    m_matrixsize;
    double ** p_ktij;
    int    *  p_imap;
    double *  p_kis;

    ATOOLS::Vec4D m_p, m_q;

    Particle_List       * p_jets;
    std::vector<double> * p_kts;

    double R2(const Vec4D &p1) const;
    double R2(const Vec4D &p1, const Vec4D &p2) const;
      
    void AddToKtlist(double );
    void AddToJetlist(const Vec4D &, int);
  public:
    static double Kt2(const Vec4D & p);    

    DIS_Algorithm(ATOOLS::Particle_Qualifier_Base * const qualifier);
    ~DIS_Algorithm();

    void   Init(int);
    bool   ConstructJets(const Particle_List *,Particle_List * ,std::vector<double> * ,double);

    double Ktmin(Vec4D *,int *,int);
  };

}

#endif
