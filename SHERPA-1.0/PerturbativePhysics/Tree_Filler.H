#ifndef Tree_Filler_H
#define Tree_Filler_H

#include "Shower_Handler.H"
#include "Final_State_Shower.H"
#include "Cluster_Partons_Base.H"
#include "Blob.H"
#include "Poincare.H"

namespace SHERPA {
  class Tree_Filler {
  private:
    Cluster_Partons_Base * p_cluster;
    int                    m_maxjetnumber, m_isrshoweron, m_fsrshoweron, m_ckkwon;
    APACIC::Tree         * p_local_tree;

    APACIC::Final_State_Shower * p_fsrshower;
 
    ATOOLS::Poincare       m_cms_boost;

    double m_iss_scale_fac, m_fss_scale_fac, m_q2_cut;

    int m_showermode;

    APACIC::Knot * Point2Knot(ATOOLS::Blob *,APACIC::Tree *,const Leg &,const ATOOLS::Vec4D &,char,ATOOLS::Particle *&mep);
    double         ColourAngle(const std::vector<APACIC::Knot *> & knots, const int i);
    bool           IsColourConnected(ATOOLS::Particle *,ATOOLS::Particle *);

    ATOOLS::Vec4D  Momentum(APACIC::Knot * mo, int & number);
  public:
    Tree_Filler(Cluster_Partons_Base *,Shower_Handler *shower,int,int);
    ~Tree_Filler();
    void FillTrees(ATOOLS::Blob *,APACIC::Tree ** ini_trees,APACIC::Tree * fin_tree);
    void FillDecayTree(APACIC::Tree * fin_tree);
    inline int    CKKWOn() const        { return m_ckkwon;   }

    void           DetermineColourAngles(const std::vector<APACIC::Knot *> & knots);
    void           EstablishRelations(APACIC::Knot * mo,APACIC::Knot * d1,APACIC::Knot * d2,
				      int mode,double x1,double x2,double scale1=0.,double scale2=0.);
    inline void   SetCKKWOn(const int &on) { m_ckkwon=on;   }
    inline void   SetCMSBoost(const ATOOLS::Vec4D &cms) 
    { m_cms_boost=ATOOLS::Poincare(cms); }
    inline ATOOLS::Poincare &CMSBoost() { return m_cms_boost; }
    inline void SetQ2Cut(const double &q2) { m_q2_cut=q2; }
    inline double Q2Cut() const { return m_q2_cut; }
  };
}

#endif
