#ifndef Tree_Filler_H
#define Tree_Filler_H

#include "Initial_State_Shower.H"
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
 
    ATOOLS::Poincare       m_cms_boost;
    Flavour_Map          * p_flmap;

    double m_iss_scale_fac, m_fss_scale_fac;

    APACIC::Knot * Point2Knot(ATOOLS::Blob *,APACIC::Tree *,const Leg &,const ATOOLS::Vec4D &,char);
    double         ColourAngle(const std::vector<APACIC::Knot *> & knots, const int i);
    bool           IsColourConnected(ATOOLS::Particle *,ATOOLS::Particle *);
    void           DetermineColourAngles(const std::vector<APACIC::Knot *> & knots);
    void           EstablishRelations(APACIC::Knot * mo,APACIC::Knot * d1,APACIC::Knot * d2,
				      int mode,double scale=0.);

    ATOOLS::Vec4D  Momentum(APACIC::Knot * mo, int & number);
  public:
    Tree_Filler(Cluster_Partons_Base *,int,int,int);
    ~Tree_Filler();
    void FillTrees(ATOOLS::Blob *,APACIC::Tree ** ini_trees,APACIC::Tree * fin_tree);
    void FillDecayTree(APACIC::Tree * fin_tree);
    inline int    CKKWOn() const        { return m_ckkwon;   }

    inline void   SetCKKWOn(const int &on) { m_ckkwon=on;   }
  };
}

#endif
