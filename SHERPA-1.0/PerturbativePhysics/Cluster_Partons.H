#ifndef ClusterPartons_h
#define ClusterPartons_h


#include "Combine_Table.H"
#include "NLL_Sudakov.H"

#include "Tree.H"
#include "Knot.H"

#include "Process_Base.H"
#include "XS_Base.H"
#include "Jet_Finder.H"

#include "Blob.H"

namespace SHERPA {
  class Matrix_Element_Handler;

  class Cluster_Partons {
    int                      maxjetnumber;
    double                   weight;
    Combine_Table          * ct;
    Combine_Table          * combi;

    NLL_Sudakov            * sud;

    //    AMEGIC::Process_Base   * lastproc;
    Matrix_Element_Handler * p_me;
    APHYTOOLS::Jet_Finder  * jf;
    APHYTOOLS::Blob        * blob;

    double scale;

    bool            FillLegs(Leg *, AMEGIC::Point *, int &, int);
    void            EstablishRelations(APACIC::Knot *,APACIC::Knot *,APACIC::Knot *,int);
    APACIC::Knot  * Point2Knot(APACIC::Tree *,const Leg &,const AMATOOLS::Vec4D &,char);
    AMATOOLS::Vec4D Momentum(APACIC::Knot *, int &);
    bool            IsColourConnected(APHYTOOLS::Parton *,APHYTOOLS::Parton *);
    int             colors[4][2];
    double          ColourAngle(const std::vector<APACIC::Knot *> &, const int);
    void            DetermineColourAngles(const std::vector<APACIC::Knot *> &);

  public :
    Cluster_Partons(Matrix_Element_Handler * me, APHYTOOLS::Jet_Finder * _jf,int _maxjetnumber) :
      p_me(me),jf(_jf),maxjetnumber(_maxjetnumber) {
      //      lastproc     = 0;
      combi        = 0;
      sud          = new NLL_Sudakov(jf->Smax(),jf->Smin());
    };

    // determine best cluster scheme, in order to initialise the PS
    //    bool ClusterConfiguration(AMEGIC::Process_Base *,APHYTOOLS::Blob *);
    bool ClusterConfiguration(APHYTOOLS::Blob *);
    // evaluate sudakov weight for matching ME to PS
    void CalculateWeight(double,double);
    // translate cluster sequence into initial or final state shower seeds 
    void FillTrees(APACIC::Tree **,APACIC::Tree *,EXTRAXS::XS_Base *);

    int  SetColours(AMATOOLS::Vec4D *, APHYTOOLS::Flavour *);

    double             Weight()      { return weight; }
    double             Scale() { return scale; }

    APHYTOOLS::Flavour Flav(int);
    AMATOOLS::Vec4D    Momentum(int); 

  };
}
#endif
