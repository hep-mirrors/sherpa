#ifndef Amegic_Apacic_Interface_h
#define Amegic_Apacic_Interface_h

#include "Perturbative_Interface.H"

#include "Hard_Interface.H"
#include "XS_Selector.H" 
#include "XS_Group.H" 
#include "XS_Base.H" 

#include "Jet_Finder.H"
#include "Combine_Table.H"
#include "NLL_Sudakov.H"

#include "Tree.H"
#include "Knot.H"

#include "Process_Base.H"
#include "XS_Base.H"
#include "Jet_Finder.H"

#include "Blob.H"


namespace SHERPA {
  class Amegic_Apacic_Interface {
  private :
    bool                     m_NLLweight;
    double                   m_ycut,m_jetscale;
    int                      m_nin,m_nout;

    APHYTOOLS::Flavour     * p_flavs;
    Combine_Table          * p_ct;
    Combine_Table          * p_combi;
    NLL_Sudakov            * p_sud;

    std::string              m_lastname;
    APHYTOOLS::Jet_Finder  * p_jf;
    EXTRAXS::XS_Base       * p_xs; 
    EXTRAXS::XS_Group      * p_two2two;
    EXTRAXS::XS_Selector   * p_xsselector;
    APACIC::Hard_Interface * p_ps;

    bool            FillLegs(Leg *, AMEGIC::Point *, int &, int);
    void            EstablishRelations(APACIC::Knot *,APACIC::Knot *,APACIC::Knot *,bool);
    APACIC::Knot  * Point2Knot(APACIC::Tree *,const Leg &,const AMATOOLS::Vec4D &,char);
    AMATOOLS::Vec4D Momentum(APACIC::Knot *, int &);
    bool            IsColourConnected(APHYTOOLS::Parton *,APHYTOOLS::Parton *);
    int             colors[4][2];
    double          ColourAngle(const std::vector<APACIC::Knot *> &, const int);
    bool            ClusterConfiguration(AMEGIC::Process_Base *,APHYTOOLS::Blob *);
    void            CalculateWeight(double,double);
    void            FillTrees(APACIC::Tree **,APACIC::Tree *,EXTRAXS::XS_Base *);
    //    int  SetColours(AMATOOLS::Vec4D * , APHTOOLS::Flavour );

  public :
    Amegic_Apacic_Interface(Matrix_Element_Handler *,Shower_Handler *);
    ~Amegic_Apacic_Interface(); 
    bool   DefineInitialConditions(APHYTOOLS::Blob *);
    double GetWeight() { return m_weight; }
  };
}
#endif







