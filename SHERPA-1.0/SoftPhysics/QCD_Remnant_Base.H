#ifndef QCD_Remnant_Base_H
#define QCD_Remnant_Base_H

#include "Remnant_Base.H"
#include "QCD_Remnant_Info.H"
#include "Parton_Finder.H"
#include "ISR_Handler.H"

namespace SHERPA {

  class QCD_Remnant_Base: public Remnant_Base {
  protected:

    std::vector<ATOOLS::Flavour> m_constit;
    std::set<ATOOLS::Particle*>  m_adjusted, m_singlet;

    std::vector<QCD_Remnant_Info*> m_sorted;

    PDF::PDF_Base *p_pdfbase;

    double       m_deltax, m_scale, m_ecms;
    unsigned int m_xscheme, m_maxtrials, m_set;
    size_t       m_initial, m_constrained;

    bool m_dupdf;

    ATOOLS::Parton_Finder m_finder;

    void SortRemnants(); 
    bool SelectCompanions(); 
    bool ConnectRemnants();

    bool TreatDipole(QCD_Remnant_Info *const cur,const size_t i);

    ATOOLS::Particle *FindConnected(ATOOLS::Particle *final,const bool anti) ;
    ATOOLS::Particle *FindClosest(ATOOLS::Particle *final,const bool anti); 
    ATOOLS::Particle *FindDisconnected(ATOOLS::Particle *final,const bool anti); 

    bool SelectCompanion(QCD_Remnant_Info *const cur);

    bool TestColours(ATOOLS::Particle *particle,unsigned int oldc,unsigned int newc,
		     bool singlet,bool force,int anti);
    bool AdjustColours(ATOOLS::Particle *particle,unsigned int oldc,unsigned int newc,
		       bool &singlet,bool force,int anti,bool foward);
    bool AdjustColours(ATOOLS::Particle *particle,unsigned int oldc,unsigned int newc,
		       bool &singlet,bool force);

  public:
    
    // constructor
    QCD_Remnant_Base(PDF::ISR_Handler *isrhandler,const unsigned int _m_beam,
		     const double _m_scale,const rtp::code type);

    // member functions
    void Clear();

    void UnDo(); 

  };// end of class QCD_Remnant_Base

}// end of namespace SHERPA

#endif
