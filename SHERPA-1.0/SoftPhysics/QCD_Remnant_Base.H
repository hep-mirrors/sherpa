#ifndef QCD_Remnant_Base_H
#define QCD_Remnant_Base_H

#include "Remnant_Base.H"
#include "ISR_Handler.H"

namespace SHERPA {

  class QCD_Remnant_Base: public Remnant_Base {
  protected:

    std::vector<ATOOLS::Flavour> m_constit;
    std::set<ATOOLS::Particle*>  m_companions, m_adjusted, m_singlet;

    std::vector<std::pair<ATOOLS::Particle*,std::pair<int,int> > > m_undo;

    PDF::PDF_Base *p_pdfbase;

    double       m_deltax, m_scale, m_ecms;
    unsigned int m_xscheme, m_maxtrials;

    ATOOLS::Particle *FindConnected(ATOOLS::Particle *cur,bool=false,int=1); 
    ATOOLS::Particle *SelectCompanion(ATOOLS::Particle *cur);

    bool TestColours(ATOOLS::Particle *particle,int oldc,int newc,
		     bool singlet,bool force,int anti);
    bool AdjustColours(ATOOLS::Particle *particle,int oldc,int newc,
		       bool &singlet,bool force,int anti,bool foward);
    bool AdjustColours(ATOOLS::Particle *particle,int oldc,int newc,
		       bool &singlet,bool force);

  public:
    
    // constructor
    QCD_Remnant_Base(PDF::ISR_Handler *isrhandler,const unsigned int _m_beam,
		     const double _m_scale, const Remnant_Base::TypeID type);

    // member functions
    void UnDo(); 

  };// end of class QCD_Remnant_Base

}// end of namespace SHERPA

#endif
