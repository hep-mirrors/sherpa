#ifndef Hadron_Remnant_H
#define Hadron_Remnant_H

#include "Remnant_Base.H"
#include "ISR_Handler.H"

namespace SHERPA {

  class Hadron_Remnant: public Remnant_Base {
  private:

    std::vector<ATOOLS::Flavour> m_constit;
    std::set<ATOOLS::Particle*>  m_companions, m_adjusted, m_singlet;

    PDF::PDF_Base *p_pdfbase;

    double       m_deltax, m_scale, m_ecms;
    unsigned int m_xscheme, m_maxtrials;

    ATOOLS::Vec4D m_hardpt;

    void GetConstituents(ATOOLS::Flavour flav);

    ATOOLS::Particle *FindConnected(ATOOLS::Particle *cur,bool=false,int=1); 
    ATOOLS::Particle *SelectCompanion(ATOOLS::Particle *cur);

    bool TestColours(ATOOLS::Particle *particle,int oldc,int newc,
		     bool singlet,bool force,int anti);
    bool AdjustColours(ATOOLS::Particle *particle,int oldc,int newc,
		       bool &singlet,bool force,int anti,bool foward);
    bool AdjustColours(ATOOLS::Particle *particle,int oldc,int newc,
		       bool &singlet,bool force);

    double GetXPDF(ATOOLS::Flavour flavour,double scale);
    void   DiceKinematics();

    bool TreatFirstQuark(ATOOLS::Particle *cur);
    bool TreatFirstGluon(ATOOLS::Particle *cur);
    bool TreatQuark(ATOOLS::Particle *cur);
    bool TreatGluon(ATOOLS::Particle *cur);

  public:
    
    // constructor
    Hadron_Remnant(PDF::ISR_Handler *isrhandler,const unsigned int _m_beam,double _m_scale);

    // member functions
    bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist);

    void UnDo(); 

  };// end of class Remnant_Base

}// end of namespace SHERPA

#endif
