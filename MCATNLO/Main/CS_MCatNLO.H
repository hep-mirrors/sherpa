#ifndef MCATNLO_Main_CS_Shower_H
#define MCATNLO_Main_CS_Shower_H

#include "PDF/Main/NLOMC_Base.H"
#include "MCATNLO/Main/CS_Cluster_Definitions.H"
#include "MCATNLO/Showers/Shower.H"
#include "ATOOLS/Phys/Blob_List.H"

namespace MCATNLO {

  class Splitting_Function_Base;
  class CS_Gamma;

  typedef std::map<size_t,std::pair<double,double> > KT2X_Map;

  class CS_MCatNLO : public PDF::NLOMC_Base {
  
  private:
 
    PDF::ISR_Handler * p_isr;
    size_t m_maxem, m_nem, m_psmode;
    double m_maxweight;
    
    Shower          * p_mcatnlo;
    All_Singlets m_allsinglets;
    CS_Cluster_Definitions *p_cluster;
    CS_Gamma *p_gamma;
    All_Singlets *p_next;

    ATOOLS::Mass_Selector     *p_ms;
    ATOOLS::Cluster_Amplitude *p_rampl;

    Singlet *TranslateAmplitude(ATOOLS::Cluster_Amplitude *const ampl,
				std::map<ATOOLS::Cluster_Leg*,Parton*> &pmap,
				std::map<Parton*,ATOOLS::Cluster_Leg*> &lmap);

    int PerformMCatNLO(const size_t &maxem,size_t &nem,const double &qfac);

    bool PrepareMCatNLO(ATOOLS::Cluster_Amplitude *const ampl);

    void CleanUp();

  public:

    // constructor 
    CS_MCatNLO(PDF::ISR_Handler *const, MODEL::Model_Base *const,
	      ATOOLS::Default_Reader *const);

    // destructor
    ~CS_MCatNLO();

    //member functions
    int GeneratePoint(ATOOLS::Cluster_Amplitude *const ampl);

    void AddRBPoint(ATOOLS::Cluster_Amplitude *const ampl);

    ATOOLS::Cluster_Amplitude *GetRealEmissionAmplitude(const int mode=0);

    inline size_t PSMode() const { return m_psmode; }

  };

}

#endif
