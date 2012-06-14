#ifndef POWHEG_Main_CS_Shower_H
#define POWHEG_Main_CS_Shower_H

#include "PDF/Main/NLOMC_Base.H"
#include "POWHEG/Main/CS_Cluster_Definitions.H"
#include "POWHEG/Showers/Shower.H"
#include "ATOOLS/Phys/Blob_List.H"

namespace POWHEG {

  class Splitting_Function_Base;
  class CS_Gamma;

  typedef std::map<size_t,std::pair<double,double> > KT2X_Map;

  class CS_POWHEG : public PDF::NLOMC_Base {
  
  private:
 
    PDF::ISR_Handler * p_isr;
    size_t m_maxem, m_nem;
    
    Shower          * p_powheg;
    All_Singlets m_allsinglets;
    CS_Cluster_Definitions *p_cluster;
    CS_Gamma *p_gamma;
    All_Singlets *p_next;

    ATOOLS::Mass_Selector     *p_ms;
    ATOOLS::Cluster_Amplitude *p_rampl;

    Singlet *TranslateAmplitude(ATOOLS::Cluster_Amplitude *const ampl,
				std::map<ATOOLS::Cluster_Leg*,Parton*> &pmap,
				std::map<Parton*,ATOOLS::Cluster_Leg*> &lmap);

    int PerformPOWHEG(const size_t &maxem,size_t &nem);

    bool PreparePOWHEG(ATOOLS::Cluster_Amplitude *const ampl);

    void CleanUp();

  public:

    // constructor 
    CS_POWHEG(PDF::ISR_Handler *const isr, MODEL::Model_Base *const model,
	      ATOOLS::Data_Reader *const dataread); 

    // destructor
    ~CS_POWHEG();

    //member functions
    int GeneratePoint(ATOOLS::Cluster_Amplitude *const ampl);

    void AddRBPoint(ATOOLS::Cluster_Amplitude *const ampl);

    PDF::SH_Pair SHSplit(const double &B,const double &Qij2,
			 const ATOOLS::RB_Data *rbd) const;

    ATOOLS::Cluster_Amplitude *GetRealEmissionAmplitude();

  };

}

#endif
