#ifndef SHERPA_Tools_Output_RootNtuple_H
#define SHERPA_Tools_Output_RootNtuple_H

#include "SHERPA/Tools/Output_Base.H"
#ifdef USING__ROOT
#include "Riostream.h"
#include "TTree.h"
#include "TFile.h"
#endif


namespace SHERPA {

  struct rntuple_evt2 {
    double weight,wgt0;
    double x1,x2,y1,y2;
    double fscale,rscale;
    size_t nparticle;
    size_t id;
    int kf1,kf2;
    int nuwgt;
    double uwgt[18]; 
  };

  class Output_RootNtuple : public Output_Base {
#ifdef USING__ROOT
    static const Int_t s_kMaxParticle = 100;
    UInt_t m_id;
    Int_t m_nparticle;
    Float_t p_px[s_kMaxParticle];
    Float_t p_py[s_kMaxParticle];
    Float_t p_pz[s_kMaxParticle];
    Float_t p_E[s_kMaxParticle];
    Int_t p_kf[s_kMaxParticle];

    Double_t m_wgt,m_wgt2,m_mewgt,m_mewgt2;
    Int_t m_nuwgt;
    Double_t p_uwgt[18];
    Double_t m_x1,m_x2,m_y1,m_y2,m_fscale,m_rscale;
    Int_t m_id1, m_id2;
    TFile* p_f;
    TTree* p_t3;
#endif
    size_t m_cnt2,m_cnt3,m_fcnt,m_evt,m_idcnt;
    size_t m_avsize,m_total;
    double m_sum,m_s2,m_s3,m_c1,m_c2;
    double m_sq,m_sq2,m_sq3,m_fsq;
    std::vector<rntuple_evt2> m_evtlist;
    std::vector<int> m_flavlist;    
    std::vector<ATOOLS::Vec4D> m_momlist;
    void StoreEvt();
  public:
    Output_RootNtuple(std::string basename,std::string ext,int precision);
    ~Output_RootNtuple();
    void ChangeFile(std::string nr);
    void Output(ATOOLS::Blob_List* blobs, double weight);
  };

}

#endif
