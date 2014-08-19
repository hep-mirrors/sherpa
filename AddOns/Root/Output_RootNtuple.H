#ifndef SHERPA_Tools_Output_RootNtuple_H
#define SHERPA_Tools_Output_RootNtuple_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Output_Base.H"
#include "ATOOLS/Org/Exception.H"
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
    double alphas;
    long int id;
    int nparticle;
    int kf1,kf2;
    int nuwgt;
    double uwgt[18]; 
    int oqcd;
    char type[2];
  };

  class Output_RootNtuple : public Output_Base,
    public ATOOLS::Terminator_Object {
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
    Short_t m_oqcd;
    Char_t m_type[2];
    Double_t m_x1,m_x2,m_y1,m_y2,m_fscale,m_rscale;
    Int_t m_id1, m_id2;
    TFile* p_f;
    TTree* p_t3;
#endif
    std::string m_basename, m_ext, m_treename;
    size_t m_cnt2,m_cnt3,m_fcnt,m_evt,m_idcnt;
    size_t m_filesize,m_avsize,m_total,m_mode;
    double m_sum,m_s2,m_s3,m_c1,m_c2;
    double m_sq,m_sq2,m_sq3,m_fsq;
    double m_alphas;
    std::vector<rntuple_evt2> m_evtlist;
    std::vector<int> m_flavlist;    
    std::vector<ATOOLS::Vec4D> m_momlist;
    void StoreEvt();
    void MPISync();
    void AddDecayProducts(ATOOLS::Particle *part,int &np);
  public:
    Output_RootNtuple(const Output_Arguments &args);
    ~Output_RootNtuple();
    void Header();
    void Footer();
    void ChangeFile();
    void Output(ATOOLS::Blob_List* blobs, double weight);
    void PrepareTerminate();
  };

}

#endif
