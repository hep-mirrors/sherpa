#ifndef HADRONS_ME_Library_Three_Body_MEs_H
#define HADRONS_ME_Library_Three_Body_MEs_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HADRONS {

  class B_Bpi_pwave : public HD_ME_Base {
    Complex m_cL, m_cR;
  public:
    B_Bpi_pwave(const ATOOLS::Flavour_Vector& flavs,
                const std::vector<int>& decayindices,
                const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
  };

  class B3_Bpi_pwave : public HD_ME_Base {
    Complex m_cL, m_cR;
  public:
    B3_Bpi_pwave(const ATOOLS::Flavour_Vector& flavs,
                 const std::vector<int>& decayindices,
                 const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
  };

  class B_Bpi_swave : public HD_ME_Base {
    Complex m_cL, m_cR;
  public:
    B_Bpi_swave(const ATOOLS::Flavour_Vector& flavs,
                const std::vector<int>& decayindices,
                const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
  };

  class B3_Bpi_dwave : public HD_ME_Base {
    Complex m_cL, m_cR;
  public:
    B3_Bpi_dwave(const ATOOLS::Flavour_Vector& flavs,
                 const std::vector<int>& decayindices,
                 const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
  };

  class B_Bphoton_M1 : public HD_ME_Base {
    double m_fac;
  public:
    B_Bphoton_M1(const ATOOLS::Flavour_Vector& flavs,
                 const std::vector<int>& decayindices,
                 const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
  };

  class B3_Bphoton_M1 : public HD_ME_Base {
    Complex m_cL, m_cR;
  public:
    B3_Bphoton_M1(const ATOOLS::Flavour_Vector& flavs,
                  const std::vector<int>& decayindices,
                  const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
  };

  class B3_Bphoton_M1_2 : public HD_ME_Base {
    Complex m_cL, m_cR;
  public:
    B3_Bphoton_M1_2(const ATOOLS::Flavour_Vector& flavs,
                    const std::vector<int>& decayindices,
                    const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
  };

  class B_Bphoton_E1 : public HD_ME_Base {
    Complex m_cL, m_cR;
  public:
    B_Bphoton_E1(const ATOOLS::Flavour_Vector& flavs,
                 const std::vector<int>& decayindices,
                 const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
  };

  class B3_Bphoton_E1 : public HD_ME_Base {
    Complex m_cL, m_cR;
  public:
    B3_Bphoton_E1(const ATOOLS::Flavour_Vector& flavs,
                  const std::vector<int>& decayindices,
                  const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
  };

  class QQ_QVQ_Spectator : public HD_ME_Base {
    ATOOLS::Flavour m_decayer;
  public:
    QQ_QVQ_Spectator(const ATOOLS::Flavour_Vector& flavs,
		     const std::vector<int>& decayindices,
		     const std::string& name);
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
    bool SetColorFlow(std::vector<ATOOLS::Particle*> outparts,
		      int n_q, int n_g, bool anti);
  };

  class QQ_PGG : public HD_ME_Base {
    double m_min_mass2;
  public:
    QQ_PGG(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& decayindices,
           const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    void SetModelParameters( GeneralModel _md );
  };

//   class P_GammaFF : public HD_ME_Base {
//     int m_phot,m_f1,m_f2;
//   public:
//     P_GammaFF(int,ATOOLS::Flavour *);
//     double operator()(const ATOOLS::Vec4D *);
//     void   operator()( 
//         const ATOOLS::Vec4D  * _p, 
//         std::vector<Complex> * _ampls_tensor, 
//         std::vector<std::pair<int,int> > * _indices,
//         int                    k0_n );
//   };
// 
//   class P_2PGamma : public HD_ME_Base {
//     int m_phot,m_p1,m_p2;
//   public:
//     P_2PGamma(int,ATOOLS::Flavour *);
//     void   operator()( 
//         const ATOOLS::Vec4D  * _p, 
//         std::vector<Complex> * _ampls_tensor, 
//         std::vector<std::pair<int,int> > * _indices,
//         int                    k0_n );
//     double operator()(const ATOOLS::Vec4D *);
//   };
// 
//   class P_P2Gamma : public HD_ME_Base {
//     int    m_phot1, m_phot2, m_p;
//     double m_mrho2, m_grho2;
//   public:
//     P_P2Gamma(int,ATOOLS::Flavour *);
//     void   operator()( 
//         const ATOOLS::Vec4D  * _p, 
//         std::vector<Complex> * _ampls_tensor, 
//         std::vector<std::pair<int,int> > * _indices,
//         int                    k0_n );
//     double operator()(const ATOOLS::Vec4D *);
//   };

//   class P_3P_DalitzDef : public HD_ME_Base {
//   private:
//     bool m_allpions, m_allsame;
//     int  m_pi0, m_pip, m_pim;
//   public:
//     P_3P_DalitzDef(int,ATOOLS::Flavour *);
//     double operator()(const ATOOLS::Vec4D *);
//   };

  class P_3P_Dalitz : public HD_ME_Base {
    // indices: P0=0, P1=1, P2=2, P3=3
    double m_const,m_liny,m_linx,m_linyphase,m_linxphase,
      m_quady,m_quadx,m_quadyphase,m_quadxphase,
      m_phaseliny,m_phaselinx,m_phasequady,m_phasequadx;
  public:
    P_3P_Dalitz(const ATOOLS::Flavour_Vector& flavs,
                const std::vector<int>& decayindices,
                const std::string& name):
      HD_ME_Base(flavs,decayindices,name) {};
    double operator()(const ATOOLS::Vec4D *);
    void Calculate(const ATOOLS::Vec4D_Vector& p, bool m_anti);
    void SetModelParameters( GeneralModel _md );
  };

}


#endif
