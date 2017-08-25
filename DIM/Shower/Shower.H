#ifndef DIM__Shower__Shower_H
#define DIM__Shower__Shower_H

#include "DIM/Tools/Parton.H"
#include "DIM/Shower/Kernel.H"
#include "DIM/Tools/Amplitude.H"

namespace ATOOLS {
  class Default_Reader;
  class Variation_Parameters;
  class Variation_Weights;
}
namespace MODEL {
  class Model_Base;
  class Running_AlphaS;
}
namespace PDF {
  class PDF_Base;
  class ISR_Handler;
}

namespace DIM {

  class Gamma;

  class Shower {
  public:

    struct Reweight_Args {
      Splitting *m_s;
      int m_acc;
      Reweight_Args(Splitting *const s,const int acc):
	m_s(s), m_acc(acc) {}
    };// end of struct Reweight_Args
    
    struct JetVeto_Args {
      ATOOLS::Cluster_Amplitude *p_ampl;
      double m_jcv;
      int m_acc;
      JetVeto_Args(ATOOLS::Cluster_Amplitude *const ampl,
		   const double &jcv):
	p_ampl(ampl), m_jcv(jcv), m_acc(0) {}
    };// end of struct JetVeto_Args
    
    typedef std::map<ATOOLS::Flavour,Kernel_Vector> SKernel_Map;

    typedef std::map<ATOOLS::Flavour,Kernel*> EKernel_Map;
    typedef std::map<ATOOLS::Flavour,EKernel_Map> SEKernel_Map;
    typedef std::map<int,SEKernel_Map> Kernel_Map;

  private:

    MODEL::Model_Base *p_model;
    PDF::PDF_Base     *p_pdf[2];

    MODEL::Running_AlphaS *p_as;

    Gamma *p_gamma;

    ATOOLS::Variation_Weights *p_vars;

    Kernel_Vector m_cks;
    SKernel_Map   m_sks;
    Kernel_Map    m_kmap;

    double m_tmin[2], m_cplfac[2], m_rsf, m_fsf;
    double m_weight, m_oef, m_pdfmin[2], m_rcf;

    int m_kfac, m_cpl;

    unsigned int m_maxem, m_maxrewem;

    Splitting m_s;

    void AddKernel(Kernel *const k);

    void AddWeight(const Amplitude &a,const double &t);

    Splitting GeneratePoint(Parton &p,const double &t,
			    const unsigned int &nem);
    Splitting GeneratePoint(const Amplitude &a,const double &t,
			    const unsigned int &nem);

    double Reweight(ATOOLS::Variation_Parameters *params,
		    ATOOLS::Variation_Weights *weights,
		    const Reweight_Args &a);
    double GetWeight(ATOOLS::Variation_Parameters *params,
		     ATOOLS::Variation_Weights *weights,
		     std::vector<double> &v);
    double VetoWeight(ATOOLS::Variation_Parameters *params,
		      ATOOLS::Variation_Weights *weights,
		      JetVeto_Args &a);

  public:

    Shower();

    ~Shower();

    bool Init(MODEL::Model_Base *const model,
	      PDF::ISR_Handler *const isr,
	      ATOOLS::Default_Reader *const read);

    void SetMS(ATOOLS::Mass_Selector *const ms);

    int Evolve(Amplitude &a,unsigned int &nem);

    double GetXPDF(const double &x,const double &Q2,
		   const ATOOLS::Flavour &fl,const int b) const;

    Kernel *GetKernel(const Splitting &s,const int mode) const;

    inline MODEL::Model_Base *Model() const { return p_model; }

    inline MODEL::Running_AlphaS *const &AlphaS() const { return p_as; }

    inline double TMin(const int i) const { return m_tmin[i]; }

    inline double CplFac(const int i) const { return m_cplfac[i]; }

    inline double PDFMin(const int i) const { return m_pdfmin[i]; }

    inline double GetWeight() const { return m_weight; }

    inline int KFactorScheme() const  { return m_kfac; }
    inline int CouplingScheme() const { return m_cpl;  }

    inline double MuR2Factor() const { return m_rsf; }

    inline const Splitting &LastSplitting() const { return m_s; }

    inline void SetGamma(Gamma *const gamma) { p_gamma=gamma; }

    inline double OEF() const { return m_oef; }

    inline void SetVariations(ATOOLS::Variation_Weights *const vars)
    { p_vars=vars; }

    inline ATOOLS::Variation_Weights *Variations() const
    { return p_vars; }

  };//end of class Shower

}// end of namespace DIM

#endif
