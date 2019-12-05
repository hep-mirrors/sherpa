#ifndef DIRE__Shower__Shower_H
#define DIRE__Shower__Shower_H

#include "DIRE/Tools/Parton.H"
#include "DIRE/Shower/Kernel.H"
#include "DIRE/Tools/Amplitude.H"

namespace ATOOLS {
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

namespace DIRE {

  class Cluster_Definitions;

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
      int m_acc, m_nlo;
      std::vector<int> m_skip;
      JetVeto_Args(ATOOLS::Cluster_Amplitude *const ampl,
		   const double &jcv,const size_t &n):
	p_ampl(ampl), m_jcv(jcv), m_acc(0), m_skip(n,0) {}
    };// end of struct JetVeto_Args
    
    typedef std::map<ATOOLS::Flavour,Kernel_Vector> SKernel_Map;

    typedef std::map<ATOOLS::Flavour,Kernel*> EKernel_Map;
    typedef std::map<ATOOLS::Flavour,EKernel_Map> SEKernel_Map;
    typedef std::map<int,SEKernel_Map> Kernel_Map;

  private:

    // set up a C array wrapper used to hold a table of cumulative integrals in
    // GeneratePoint
    // rows: kernels allowed by the splitter
    // cols: possible spectators
    // NOTE: by using (encapsulated) C arrays we can avoid using vectors of
    // vectors in the hot-spot function GeneratePoint
    struct CumulativeIntegralTable {
      CumulativeIntegralTable() : nrows {0}, ncols {0} {};
      CumulativeIntegralTable(int _nrows, int _ncols) :
        nrows{_nrows}, ncols {_ncols}
      {
        sumsizes = new int[nrows];
        spects = new size_t[nrows*ncols];
        sums = new double[nrows*ncols];
      }
      CumulativeIntegralTable& operator=(CumulativeIntegralTable&& other)
      {
        if (this == &other)
          return *this;
        if (nrows > 0) {
          delete[] sumsizes;
          delete[] spects;
          delete[] sums;
        }
        nrows = other.nrows;
        ncols = other.ncols;
        sumsizes = other.sumsizes;
        spects = other.spects;
        sums = other.sums;
        other.nrows = 0;
        return *this;
      }
      ~CumulativeIntegralTable() {
        if (nrows > 0) {
          delete[] sumsizes;
          delete[] spects;
          delete[] sums;
        }
      }
      void Clear(int row) { sumsizes[row] = 0; }
      int Size(int row) { return sumsizes[row]; }
      double Sum(int row, int col) { return sums[row*ncols + col]; }
      size_t Spect(int row, int col) { return spects[row*ncols + col]; }
      double LastSum(int row) { return sums[row*ncols + sumsizes[row] - 1]; }
      double LastSpect(int row) { return spects[row*ncols + sumsizes[row] - 1]; }
      void AppendSumAndSpect(int row, double sum, size_t spect) {
        const auto idx = row*ncols + sumsizes[row];
        sums[idx] = sum;
        spects[idx] = spect;
        ++(sumsizes[row]);
      }
      int nrows, ncols;
      int* sumsizes;
      size_t* spects;
      double* sums;
    };

    CumulativeIntegralTable m_sums;

    MODEL::Model_Base *p_model;
    PDF::PDF_Base     *p_pdf[2];

    MODEL::Running_AlphaS *p_as;

    Cluster_Definitions *p_cluster;

    ATOOLS::Variation_Weights *p_vars;

    Kernel_Vector m_cks;
    SKernel_Map   m_sks;
    Kernel_Map    m_kmap;

    double m_tmin[2], m_rewtmin, m_cplfac[2], m_rsf, m_fsf;
    double m_weight, m_oef, m_pdfmin[2], m_rcf, m_tcef;

    int m_kin, m_kfac, m_cpl, m_mec;

    unsigned int m_maxem, m_maxpart, m_maxrewem;

    void AddKernel(Kernel *const k);

    void AddWeight(const Amplitude &a,const double &t);

    void SetKinSpects(Amplitude &a);

    Splitting GeneratePoint(Parton &p,const double &t,
			    const int &cm,const unsigned int &nem);
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

    bool Init(MODEL::Model_Base *const,
	      PDF::ISR_Handler *const);

    void SetMS(const ATOOLS::Mass_Selector *const ms);

    int Evolve(Amplitude &a,double &w,unsigned int &nem);

    double GetXPDF(const double &x,const double &Q2,
		   const ATOOLS::Flavour &fl,const int b) const;

    Kernel *GetKernel(const Splitting &s,const int mode) const;

    int RemnantTest(Parton *const c,const ATOOLS::Vec4D &p);

    inline MODEL::Model_Base *Model() const { return p_model; }

    inline MODEL::Running_AlphaS *const &AlphaS() const { return p_as; }

    inline double TMin(const int i) const { return m_tmin[i]; }

    inline double CplFac(const int i) const { return m_cplfac[i]; }

    inline double PDFMin(const int i) const { return m_pdfmin[i]; }

    inline double GetWeight() const { return m_weight; }

    inline int KFactorScheme() const  { return m_kfac; }
    inline int CouplingScheme() const { return m_cpl;  }

    inline int KinematicsScheme() const { return m_kin; }

    inline int MECorrection() const { return m_mec; }

    inline double MuR2Factor() const { return m_rsf; }

    inline void SetVariations(ATOOLS::Variation_Weights *const vars)
    { p_vars=vars; }

    inline ATOOLS::Variation_Weights *Variations() const
    { return p_vars; }

  };//end of class Shower

}// end of namespace DIRE

#endif
