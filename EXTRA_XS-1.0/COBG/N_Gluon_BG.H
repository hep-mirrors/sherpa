#ifndef N_Gluon_BG_H
#define N_Gluon_BG_H

#include "BG_Amplitude.H"

namespace EXTRAXS {

  typedef std::vector<Complex> Complex_Vector;

  class N_Gluon_BG {
  private:

    size_t m_nin, m_nout;

    BG_Amplitude m_ampl;

    std::vector<ATOOLS::Vec4D> m_moms;
  
    double m_as, m_sf;
    int    m_mode, m_tests;

    std::map<Int_Vector,size_t> m_omap;
    std::vector<Complex_Vector> m_colfs, m_cress;
    std::vector<Int_Vector>     m_cords, m_evals;

    double Factorial(const double &n);

    void ReOrder(Int_Vector &ci,Int_Vector &cj);

    void Construct(Int_Vector ic,Int_Vector left);
    void Construct();

  public:

    // constructor
    N_Gluon_BG(const size_t nin,const size_t nout,
	       const std::vector<ATOOLS::Flavour> &flavs);

    // destructor
    ~N_Gluon_BG();

    // member functions
    bool GenerateConfiguration();

    double Differential(const std::vector<ATOOLS::Vec4D> &momenta);

    bool GaugeTest(std::vector<ATOOLS::Vec4D> momenta);

    // inline functions
    inline void SetMode(const int &mode)    { m_mode=mode;   }
    inline void SetTests(const int &tests)  { m_tests=tests; }

  };// end of class N_Gluon_BG

}// end of namespace EXTRAXS

#endif

