#ifndef N_Parton_CDBG_T_H
#define N_Parton_CDBG_T_H

#include "CDBG_Amplitude.H"

namespace ATOOLS { class Expression; }

namespace EXTRAXS {

  typedef std::vector<Complex> Complex_Vector;

  class N_Parton_CDBG_T {
  private:

    size_t m_nin, m_nout;

    CDBG_Amplitude m_ampl;

    Vec4D_Vector   m_moms;
    Flavour_Vector m_flavs;
    Int_Vector     m_imap;
  
    double m_sf, m_norm;
    int    m_mode, m_tests, m_ng;

    std::map<Int_Vector,size_t> m_omap;
    std::vector<Complex_Vector> m_colfs, m_cress;
    std::vector<Int_Vector>     m_cords, m_evals;

    double Factorial(const double &n);

    bool CheckOrder(const Int_Vector &co);
    bool CheckOrder(const Int_Vector &cor,const Int_Vector &coa);
    void ReOrder(Int_Vector &ci,Int_Vector &cj);

    bool AddColors(ATOOLS::Expression *const ex,const bool &b,
		   const std::vector<int> &cor);
    Complex ColorFactor(const std::vector<int> &cor,
			const std::vector<int> &coa);
    void Construct(Int_Vector ic,Int_Vector left);
    void Construct();

    void OrderFlavours();

  public:

    // constructor
    N_Parton_CDBG_T(const size_t nin,const size_t nout,
		    const std::vector<ATOOLS::Flavour> &flavs,
		    const std::vector<std::string> &models);

    // destructor
    ~N_Parton_CDBG_T();

    // member functions
    double Differential(const std::vector<ATOOLS::Vec4D> &momenta);

    bool GaugeTest(std::vector<ATOOLS::Vec4D> momenta);

    // inline functions
    inline void SetMode(const int &mode)    { m_mode=mode;   }
    inline void SetTests(const int &tests)  { m_tests=tests; }

    inline Flavour_Vector Flavs() const { return m_flavs; }

  };// end of class N_Parton_CDBG_T

}// end of namespace EXTRAXS

#endif

