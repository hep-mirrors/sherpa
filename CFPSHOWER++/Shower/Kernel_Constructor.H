#ifndef CFPShower_Shower_Kernel_Constructor_H
#define CFPShower_Shower_Kernel_Constructor_H

#include "CFPSHOWER++/Shower/Kernels.H"

namespace MODEL {
  class Model_Base;
  class Running_AlphaS;
  class Running_AlphaQED;
  class Single_Vertex;
}

namespace PDF {
  class PDF_Base;
  class ISR_Handler;
}


namespace CFPSHOWER {
  class Shower;
  
  class Kernel_Constructor {
  private:
    Shower * p_shower;
    double   m_asfactor[2], m_muR2factor, m_muF2factor, m_pdfmin[2];
    int      m_kinscheme, m_kfactor, m_cplscheme, m_MEcorrs, m_SForder;
    
    MODEL::Model_Base           * p_model;
    PDF::PDF_Base               * p_pdf[2];
    MODEL::Running_AlphaS       * p_as;
    MODEL::Running_AlphaQED     * p_alpha;
    
    std::map<ATOOLS::Flavour, Kernels *> * p_kernels[5];
    std::vector<std::list<std::vector<size_t> > > m_permutations;

    bool InitializeKernels(MODEL::Model_Base * const model);
    void MakeKernelsFromVertices(MODEL::Model_Base * const model);
    void MakeKernelsFromKernels(const kernel_type::code & type);
    void MakeKernels(ATOOLS::Flavour_Vector & flavs,
		     const kernel_type::code & type);
    void DeleteAllKernelsWithWrongOrder();
    void MakePermutations();
    void GeneratePermutation(std::vector<size_t> perm,const size_t & pos,
			     const size_t & length);
  public:
    Kernel_Constructor(Shower * shower);
    ~Kernel_Constructor();

    inline std::map<ATOOLS::Flavour, Kernels *> * operator()(const size_t & i) {
      return p_kernels[i];
    }
    bool Init(MODEL::Model_Base * const model,PDF::ISR_Handler * const isr);
    void PrintKernels(const bool & full=false);
  };
}
#endif
