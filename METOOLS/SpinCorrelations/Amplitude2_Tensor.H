#ifndef METOOLS_SpinCorrelations_Amplitude2_Tensor_H
#define METOOLS_SpinCorrelations_Amplitude2_Tensor_H

#include <vector>
#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Matrix.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Phys/Particle.H"

namespace PHASIC {
  class Color_Function_Decay;
  typedef 
  std::vector<std::pair<METOOLS::Spin_Amplitudes*,Color_Function_Decay*> >
  DiagColVec;
}

namespace METOOLS {
  class Decay_Matrix;


  class Amplitude2_Tensor {
    std::vector<Amplitude2_Tensor*>* p_next;
    Complex m_value;
    ATOOLS::Particle* p_part;
    size_t m_nhel;

    Complex ContractRemaining(const std::vector<ATOOLS::Particle*> &parts,
                              const std::vector<int> &permutation,
                              size_t level,
                              const std::vector<Spin_Amplitudes> &diagrams,
                              const std::vector<std::vector<Complex> > &cmatrix,
                              std::vector<int> &spin_i, std::vector<int> &spin_j,
                              double factor) const;

    public:
    Amplitude2_Tensor(const std::vector<ATOOLS::Particle*>& parts, size_t level);
    Amplitude2_Tensor(const std::vector<ATOOLS::Particle*>& parts,
                      size_t level,
                      const PHASIC::DiagColVec& diagrams,
                      const std::vector<std::vector<Complex> >& colormatrix,
                      std::vector<int>& spin_i,
                      std::vector<int>& spin_j);
    Amplitude2_Tensor(const std::vector<ATOOLS::Particle*>& parts,
                      const std::vector<int>& permutation,
                      size_t level,
                      const std::vector<Spin_Amplitudes>& diagrams,
                      const std::vector<std::vector<Complex> >& cmatrix,
                      std::vector<int>& spin_i, std::vector<int>& spin_j);
    ~Amplitude2_Tensor();

    void Contract(const Amplitude2_Matrix* D);
    Amplitude2_Matrix ReduceToMatrix(const ATOOLS::Particle* leftover) const;
    void Add(const Amplitude2_Tensor* amp, const Complex& factor);
    void Multiply(const Complex& factor);
    Complex Trace() const;
    bool Contains(const ATOOLS::Particle* part) const;

    inline const std::vector<Amplitude2_Tensor*> Next() const {return *p_next;}

    void Print(std::ostream& ostr, std::string label) const;
    friend std::ostream& operator<<(std::ostream&, const Amplitude2_Tensor&);

    static bool SortCrit(const std::pair<ATOOLS::Particle*, size_t>& p1,
                         const std::pair<ATOOLS::Particle*, size_t>& p2);
  };


}


#endif
