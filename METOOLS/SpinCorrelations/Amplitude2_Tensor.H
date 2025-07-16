#ifndef METOOLS_SpinCorrelations_Amplitude2_Tensor_H
#define METOOLS_SpinCorrelations_Amplitude2_Tensor_H

#include <vector>
#include <memory>
#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Matrix.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Vector.H"

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
                              std::vector<int> &spin_i, std::vector<int> &spin_j,
                              double factor) const;

  public:
    Amplitude2_Tensor(const std::vector<ATOOLS::Particle*>& parts, size_t level);
    Amplitude2_Tensor(const std::vector<ATOOLS::Particle*>& parts,
                      size_t level,
                      const std::vector<Spin_Amplitudes*>& diagrams,
                      std::vector<int>& spin_i,
                      std::vector<int>& spin_j);
    Amplitude2_Tensor(const std::vector<ATOOLS::Particle*>& parts,
                      const std::vector<int>& permutation,
                      size_t level,
                      const std::vector<Spin_Amplitudes>& diagrams,
                      std::vector<int>& spin_i, std::vector<int>& spin_j);
    Amplitude2_Tensor(const Amplitude2_Tensor& other);
    ~Amplitude2_Tensor();

    void Contract(const Amplitude2_Matrix* D);
    Amplitude2_Matrix ReduceToMatrix(const ATOOLS::Particle* leftover) const;
    void Add(const Amplitude2_Tensor* amp, const Complex& factor);
    void Multiply(const Complex& factor);
    /**
    * @param D The (pointer to) the Amplitude2_Matrix with which the Amplitude2_Tensor should be multiplied with.
    *
    * method multiplies the entries of an Amplitude2_Matrix D to all entries of the Ampltitude2_Tensor with the same
     * helicity
    */
    void Multiply(const Amplitude2_Matrix* D);
    /**
    * @return unpolarized matrix element squared.
    *
    * method for summing up all entries of the Amplitude2_Tensor to receive the unpolarized result
    */
    Complex Sum();
    Complex Trace() const;
    bool Contains(const ATOOLS::Particle* part) const;
    /**
    * @param part_number number of the searched particle according to Sherpa's internal particle ordering
    * @param level level of Amplitude2_Tensor which describes the searched particle, only for recursive calls;
     * default: 1
    *
    * @return Pair of level of and the pointer to searched particle; if particle is not found, pair of last level and
    *         an empty pointer is returned
    *
    * method for searching, whether (and where) a particle with a specific particle number is included in the Amplitude2_Tensor
    */
    std::pair<const int, const ATOOLS::Particle*> Search(int part_number, int level=1) const;
    /**
    * @param num current number of particles, necessary for recursive call of the function; default: 0
    *
    * @return The number of particles in the Amplitude2_Tensor
    *
    * method for determining the number of particles in the Amplitude2_Tensor
    */
    int NumberParticles(int num=0) const;
    /**
    * @return bool; false if the end of an Amplitude2_Tensor branch is reached
    * method for determining whether the end of an Amplitude2_Tensor branch is reached
    *
    */
    bool IsP_Next() const;
    void UpdateParticlePointers(const std::map<ATOOLS::Particle*,ATOOLS::Particle*>& pmap);

    inline const std::vector<Amplitude2_Tensor*> Next() const {return *p_next;}
    inline const ATOOLS::Particle CurrentParticle() const { return *p_part; }
    inline const Complex Value() const { return m_value; }
    inline const size_t SpinDegreesofFreedom() const { return m_nhel; }

    void Print(std::ostream& ostr, std::string label) const;
    friend std::ostream& operator<<(std::ostream&, const Amplitude2_Tensor&);

    static bool SortCrit(const std::pair<ATOOLS::Particle*, size_t>& p1,
                         const std::pair<ATOOLS::Particle*, size_t>& p2);

    /**
    * @param coeff vector containing the transformation coefficients,
     * first level of vectors contains the coefficients for the different particles, second level for the different
     * final helicities and third level for the different old helicities in dependency of the final ones one level
     * above; first level should be ordered according to the particle ordering in Amplitude2_Tensor which should be
     * transformed; second level and third level according to its helicity ordering
    * @param conj_coeff vector containing the coefficients for the transformation of the complex conjugate
     * matrix elements, form analog to coeff
    * @param level level describes the level of the Amplitude2_tensor the method is working on at the moment,
     * necessary for recursive call of the method
     * default: 0
    * @param coeff_tmp vector containing the transformation coefficients,
     * form analog to coeff, only for recursive calls to save selected transformation coefficients necessary for
     * calculating the transformation of single entries in Amplitude2_Tensor
     * default: std::vector<std::vector<std::vector<Complex>>>()
    * @param conj_coeff_tmp vector containing the selected coefficients for the transformation of the complex conjugate
     * matrix elements, form analog to coeff, for recursive calls only
     * default: std::vector<std::vector<std::vector<Complex>>>()
    * @param old_amps  pointer to Amplitude2_Tensor before transformation, default: NULL
    *
    * method for transforming the Amplitude2_Tensor to another basis for polarization definition,
    * for that transformation coefficients for the matrix elements and their complex conjugate must be handed in;
    * the method iterates two times through an Amplitude2_Tensor, first time to reach the entries which should be
    * changed and to select the coefficients which are important for calculation of the current new entry
    * (contained in second level of coeff); during the second time the new entry is calculated by multiplying
    * every entry of the old tensor by the corresponding coefficients
    */
    void PolBasisTrafo(const std::vector<std::vector<std::vector<Complex> > >& coeff,
                       const std::vector<std::vector<std::vector<Complex> > >& conj_coeff, int level=0,
                       std::vector<std::vector<std::vector<Complex> > > coeff_tmp=std::vector<std::vector<std::vector<Complex> > >(),
                       std::vector<std::vector<std::vector<Complex> > > conj_coeff_tmp=std::vector<std::vector<std::vector<Complex> > >(),
                       Amplitude2_Tensor* old_amps=NULL);
  };

  using Amplitude2_Tensor_SP = std::shared_ptr<METOOLS::Amplitude2_Tensor>;
}


#endif
