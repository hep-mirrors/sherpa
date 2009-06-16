#ifndef ATOOLS_Phys_Spin_Correlation_Tensor_H
#define ATOOLS_Phys_Spin_Correlation_Tensor_H

#include <vector>
#include <set>
#include "ATOOLS/Math/MyComplex.H"
#include <utility>
#include "ATOOLS/Org/Node.H"
#include "ATOOLS/Phys/Spin_Density_Matrix.H"


namespace ATOOLS {
  
  class scmode {
  public:
    enum code {
      None     = 0,
      Diagonal = 1,
      Full     = 2 };
  };

  class Spin_Correlation_Tensor
  {
  protected:
    //! The gauge used in the Helicity Formalism.
    static int m_k0_n;
    //! A list of particles which can contribute indices to a Spin_Correlation_Tensor.
    static std::set<int> m_possible_particles;
    //! Level of spin correlations used in the current run.
    static scmode::code m_mode;
    //! Flag if a Spin_Correlation_Tensor already exists/existed in the current run. If true, changing m_mode becomes impossible.
    static bool m_created;

    // Member-variables for each node.
  protected:
    //! In end-nodes, the complex value of the index is stored here.
    Complex m_value;
    //! Number of the particle the node belongs to; -1 for end-nodes.
    int m_particle;
    //! Pointers to the follow-up nodes. NULL for end-nodes.
    std::vector<Spin_Correlation_Tensor*> *p_next;

    // Internal helper-methods.
  protected:
    //! Return the maximum spin; usefull if m_mode == scmode::Full.
    size_t GetIdxRange();
    //! Removes a node without removing the follow-up nodes
    void Soft_Delete();
    //! Adds corresponding entries from SCT to current object.
    Spin_Correlation_Tensor& operator+=(Spin_Correlation_Tensor& SCT);
    //! Return a pointer to entry "number".
    Spin_Correlation_Tensor* operator()(const size_t &number);
    //! Multiplies all follow-up elements by a Number
    Spin_Correlation_Tensor& operator*=(Complex& c);


    // Methods for accessing the global spin-correlation settings.
  public:

    static void AddPossibleParticles( std::set<kf_code> * flavs  );
    static void AddPossibleParticle( kf_code flav );
    static bool PossibleParticle( kf_code flav );
    static void PrintPossibleParticles();

    static void SetMode(scmode::code mode);
    static scmode::code Mode();

    //! Tells the Spin_Correlation_Tensor which gauge was used for computing the entries.
    static void Set_k0(int n) {m_k0_n = n;}
    //! Returns the gauge used for computing the entries.
    static int Get_k0_n() { return m_k0_n; }


    // Methods for accessing an individual Spin_Correlation_Tensor.
  public:
    //! Default constructor; sets p_next to NULL to avoid crashes.
    Spin_Correlation_Tensor():p_next(NULL){};
    //! Constructor used in the HADRONS module.
    Spin_Correlation_Tensor(
        std::vector<std::pair<int,int> >* particles, 
        std::vector<Complex>* Amplitudes,
        size_t pPos=0, 
        size_t aPos1=0, 
        size_t aPos2=0,
        size_t add=1 );
    virtual ~Spin_Correlation_Tensor();

    //! Returns the trace over the tensor structure (contracting with unit-matrices).
    Complex Trace();
    //! Contracts the particle index in the top-node with sigma0 and then returns the trace over the remaining indices.
    Complex Trace( Spin_Density_Matrix * sigma0 ); 

    //! Returns the Spin_Density_Matrix for particle i.
    Spin_Density_Matrix GetSigma(int i);
    //! Contracts the top-node particle with sigma0 and returns the resulting Spin_Density_Matrix for particle i.
    Spin_Density_Matrix GetSigma(int i, Spin_Density_Matrix * sigma0 ); 

    //! Contract particle i with SDM. This is a hard contract; particle i will no longer contribute indices to the tensor.
    void Contract(int i, Spin_Density_Matrix* SDM);

    //! An outstream operator used for debugging.
    friend std::ostream &operator<<
      (std::ostream &ostr, Spin_Correlation_Tensor &sct);

    //! returns the number of particles the spin-correlation tensor currently administrates. Call without an argument!
    size_t GetDepth( size_t i=0 );
  };
  /*!
    \class Spin_Correlation_Tensor
    \brief Used for transporting spin informations through the different stages of event generation.

    Spin_Correlation_Tensor is the main class for implementing the spin correlation 
    algorithm proposed in hep-ph/0110108 in SHERPA. The class represents the correlation
    tensor \f$T_{(ij)_1\dots (ij)_N}=A_{i_1\dots i_N}A^*_{j_1\dots j_N}\f$
    where \f$ A_{i_1 \dots i_N} \f$ is the amplitude for for the case that particle 
    \f$n\f$ has a spin orientation \f$i_n\f$ and the number of particles is N.

    The tensor structure is implemented as a tree structure. Each level of the tree 
    corresponds to the indices \f$(ij)_n\f$ for one particle \f$n\f$. The end-nodes 
    (leafs) do not correspond to specific particles. They only store the value the 
    tensor has for the given combination of indices.

    To save memory and CPU-time, spin informations are not created for all particles in
    the event. Instead, there is a global positive-list 
    Spin_Correlation_Tensor::m_possible_particles which stores all particle types for which
    spin informations are to be created. The list is filled by SHERPA. Interesting particles
    are unstable particles created by the hard process which will decay in a later stage
    of the program (in the HADRONS module), namely the tau-fermion.
  */
  /*!
    \fn  static void Spin_Correlation_Tensor::AddPossibleParticles( std::set<kf_code> * flavs  )
    \brief Stores the list of kf-codes given as particle types for which spin informations are created.
  */
  /*!
    \fn static void Spin_Correlation_Tensor::AddPossibleParticle( kf_code flav )
    \brief Enables creation of spin information for the particle type flav.
  */
  /*!
    \fn static bool Spin_Correlation_Tensor::PossibleParticle( kf_code flav )
    \brief Tells the program if spin informations for particle type flav are created.
  */
  /*!
    \fn static void Spin_Correlation_Tensor::PrintPossibleParticles()
    \brief Prints the list of all particle types for which spin informations are stored by an SCT.
  */
  /*!
    \fn static void Spin_Correlation_Tensor::SetMode(scmode::code mode)
    \brief Sets the mode of spin information to be used by the program.

    Creating full spin informations for a hard process and consequently using it in the
    HADRON decays is a very time consuming process. Using the full information can easily
    increase the runtime of the program by a factor of 2-3. Additionally, many observables
    are not or only marginally sensitive to spin. 
    Therefore, Spin_Correlation_Tensor 
    has a static variable Spin_Correlation_Tensor::m_mode which can reduce the amount of
    spin informations to be used by the program. The three possible values to set m_mode
    to are:
    - scmode::Off - No spin informations are generated at all. This is the default value.
    - scmode::Diagonal - Only elements for which \f$i_n = j_n\f$ are used in the correlation
    tensor. This mode is sufficient for most observables.
    - scmode::Full - The whole tensor is created and used. This mode is the most memory- and
    CPU-consuming one.

    The mode being used is set by SHERPA before event generation. It is determined by the
    value of the command-line parameter "SPIN_CORRELATIONS".
  */
  /*!
    \fn static scmode::code Spin_Correlation_Tensor::Mode()
    \brief Returns the mode of spin correlations currently used.
  */

} // end of namespace ATOOLS
#endif
