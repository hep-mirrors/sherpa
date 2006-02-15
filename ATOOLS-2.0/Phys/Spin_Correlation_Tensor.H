#ifndef Spin_Correlation_Tensor_h
#define Spin_Correlation_Tensor_h

#include <vector>
#include <set>
#include <MyComplex.H>
#include <utility>
#include "Node.H"
#include "Spin_Density_Matrix.H"


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
    // Statics; global parameters.
  protected:
    static int m_k0_n;
    static std::set<int> m_possible_particles;
    static scmode::code m_mode;
    static bool m_created;

    // Member-variables for each node.
  protected:
    Complex m_value;
    int m_particle;
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

    static void AddPossibleParticles( std::set<ATOOLS::kf::code> * flavs  );
    static void AddPossibleParticle( ATOOLS::kf::code flav );
    static bool PossibleParticle( ATOOLS::kf::code flav );
    static void PrintPossibleParticles();

    static void SetMode(scmode::code mode);
    static scmode::code Mode();

    static void Set_k0(int n) {m_k0_n = n;}
    static int Get_k0_n() { return m_k0_n; }


    // Methods for accessing an individual Spin_Correlation_Tensor.
  public:
    //! Default constructor; sets p_next to NULL to avoid crashes.
    Spin_Correlation_Tensor():p_next(NULL){};
    //! constructor for HADRONS
    Spin_Correlation_Tensor(
        std::vector<std::pair<int,int> >* particles, 
        std::vector<Complex>* Amplitudes,
        size_t pPos=0, 
        size_t aPos1=0, 
        size_t aPos2=0,
        size_t add=1 );
    virtual ~Spin_Correlation_Tensor();

    //! Returns the trace; equivalent to soft-contracting with unit-matrices.
    Complex Trace();
    Complex Trace( Spin_Density_Matrix * sigma0 ); 

    //! Return the SDM for particle i.
    Spin_Density_Matrix GetSigma(int i);
    Spin_Density_Matrix GetSigma(int i, Spin_Density_Matrix * sigma0 ); 

    //! Contract particle i with SDM/DM
    void Contract(int i, Spin_Density_Matrix* SDM);

    friend std::ostream &operator<<
      (std::ostream &ostr, Spin_Correlation_Tensor &sct);

    //! returns the number of particles the spin-correlation tensor currently administrates. Call without an argument!
    size_t GetDepth( size_t i=0 );
  };
  
  std::ostream &operator<<(std::ostream &ostr,
			   Spin_Correlation_Tensor &sct);
} // end of namespace ATOOLS
#endif
