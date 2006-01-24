#ifndef Spin_Correlation_Tensor_h
#define Spin_Correlation_Tensor_h

#include <vector>
#include <MyComplex.H>
#include <utility>
#include "Node.H"
#include "Spin_Density_Matrix.H"

namespace ATOOLS {

  class Spin_Correlation_Tensor
  {
  protected:

    static int m_k0_n;

    //! For an endnode: The complex value stored for this spin-combination
    Complex m_value;
    //! The particle the indices this node stores; -1 for endnodes which are only a complex number.
    int m_particle;

    //! pointers the follow-up nodes each corresponding to a different spin-combination for the current particle.
    std::vector<Spin_Correlation_Tensor*> *p_next;

    //! Returns the index-range for the current node (for example two for spin-1/2 particles).
    size_t GetIdxRange();
  public:
    //! Default constructor; sets p_next to NULL to avoid crashes.
    Spin_Correlation_Tensor():p_next(NULL){};
    virtual ~Spin_Correlation_Tensor();
     
    //! constructor
    Spin_Correlation_Tensor(
        std::vector<std::pair<int,int> >* particles, 
        std::vector<Complex>* Amplitudes,
        size_t pPos=0, 
        size_t aPos1=0, 
        size_t aPos2=0 );

    //! Test-constructor; creates a product of delta(i,j).
    Spin_Correlation_Tensor(double x, std::vector<int>* creation_list, size_t pos);

    //! Returns the trace; equivalent to soft-contracting with unit-matrices.
    Complex Trace();
    Complex Trace( Spin_Density_Matrix * sigma0 ); 

    //! Return the SDM for particle i.
    Spin_Density_Matrix GetSigma(int i);
    Spin_Density_Matrix GetSigma(int i, Spin_Density_Matrix * sigma0 ); 

    //! Contract particle i with SDM/DM
    void Contract(int i, Spin_Density_Matrix* SDM);

    //! Return a pointer to entry "number".
    Spin_Correlation_Tensor* operator()(const size_t &number);

    //! Removes a node without removing the follow-up nodes
    void Soft_Delete();

    //! Adds corresponding entries from SCT to current object.
    Spin_Correlation_Tensor& operator+=(Spin_Correlation_Tensor& SCT);

    //! Multiplies all follow-up elements by a Number
    Spin_Correlation_Tensor& operator*=(Complex& c);


    void Set_k0(int n) {m_k0_n = n;}
    int Get_k0_n() { return m_k0_n; }

    friend std::ostream &operator<<
      (std::ostream &ostr, Spin_Correlation_Tensor &sct);

    //! returns the number of particles the spin-correlation tensor currently administrates. Call without an argument!
    size_t GetDepth( size_t i=0 );
  };
  
  std::ostream &operator<<(std::ostream &ostr,
			   Spin_Correlation_Tensor &sct);
} // end of namespace ATOOLS
#endif
