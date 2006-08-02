#ifndef Current_Base_H
#define Current_Base_H

#include <iostream>
#include <string>
#include "MathTools.H"
#include "Vector.H"
#include "Flavour.H"
#include "Spin_Correlation_Tensor.H"
#include "ComplexVector.H"
#include "Getter_Function.H"
#include "Tools.H"

#ifndef SQRT_05
#define SQRT_05 0.70710678118654757
#endif

using namespace ATOOLS;

namespace HADRONS {

  class Current_Base {
    public:
    ComplexVec4D    * p_results;
    int               m_n, m_spin_combinations, m_k0n;
    ATOOLS::Flavour * p_flavs;
    ATOOLS::Vec4D   * p_moms;
    double          * p_masses,* p_masses2;
    int             * p_intspins; // integer spins
    int             * p_decayindices;  // needed later, when momenta are inserted
    
    std::string       m_name;
    std::string       m_path;

    const ComplexVec4D& operator[] (int i) const;
    virtual void Calc() = 0;
    void SetMomsAndk0n( const Vec4D * moms, const int k0n );

    Current_Base(const ATOOLS::Flavour* flavs, int n, int* decayindices, std::string name);
    virtual ~Current_Base();
    virtual void SetModelParameters( struct GeneralModel _md ) {};
    void SetPath( std::string s ) { m_path=s; }
  };
  
  void ContractCurrents(Current_Base* c1,
                        Current_Base* c2,
                        const Vec4D * moms,
                        const int k0n,
                        double mefactor,
                        std::vector<Complex>* ampls,
                        std::vector<std::pair<int,int> >* sctindices);

  typedef ATOOLS::Getter_Function<Current_Base,Flavour_Info> Current_Getter_Function;
}

std::ostream& operator<<(std::ostream& s, const HADRONS::Current_Base& cb);

#endif
