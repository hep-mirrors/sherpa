#ifndef HADRONS_Current_Library_Current_Base_H
#define HADRONS_Current_Library_Current_Base_H

#include <iostream>
#include <string>
#include <iomanip>
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "HADRONS++/Main/Tools.H"
#include "METOOLS/Main/Spin_Structure.H"

#define DEFINE_CURRENT_GETTER(CLASS,NAME,TAG)                                   \
DECLARE_GETTER(NAME,TAG,Current_Base,Flavour_Info);                             \
Current_Base* NAME::operator()(const Flavour_Info &parameters) const            \
{ return new CLASS(parameters.flavs, parameters.n, parameters.indices,TAG); }

namespace METOOLS {
  typedef Spin_Structure<ATOOLS::Vec4C> Spin_Currents;
}

namespace HADRONS {

  class Current_Base : public METOOLS::Spin_Currents {
  protected:
    int               m_n;
    bool              m_anti;
    ATOOLS::Flavour * m_flavs;
    double          * p_masses;
    int             * p_i;  // saving them for later
    
    std::string       m_name;
    std::string       m_path;
    
    virtual void Calc(const ATOOLS::Vec4D * moms) = 0;
  public:
    Current_Base(ATOOLS::Flavour* flavs, int n, int* decayindices, std::string name);
    virtual ~Current_Base();
    virtual void SetModelParameters( struct GeneralModel _md ) {};
    void SetAnti(bool anti) { m_anti = anti; }
    
    std::string Name() const              { return m_name; }
    void        SetName(std::string name) { m_name = name; }
    std::string Path() const              { return m_path; }
    void        SetPath(std::string path) { m_path = path; }
    int         GetN() const              { return m_n; }
    int       * DecayIndices() const      { return p_i; }

    friend void ContractCurrents( Current_Base* c1,
				  Current_Base* c2,
                                  const ATOOLS::Vec4D * moms,
				  Complex mefactor,
				  METOOLS::Spin_Amplitudes* amps);
    friend std::ostream& operator<<(std::ostream& s, const Current_Base& cb);
  };
  

  typedef ATOOLS::Getter_Function<Current_Base,Flavour_Info> Current_Getter_Function;
}


#endif
