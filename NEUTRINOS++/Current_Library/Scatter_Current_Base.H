#ifndef NEUTRINOS_Current_Library_Scatter_Current_Base_H
#define NEUTRINOS_Current_Library_Scatter_Current_Base_H

#include "NEUTRINOS++/Tools/Form_Factor_Library.H"
#include "NEUTRINOS++/Tools/Propagator_Library.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

/*
#define DEFINE_CURRENT_GETTER(CLASS,TAG)				\
  DECLARE_GETTER(CLASS,TAG,Scatter_Current_Base,ME_Parameters);			\
  Scatter_Current_Base *ATOOLS::Getter<Scatter_Current_Base,ME_Parameters,CLASS>::	\
  operator()(const ME_Parameters &parameters) const			\
  { return new CLASS(parameters.flavs, parameters.indices, TAG); }
*/
namespace NEUTRINOS {
  /////////////////////////////////////////////////////////////////////////
  // This is pretty much copied from the HADRONS++ module.
  // The Scatter_Current_Base represents a vector current J^\mu in a scattering
  // process, which can be either weak or electromagnetic.  To store its
  // calculated values for each helicity combination, it inherits from
  // METOOLS::Spin_Structure.
  // Each class derived from this has to implement the Scatter_Current_Base::Calc
  // method for the actual calculation of the current.
  /////////////////////////////////////////////////////////////////////////
  class Scatter_Current_Base : public METOOLS::Spin_Structure<ATOOLS::Vec4C> {
  protected:
    const ATOOLS::Flavour_Vector& m_flavs;
    std::vector<double> m_masses;
    std::vector<int>    m_indices;
    std::string         m_name;
    std::string         m_path;
    std::map<std::string, Form_Factor_Base *> m_ffs;
    std::map<std::string, Propagator_Base *> m_ffprops;
    //std::map<std::string,std::map<size_t, Vec4C > > ProcessCurrents;
  public:
    Scatter_Current_Base(const ATOOLS::Flavour_Vector& flavs,
			 const std::vector<int>& indices,
			 const std::string& name);
    virtual ~Scatter_Current_Base() {
      while (!m_ffs.empty()) { delete m_ffs.begin()->second; m_ffs.erase(m_ffs.begin()); }
      while (!m_ffprops.empty()) { delete m_ffprops.begin()->second; m_ffprops.erase(m_ffprops.begin()); }
      //while (!ProcessCurrents.empty()) { ProcessCurrents.erase(ProcessCurrents.begin()); }
    }
    
    virtual void Calc(const ATOOLS::Vec4D_Vector& moms, METOOLS::XYZFunc * F) = 0;

    // virtual void Insert_ProcessType(string diagram, Vec4C value, std::vector<std::pair<int,int> >& spins)
    // {
    //   ProcessCurrents[diagram][GetNumber(spins)] = value;
    // }

    // virtual Vec4C Get_ProcessType(string diagram, const std::vector<int>& spins)
    // {
    //   return ProcessCurrents[diagram][GetNumber(spins)];
    // }

    inline std::string Name()                     const { return m_name; }
    inline void        SetName(std::string name)        { m_name = name; }
    inline int         GetN()                     const { return m_indices.size(); }
    inline const std::vector<int>& DecayIndices() const { return m_indices; }
  };
  
  //typedef ATOOLS::Getter_Function<Scatter_Current_Base,ME_Parameters>
  //Current_Getter_Function;
}

#endif
