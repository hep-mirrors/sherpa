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

#define DEFINE_CURRENT_GETTER(CLASS,NAME,TAG)                                  \
DECLARE_GETTER(NAME,TAG,Current_Base,ME_Parameters);                           \
Current_Base* NAME::operator()(const ME_Parameters &parameters) const          \
{ return new CLASS(parameters.flavs, parameters.indices, TAG); }

namespace HADRONS {

  /** 
   * Represents a vector current in a decay process.
   * To store its calculated values for each helicity combination, it inherits
   * from METOOLS::Spin_Structure.
   * Each class derived from this has to implement the Current_Base::Calc
   * method for the actual calculation of the current.
   */
  class Current_Base : public METOOLS::Spin_Structure<ATOOLS::Vec4C> {
  protected:
    /// Reference to all flavours involved in the decay.
    /// This usually contains more flavours then involved in the current itself
    /// and thus needs Current_Base::p_i for a mapping
    const ATOOLS::Flavour_Vector& m_flavs;
    /// Array of masses of the flavours involved in the current
    /// (in the internal ordering)
    double          * p_masses;
    /// Index mapping from external flavours/momenta to current internal ones
    std::vector<int> p_i;

    /// name of the current
    std::string       m_name;
    /// Path to the directory which contains the decay channel's file
    std::string       m_path;
    
  public:
    /** 
     * Constructor which initialises all fields of the class.
     * 
     * @param flavs All flavours involved in the decay (in external order).
     * @param decayindices Index mapping from external flavours/momenta to
     * current internal ones (cf. Current_Base::p_i)
     * @param name Name of the current
     */
    Current_Base(const ATOOLS::Flavour_Vector& flavs,
                 const std::vector<int>& decayindices,
                 const std::string& name);
    virtual ~Current_Base();

    /** 
     * Pure virtual function, which requires the implementation of the
     * current calculation in all derived classes.
     * When implementing such a calculation, one can use Spin_Structure::Insert
     * to store the calculated values for all helicities combination.
     * 
     * @param momenta All momenta in the decay in external order
     * @param anti Whether to consider the charge conjugated process
     */
    virtual void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti) = 0;

    /** 
     * This method sets the parameters for the decay current.
     * Every subclass of Current_Base \b must have this method if it needs
     * parameters that are written in the decay channel file.
     * Basically, this method gets called after reading the corresponding DC
     * file and it is its job to set the parameters in the current class to the
     * values written in the DC file.
     * 
     * @param _md Has to contain the values and names of all parameters
     * relevant for the current
     */
    virtual void SetModelParameters( struct GeneralModel _md ) {};

    //@{
    /// Getter/setter method
    inline std::string Name() const              { return m_name; }
    inline void        SetName(std::string name) { m_name = name; }
    inline std::string Path() const              { return m_path; }
    inline void        SetPath(std::string path) { m_path = path; }
    inline int         GetN() const              { return p_i.size(); }
    inline const std::vector<int>& DecayIndices() const { return p_i; }
    //@}

    /** 
     * Stream output operator giving the name of the current and its values.
     */
    friend std::ostream& operator<<(std::ostream& s, const Current_Base& cb);
  };
  

  typedef ATOOLS::Getter_Function<Current_Base,ME_Parameters>
  Current_Getter_Function;
}


#endif
