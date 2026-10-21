#ifndef METOOLS_Current_Library_Current_Base_H
#define METOOLS_Current_Library_Current_Base_H

#include <iostream>
#include <string>
#include <iomanip>
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "METOOLS/HadronCurrents/Photon_Emission.H"
#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"


#define DEFINE_CURRENT_GETTER(CLASS,TAG)				\
  DECLARE_GETTER(CLASS,TAG,METOOLS::Current_Base,METOOLS::ME_Parameters);		\
  METOOLS::Current_Base *ATOOLS::Getter<METOOLS::Current_Base,METOOLS::ME_Parameters,CLASS>:: \
  operator()(const ME_Parameters &parameters) const		\
  { return new CLASS(parameters.flavs, parameters.indices, TAG); }

namespace METOOLS {

  /** 
   * Represents a vector current in a decay process.
   * To store its calculated values for each helicity combination, it inherits
   * from METOOLS::Spin_Structure.
   * Each class derived from this has to implement the Current_Base::Calc
   * method for the actual calculation of the current.
   */
  class Current_Base : public Spin_Structure<ATOOLS::Vec4C> {
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

    /// How many of the EXTERNAL legs are incoming: 1 for a decay, where
    /// leg 0 is the decaying particle, 2 for a scattering process. A leg
    /// index i is incoming iff i < m_nin.
    ///
    /// Defaults to the decay convention, so a current built by the legacy
    /// Current_ME path behaves exactly as it always has. Exchange_ME sets
    /// it explicitly from the context it was built in.
    size_t            m_nin;

    /// Conventions and coupling for a real photon emitted from inside this
    /// current. Unused unless the current is handed a photon leg; set by
    /// Exchange_ME. See Photon_Emission.H.
    Photon_Emission   m_photon;

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
     * @param moms All momenta in the decay in external order
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

    /**
     * Check that the tag this current was built under matches the topology
     * it was actually handed. Called by SetNIn, i.e. once the current knows
     * which of its legs are incoming -- which it cannot know at
     * construction. The default does nothing; currents with a meaningful
     * in/out distinction override it.
     */
    virtual void CheckTopology() const {}

    /// True if the current is identically zero for its chosen model. Lets
    /// Exchange_ME stop a process whose every exchange vanishes, while a
    /// single vanishing exchange (e.g. the photon on a pointlike neutron
    /// next to a Z) is allowed.
    virtual bool Vanishes() const { return false; }

    /// Tell the current how many external legs are incoming, then validate.
    inline void SetNIn(const size_t nin) { m_nin = nin; CheckTopology(); }

    inline Photon_Emission & Photon() { return m_photon; }

    /// Whether external leg index i is an incoming leg.
    inline bool IsIncoming(const int i) const { return i<int(m_nin); }

    /// How many of THIS current's legs are incoming.
    inline size_t NIncomingLegs() const
    {
      size_t n(0);
      for (size_t k(0);k<p_i.size();++k) if (IsIncoming(p_i[k])) ++n;
      return n;
    }

    //@{
    /// Getter/setter method
    inline std::string Name() const              { return m_name; }
    inline void        SetName(std::string name) { m_name = name; }
    inline int         GetN() const              { return p_i.size(); }
    inline const std::vector<int>& DecayIndices() const { return p_i; }
    inline const ATOOLS::Flavour_Vector& Flavs()  const { return m_flavs; }
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
