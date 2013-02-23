#ifndef HADRONS_ME_Library_HD_ME_Base_H
#define HADRONS_ME_Library_HD_ME_Base_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include <string>
#include <utility>
#include "METOOLS/Main/Spin_Structure.H"
#include "HADRONS++/Main/Tools.H"
#include "ATOOLS/Org/Getter_Function.H"

#define DEFINE_ME_GETTER(CLASS,TAG)				       \
  DECLARE_GETTER(CLASS,TAG,HD_ME_Base,ME_Parameters);		       \
  HD_ME_Base* ATOOLS::Getter<HD_ME_Base,ME_Parameters,CLASS>::	       \
  operator()(const ME_Parameters &parameters) const		       \
  { return new CLASS(parameters.flavs, parameters.indices,TAG); }

namespace HADRONS {

  /** 
   * Represents an amplitude for a decay process.
   * To store its calculated values for each helicity combination, it inherits
   * from METOOLS::Spin_Amplitudes.
   * Each class derived from this has to implement the HD_ME_Base::Calculate
   * method for the actual amplitude calculation.
   */
  class HD_ME_Base : public METOOLS::Spin_Amplitudes {
  protected:
    /// Name of the amplitude, e.g. B_Bpi_pwave
    std::string       m_name;
    /// Reference to flavours involved in the decay (in external order)
    const ATOOLS::Flavour_Vector& m_flavs;
    /// Array of masses of the external particles
    double          * p_masses,* p_masses2;
    /// Path to the directory which contains the decay channel's file
    std::string       m_path;
    /// Index mapping from external flavours/momenta to amplitude internal ones
    std::vector<int>  p_i;
    /// Pre-factor for amplitude
    Complex           m_factor;

    virtual double lambdaNorm(const double,const double,const double);
  public:
    /** 
     * Constructor which initialises all fields of the class.
     * 
     * @param flavs Flavours involved in the decay (in external order)
     * @param decayindices Index mapping from external flavours/momenta to
     * amplitude internal ones (cf. HD_ME_Base::p_i)
     * @param name Name of the amplitude
     */
    HD_ME_Base(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& decayindices,
               const std::string& name);
    virtual ~HD_ME_Base();

    /** 
     * Pure virtual function, which requires the implementation of the
     * amplitude calculation in all derived classes.
     * When implementing such a calculation, one has to use HD_ME_Base::Insert
     * to store the calculated values for all helicities combination.
     * 
     * @param momenta Momenta in external order
     * @param anti Whether to consider the charge conjugated process
     */
    virtual void Calculate(const ATOOLS::Vec4D_Vector& momenta,bool anti)=0;

    /** 
     * Redefinition of the Spin_Amplitudes::Insert method, to take into
     * account the pre-factor HD_ME_Base::m_factor.
     * 
     * @param value Amplitude value (w/o pre-factor)
     * @param spins List of helicities for each particle in the decay,
     * describing the helicity combination for this <var>value</var>
     */
    inline void Insert(Complex value, std::vector<std::pair<int,int> >& spins)
    {
      (*this)[GetNumber(spins)]=m_factor*value;
    }

    /** 
     * This method sets the parameters for the decay amplitude.
     * Every subclass of HD_ME_Base \b must have this method if it needs
     * parameters that are written in the decay channel file.
     * Basically, this method gets called after reading the corresponding DC
     * file and it is its job to set the parameters in the ME class to the
     * values written in the DC file.
     * 
     * @param _md Has to contain the values and names of all parameters
     * relevant for the amplitude
     */
    virtual void   SetModelParameters( GeneralModel _md ) {};

    /** 
     * Tries to set the color flow for the particles given. This can be needed
     * for partonic decay channels, where the color flow is not trivial.
     * 
     * @param outparts Particles to set color flow for
     * @param n_q Number of gluons in outparts
     * @param n_g Number of quarks in outparts
     * @param anti Set the color flow for the charge conjugated decay.
     * 
     * @return true if setting the color flow succeeded, false if not.
     */
    virtual bool   SetColorFlow(std::vector<ATOOLS::Particle*> outparts,
                                int n_q, int n_g, bool anti);

    //@{
    /// Getter/setter method.
    void SetPath( std::string s ) { m_path=s; }
    virtual std::string Name() { return m_name; }
    inline void SetFactor(const Complex& factor) { m_factor=factor; }
    //@}
  };

  typedef ATOOLS::Getter_Function<HD_ME_Base,ME_Parameters>
  HD_ME_Getter_Function;
}


#endif
