#ifndef HADRONS_Main_Tools_H
#define HADRONS_Main_Tools_H

#ifndef SQRT_05
#define SQRT_05 0.70710678118654757
#endif

#include <map>
#include <iostream>
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Matrix.H"
#include "PHASIC++/Decays/Decay_Table.H"

namespace HADRONS
{

  /** 
   * Object containing all parameters needed to initialise a HD_ME_Base
   * or Current_Base. This is necessary as argument for the Getter_Function
   * of these objects.
   */
  struct ME_Parameters {
    const ATOOLS::Flavour_Vector& flavs;
    const std::vector<int>& indices;
    ME_Parameters(const ATOOLS::Flavour_Vector& _flavs,
                  const std::vector<int>& _indices) :
      flavs(_flavs), indices(_indices)
    {}
  };

  /**
   * A map with all relevant parameters for a certain decay amplitude.
   * It is a map that both contains the name of the parameter (as it is written
   * in the decay channel file) and its (double precision) value.
   */
  struct GeneralModel: public std::map<std::string,double> {
    std::map<kf_code, kf_code> m_aliases;
    
  public:
    /** 
     * Returns the value of the parameter with the name <var>tag</var>.
     * If such a parameter does not exist, it returns the default value given
     * by <var>def</var>.
     */
    inline double operator()(const std::string &tag,const double &def) const 
    {
      std::map<std::string,double>::const_iterator fit(find(tag));
      return fit!=end()?fit->second:def;
    }
  };

  /** 
   * Utility class providing some regularly needed tools. They are all
   * implemented as static objects, such that they can be used without
   * instantiation.
   */
  class Tools {
  public:
    static PHASIC::Decay_Table * partonic_b, * partonic_c;
    /** 
     * Phase space function \f$\lambda = (a-b-c)^2-4bc\f$;
     * 
     * @param a 
     * @param b 
     * @param c 
     * 
     * @return \f$ \lambda \f$
     */

    static double   Lambda(double a,double b,double c);
    static Complex  BreitWigner( double s, double Mass2, double MassWidth );
    static Complex  BreitWignerFix( double s, double Mass2, double MassWidth );
    
    static double   OffShellMassWidth( double s, double Mass2, 
				       double Width, double ms );
    static double   OffShellMassWidth( double s, double Mass2, 
				       double Width, double ms1, double ms2 );

    /** 
     * Extracts the kfcodes from a string like {-211,111,22}
     * 
     * @param kfc The returned kfcodes will be inserted here.
     * @param tag The string to be parsed.
     * 
     * @return bool indicating whether parsing succeeded.
     */
    static bool     ExtractFlavours(std::vector<int>& kfc,std::string tag);
    static bool     ExtractSpecs(std::string entry,std::vector<int> & specs,
				 std::vector<double> & specweights);
    /** 
     * Extracts branching ratio information from a string like 0.5(0.01)[PDG]
     * 
     * @param entry The string to be parsed.
     * @param sbr Extracted branching ratio (e.g. 0.5)
     * @param sdbr Extracted error on branching ratio (e.g. 0.01)
     * @param sorigin Extracted source of the branching ratio (e.g. PDG)
     */
    static void     ExtractBRInfo(std::string entry, double & sbr, 
				  double & sdbr,std::string & sorigin );

    //@{
    /** CKM matrix element. */
    static const double Vud;
    static const double Vus;
    static const double Vub;
    static const double Vcd;
    static const double Vcs;
    static const double Vcb;
    static const double Vtd;
    static const double Vts;
    static const double Vtb;
    //@}
  }; 
}

#endif
