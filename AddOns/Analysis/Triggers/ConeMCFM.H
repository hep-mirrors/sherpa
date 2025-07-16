/* Feb 6th 2009
 * 
 * Close copy of MCFM-5.3/src/User/genclust_cone.f cone algorithm implementation.
 * Implementation for up to 4 jets.
 * Uses sherpa infrastructure.
 */

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Math/Poincare.H"

namespace ATOOLS {
  class ConeMCFM {
    double m_Rmin, m_Rsep;

    void MergeSplit(std::vector<Vec4D> & qjet ,std::vector<Vec4D> & protoq, std::vector< std::vector<int> > & protoc, std::vector<Vec4D> & qfinal);

    double deltarq(Vec4D p1,Vec4D p2);
    double etarap(Vec4D p1);
    double getet(Vec4D p1);

  public:
    std::vector<Vec4D> m_pjets; 						/* jet momenta */
    void ConstructJets(std::vector<Vec4D> & p); 				/* constructs jet momenta stored in m_pjets.  */	
    ConeMCFM(double Rmin, double Rsep);     					/* constructs on jet momenta stored.         */
    ~ConeMCFM() {}

  };
}
