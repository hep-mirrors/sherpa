#ifndef Decay_Channel_h
#define Decay_Channel_h

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle.H"
#include <set>

namespace METOOLS {
  class Spin_Amplitudes;
  class Amplitude2_Tensor;
  class Spin_Density;
}

namespace PHASIC {
  class Multi_Channel;
  class Single_Channel;

  class Decay_Channel {
  private:
    static std::map<std::string, double> s_kinmaxfails;
    
    /** 
     * This method makes the phase space generator choose a phase space point
     * and calculate the corresponding ME and PS weights.
     * The ME weights (\b not squared) for each helicity saved in
     * Decay_Channel::m_diagrams are then summed for all available diagrams and
     * helicities taking into account the symmetry factor and
     * averaging over the incoming helicity states.
     * So far, this Differential method does not care about spin correlations
     * it merely sums up all absolute squared entries of those
     * in order to obtain a value for \f$d\Gamma\f$, and averages over incoming
     * spin.
     *
     * @param momenta Reference to phase space point, which has size=NOut()+1
     * and the incoming momentum already set.
     * @param anti Whether to consider the charge conjugated matrix elements.
     * @param sigma If not NULL, this spin density matrix is considered for the
     * incoming particle.
     * 
     * @return \f$d\Gamma(\f$
     */
    double Differential(ATOOLS::Vec4D_Vector& momenta, bool anti,
                        METOOLS::Spin_Density* sigma=NULL,
                        const std::vector<ATOOLS::Particle*>& p=
                          std::vector<ATOOLS::Particle*>());
    double ME2(const ATOOLS::Vec4D_Vector& momenta, bool anti,
               METOOLS::Spin_Density* sigma=NULL,
               const std::vector<ATOOLS::Particle*>& p=
               std::vector<ATOOLS::Particle*>());

  protected :
    double      m_width, m_deltawidth, m_minmass, m_max, m_symfac;

    /// integrated width (as opposed to the one given by BR)
    double m_iwidth;
    /// integrated deltawidth (as opposed to the one given by BR)
    double m_ideltawidth;

    /**
     * m_active[n] specifies the status of the decay channel for the decay
     * of the <n>'th particle of this type in the event.
     * It can take the following values: m_active[n]=
     *   0: Disabled decay channel, but contributes to total width
     *   1: Enabled decay channel
     *   2: Forced decay channel
     *   -1: Disabled decay channel, and does not contribute to width
     */
    std::vector<int> m_active;

    ATOOLS::Flavour_Vector m_flavours;
    std::vector<METOOLS::Spin_Amplitudes*> m_diagrams;
    PHASIC::Multi_Channel* p_channels;
    METOOLS::Amplitude2_Tensor* p_amps;
    /**
     * The Mass_Selector knows which type of mass is to be used for all decays
     * associated with this decay table. This can be the HadMass() if the decays
     * belong to the non-perturbative regime of the event, or the Mass() if
     * they belong to the perturbative part. 
     */
    const ATOOLS::Mass_Selector* p_ms;
    
    double Lambda(const double& a, const double& b, const double& c) const;
    double MassWeight(const double& s, const double& sp,
                      const double& b, const double& c) const;

  public :
    Decay_Channel(const ATOOLS::Flavour &, const ATOOLS::Mass_Selector* ms);
    virtual ~Decay_Channel();
    void AddDecayProduct(const ATOOLS::Flavour & flout,const bool & sort=true);
    inline void SetWidth(const double & width) { m_width =width; }
    inline void SetDeltaWidth(const double & delta) { m_deltawidth=delta; }

    std::string Name() const;
    std::string IDCode() const;
    static std::string IDCode(const ATOOLS::Flavour& decayer,
                              const ATOOLS::Flavour_Vector& daughters);
    inline const double& Width() const { return m_width; }
    inline const double& DeltaWidth() const { return m_deltawidth; }
    inline const double& IWidth() const { return m_iwidth; }
    inline void SetIWidth(const double& width) { m_iwidth=width; }
    inline const double& IDeltaWidth() const { return m_ideltawidth; }
    inline void SetIDeltaWidth(const double& delta) { m_ideltawidth=delta; }
    inline const double& Max() const { return m_max; }
    inline void SetMax(const double& max) { m_max=max; }
    inline const double& MinimalMass() const { return m_minmass; }
    inline int Active(const size_t& counter) const { return counter<m_active.size() ? m_active[counter] : m_active[m_active.size()-1]; }
    inline void SetActive(std::vector<int> active) { m_active=active; }
    inline void SetActive(size_t counter, int active) { m_active[counter]=active; }
    inline void SetActiveAll(int active) { for (size_t i=0; i<m_active.size(); ++i) m_active[i]=active; }

    double GenerateMass(const double& max, const double &width) const;
    double SymmetryFactor();

    inline int NOut() const { return m_flavours.size()-1; }
    inline const std::vector<ATOOLS::Flavour>& Flavs() const{return m_flavours;}
    inline const ATOOLS::Flavour& GetDecaying() const { return m_flavours[0]; }
    inline const ATOOLS::Flavour& GetDecayProduct(size_t i) const {
      return m_flavours[i+1];
    }
    inline const std::vector<METOOLS::Spin_Amplitudes*>& GetDiagrams() const
    { return m_diagrams; }
    void AddDiagram(METOOLS::Spin_Amplitudes* amp);
    inline PHASIC::Multi_Channel* Channels() const { return p_channels; }
    inline void SetChannels(PHASIC::Multi_Channel* chans) { p_channels=chans; }
    void AddChannel(PHASIC::Single_Channel* chan);

    inline void ResetDiagrams() { m_diagrams.clear(); }
    void ResetChannels();

    void Output() const;

    void CalculateWidth(double acc=0.01, double ref=0.0, int iter=2500);
    void GenerateKinematics(ATOOLS::Vec4D_Vector& momenta, bool anti,
			    METOOLS::Spin_Density* sigma=NULL,
			    const std::vector<ATOOLS::Particle*>& p=
			    std::vector<ATOOLS::Particle*>());
    inline METOOLS::Amplitude2_Tensor* Amps() { return p_amps; }

    friend std::ostream &operator<<(std::ostream &os, const Decay_Channel &dt);

    static void SortFlavours(ATOOLS::Flavour_Vector& flavs);
    static bool FlavourSort(const ATOOLS::Flavour &fl1,
                            const ATOOLS::Flavour &fl2);
    static void PrintMaxKinFailStatistics(std::ostream &str);
  };
}
#endif
