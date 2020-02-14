#ifndef EXTAMP_Main_Dipole_Wrapper_Process_H
#define EXTAMP_Main_Dipole_Wrapper_Process_H

#include "EXTAMP/CS_Dipole.H"
#include "EXTAMP/RS_Process.H"

#include "PHASIC++/Process/Single_Process.H"

namespace EXTAMP {

  /**
     A wrapper around a EXTAMP::CS_Dipole that basically adds zero
     functionality to it, apart from making it look like a
     PHASIC::Single_Process with a 'Partonic' method that returns the
     value of the dipole. This is necessary in the MC@NLO
     implementation due to AMEGIC's legacy structure of representing
     dipole terms by processes.
  */
  
  class Dipole_Wrapper_Process : public PHASIC::Single_Process {

  public:

    Dipole_Wrapper_Process(const RS_Process& rsproc,
			   EXTAMP::CS_Dipole* dipole,
			   BEAM::Beam_Spectra_Handler* beam,
			   PDF::ISR_Handler* isr);

    /* Calc dipole kinematics without evaluating dipole. Also perform
       trigger and set corresponding flag of NLO_subevt */
    void CalcKinematics(const ATOOLS::Vec4D_Vector& p);

    /* After having called CalcKinematics(), evaluate the dipole. */
    double Calc(ATOOLS::NLO_subevt* sub);

    void CalculateScale(const ATOOLS::Vec4D_Vector& real_p,
			const ATOOLS::Vec4D_Vector& born_p,
			ATOOLS::NLO_subevt* const evt);

    /* Inherited from PHASIC::Single_Process. Calls CalcKinematics and
       Calc, then returns the result. This is used in by Sherpa's
       MC@NLO implementation for calculating DA and DA-DS */
    double Partonic(const ATOOLS::Vec4D_Vector &p, int mode);

    /* A sign implementing DA and DA-DS as defined in arXiv:1111.1200 */
    int MCModeSign(ATOOLS::NLO_subevt* const evt) const;
    double GetMaxKT2ForDA() const;
    double GetKT2ofSplitting(const EXTAMP::Dipole_Kinematics& kin) const;

    bool Combinable(const size_t &idi,
		    const size_t &idj);

    void FillProcessMap(PHASIC::NLOTypeStringProcessMap_Map *apmap);

    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);

    /* Set properties of NLO_subevt according to this dipole */
    void SetSubEventProperties(ATOOLS::NLO_subevt& sub);

    /* Set a member pointer to the subevent passed */
    void AssignSubEvent(ATOOLS::NLO_subevt* sub);

    void SetScaleSetter(PHASIC::Scale_Setter_Base* scl) { p_scale = scl; }

    void SetNLOMC(PDF::NLOMC_Base *const nlomc);

    /* Set p_scale to NULL in order to avoid deleting it in destructor
       of PHASIC::Process_Base. It gets deleted in the RS_Process,
       which owns it. */
    ~Dipole_Wrapper_Process() { p_scale = NULL; }

  private:

    /* Flavours of this process have to be those of real emission
       config. NLO_subevts are assigned the born configuration,
       however. Store corresponding information in these members: */
    std::string m_born_name;
    PHASIC::Process_Info m_born_procinfo;
    ATOOLS::Flavour_Vector m_born_flavs;
    const ATOOLS::Flavour_Vector& BornFlavours() const { return m_born_flavs; }

    static PHASIC::Process_Info ConstructBornProcessInfo(const PHASIC::Process_Info& rsinfo,
							 size_t i, size_t j,
							 const ATOOLS::Flavour& flav_ij) ;

    /* This pointer (passed in constructor) is assumed to be valid
       throughout, associated memory is not managed. */
    EXTAMP::CS_Dipole* p_dipole;
    EXTAMP::CS_Dipole* Dipole() const { return p_dipole; }

    /* RS_Process is supposed to initialize and manages the
       subevent. */
    ATOOLS::NLO_subevt* p_subevent;

    /* Real emission momentum configuration with momenta of initial
       state particles reversed. This convention is used in
       NLO_subevts. */
    ATOOLS::Vec4D_Vector m_moms;
    const ATOOLS::Vec4D_Vector& Momenta() const { return m_moms; }

    /* Emitter, emitted, spectator index in real emission flavour
       vector. Follow Sherpa conventions: i<j */
    size_t I() const { return std::min(Dipole()->I(),Dipole()->J()); }
    size_t J() const { return std::max(Dipole()->I(),Dipole()->J()); }
    size_t K() const { return Dipole()->K(); }

    /* Position of combined flavour ij and spectator k in born
       configuration */
    size_t BornIJ() const { return m_inversemap[Dipole()->BornIJ()]; }
    size_t BornK()  const { return m_inversemap[Dipole()->BornK() ]; }

    /* Same as in EXTAMP::Process */
    double m_norm;
    const double& NormFac() const { return m_norm; }

    /* Mapping from ordering in p_dipole to ordering in this class
       (has to comply with Sherpa ordering conventions):
       Flavours[i] == Dipole()->Flavours[m_indexmap[i]] 
       Flavours[m_inversemap] == Dipole()->Flavours[i] */
    std::vector<size_t> m_indexmap;
    std::vector<size_t> m_inversemap;
    static std::vector<size_t> ConstructIndexMapping(const ATOOLS::Flavour_Vector& dipole_flavs,
						     const ATOOLS::Flavour_Vector& process_flavs,
						     size_t nin);
    static std::vector<size_t> InvertIndexMapping(const std::vector<size_t>& map);


    /* ID vector encoding the id's of particles as they are ordered in
       the born flavour vector. This is needed for
       ATOOLS::NLO_Subevents, which own a pointer to such a vector.
       id's are in binary encoding, meaning the i'th particle has id
       1<<i and the combined particles i and j have inded
       (1<<i|1<<j). */
    std::vector<size_t> m_id_vector;
    const std::vector<size_t>& IDVector() const {return m_id_vector; }
    std::vector<size_t> ConstructIDVector() const;

    std::map<size_t, ATOOLS::Flavour_Vector> m_cluster_flav_map;

  };

}

#endif
