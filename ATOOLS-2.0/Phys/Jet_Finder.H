#ifndef Jet_Finder_h
#define Jet_Finder_h

#include "Particle_List.H"
#include "Selector.H"
#include "Poincare.H"

namespace ATOOLS {
  class Jet_Finder : public Selector_Base {
    double m_ycut,m_delta_r;
    int    m_type;
    bool   m_pt_def;
    double m_ene, m_s, m_sprime; //, m_smin, m_smax;

    double m_shower_pt2,m_value;

    std::vector<std::set<int> > m_combs;

    /*
      The vectors related to the frame
    */
    Vec4D *  p_frame;
    Poincare m_cms_boost;
    Poincare m_zrot;

    void FillCombinations(const std::string &name,int &cp);
    void FillCombinations();

    void BoostInFrame(Vec4D *);
    void BoostBack(Vec4D *);
    void BoostInFrame(std::vector<Vec4D>);
    void BoostBack(std::vector<Vec4D>);
    void BoostInFrame(Vec4D &);
    void BoostBack(Vec4D &);

    bool ConstructJetSystem(Vec4D *,Flavour *,std::vector<Vec4D> &,
			    std::vector<int> jets,std::vector<double> & lastys);

    double YminKt(Vec4D *,int &,int &);
    double YminKt(Vec4D *,Flavour *,std::vector<Vec4D>,int &,int &);
    double Yij(Vec4D,Vec4D);

    double DEta12(Vec4D &,Vec4D &);
    double CosDPhi12(Vec4D &,Vec4D &);
    double DCos12(Vec4D &,Vec4D &);
    double Coshyp(double);
    bool   ConstructJets(Particle_List *, double, int, bool );
  public:
    Jet_Finder(const double,const int=1,const bool=true);
    Jet_Finder(const int,Flavour *,const double,const int,const bool=true);

    ~Jet_Finder();

    void   Init(const Vec4D *);
    bool   ConstructJets(const Particle_List *,
			 const std::vector<int> &,std::vector<double> & ,bool);
    bool   ConstructJets(Particle_List *, double, bool );
    bool   ConstructJets(Particle_List *, int, bool );


    bool   Trigger(const Vec4D *);

    bool   TwoJets(const Vec4D &,const bool=false);
    bool   TwoJets(const Vec4D &, const Vec4D &,const bool=false);
    bool   TwoJets(double &,double &,double &,bool);
    double MTij2(Vec4D,Vec4D);
    void   BuildCuts(Cut_Data *);
    void   UpdateCuts(double,double,Cut_Data *);

    double   ActualValue() const;

    inline void SetType(const int mode) { m_type=mode;   }
    inline int  Type() const            { return m_type; }

    int      IsConditional()                const { return m_type==1; }
    int      NeedUpdate()                   const { return m_type==1; }
    double   ShowerPt2()                    const { return m_shower_pt2; }
    double   Smin()                         const { return m_smin; }
    double   Smax()                         const { return m_smax; }
    double   DeltaR()                       const { return m_delta_r; }
    double   Ycut()                         const { return m_ycut; }
    void     SetShowerPt2(const double pt2)       { m_shower_pt2=pt2; }
    void     SetDeltaR(double dr);
  };

  inline bool Jet_Finder::ConstructJets(Particle_List * pl, double y_cut, bool final_only) { 
    return ConstructJets(pl,y_cut,0,final_only); 
  }
  inline bool Jet_Finder::ConstructJets(Particle_List * pl, int number, bool final_only) { 
    return ConstructJets(pl,1.0,number,final_only); 
  }






/*! 
  \class Jet_Finder
  \brief A general jet finder for Durham (k_t)-jets.

  The jet definition is managed by the Jet_Finder 
  giving a flag for whether it's a Hadron-Hadron(4), 
  Photon-Hadron(3), Lepton-Hadron(2), or a Lepton-Lepton(1) 
  collision, or whether it is an ordinary decay (0).
  For the usage in ME'S and in the interface to the
  parton shower, the number of particles and their flavours
  are also important.
  For the k_t algorithm the following frames have to
  be employed :
  
\verbatim
  Collision type           | Frame
  hadron-hadron            | c.m. of the hadrons
  lepton-lepton            | c.m. frame of the leptons
  DIS (electro-production) | Breit-frame
  DIS (photon-production)  | c.m. of photon-hadron
\endverbatim  

  So far no thought was lost on photon-photon collisions, since 
  they might be of basically any type. Presumably a new type has to
  be set up with an awful lot of conditions in it due to the
  fully, partly or un-resolved nature of this collision type.
  Therefore, before applying the k_t algorithm, a boost has 
  to be performed to the relevant frame.
*/

/*!
  \var   int Jet_Finder::m_type
  \brief The general type of process - it affects the proper jet definition.
  
  Decay = 0, 
  Lepton-Lepton = 1, Lepton-Hadron = 2, 
  Photon-Hadron = 3, Hadron-Hadron = 4

  Note that choices 1 and 2 are very similar, the appropriate reference frame
  is the rest-frame of the intermediate boson, whereas choices 3 and 4 again are very 
  similar - 3 looks into the hadronic structure of a photon.

  \todo I guess we'll have to work this out in more detail.
*/

/*!
  \var   double Jet_Finder::m_ycut
  \brief The jet resolution parameter.
*/

/*
  \var   double E, s, sprime, smin, smax;
  \brief General scales : energy, s, the actual s, its limits.
*/

/*!
  \var   Vec4D * Jet_Finder::m_frame;
  \brief The vectors related to the frame.

  These vectors span the corresponding frame in which the jet-algorithm is
  performed. It is the c.m. frame of the vector boson (i.e. the vectors of the
  leptons) or the c.m. frame of the hadrons.
*/

/*!
  \fn    double Jet_Finder::YminKt(Vec4D *,int &,int &);
  \brief Find the minimal y for a number of vectors, fill in the appropriate 
         numbers.

  This method scans over all momenta (being constrained by n) and
  determines the minimal y in the k_t scheme. This should be
  used for matrix elements only.
*/

/*!
  \fn    Jet_Finder::Jet_Finder(double,int);
  \brief The explicit constructor for the usage in the Parton shower.

  This explicit constructor just sets the ycut of the jet criterion,
  and it has to be stgressed again, this is a "quick'n'dirty" jet criterion
  that works SOLELY within the Parton shower. 
*/

/*!
  \fn    Jet_Finder::Jet_Finder(int,Flavour *,double,int,int);
  \brief The explicit constructor for the usage in matrix elements and their 
         interface to the Parton shower.

  This explicit constructor sets the flavours of the process
  (here, care has to be taken to extract "good" flavours, i.e.
  the beam flavour), and sets up a default frame.
  For hadron-hadron collisions and for decays the boosts are
  already set up, because it is implicitly assumed that the 
  energy of the c.m. frame is fixed. In other words, phenomena
  like the smearing of hadron momenta due to collider effects
  or of decaying particles due to finite width effects are 
  discarded.
*/

/*!
  \fn    void Jet_Finder::Init(Vec4D *);
  \brief Setting up the reference system for a number of vectors.

  In this method, only the initial momenta are analysed and
  boosts into the frame of the jet_finder and back are set up.
  This method is important for : 
  - lepton-lepton collisions due to varying c.m. systems of the 
    lepton pair as a consequence of ISR effects.
  - photoproduction and electroproduction because there 
    varying photon momenta are unavoidable.
  - photon-photon collisions.

  \todo  The last two cases have not been implemented yet.
*/

/*!
  \var   bool Jet_Finder::Trigger(Vec4D *);
  \brief Checks for pairs of vectors to belong to one jet only.

  This method yields a zero if any pair of momenta of strong 
  interacting particles in the set of momenta is not
  separated according to the jet-measure.
*/

/*!
  \var   bool Jet_Finder::TwoJets(Vec4D &,Vec4D &);
  \brief Checks, if the two vectors belong to one jet.

  This method returns a zero if the two vectors do not fulfil the 
  jet-criterion. This is to be used basically everywhere,
  in matrix elements and the parton shower.
*/

/*!
  \var   bool Jet_Finder::TwoJets(double &, double &, double &, bool)
  \brief Checks, if the two vectors belong to one jet.

  This method returns a zero if the two vectors do not fulfil the 
  jet-criterion. This is to be used in the parton shower only,
  and the flag decides whether two final state (0) or one inital and
  one final state particle (1) is under consideration.
*/

/*!
  \var   double Jet_Finder::PTij(Vec4D,Vec4D);
  \brief Returns the actual p_t of two vectors.

  Depending on the set-up of the jet-finder the actual value of
  p_t of two vectors is calculated and returned.
*/

/*!
  \var   double Jet_Finder::PTi(Vec4D);
  \brief Returns the actual p_t of a vector w.r.t. the beam.

  Depending on the set-up of the jet-finder the actual value of
  p_t of the vectors is calculated and returned.
*/

/*!
  \var   void Jet_Finder::BuildCuts(Cut_Data *);
  \brief Creates cuts to be used in the evaluation of the ME.
*/

}

#endif








