#ifndef Hard_Decay_Handler_h
#define Hard_Decay_Handler_h

#include "Matrix_Element_Handler.H"
#include "Amegic.H"
#include "Process_Base.H"
#include "Model_Base.H"
#include "Decay_Table.H"
#include "Kinematics.H"
#include <map>

namespace SHERPA {
  class Decay_Table_Sorting_Criterion {
  public :
    bool operator() (ATOOLS::Decay_Table * dt1,ATOOLS::Decay_Table * dt2) {
      if (dt1->Flav().Mass()>dt2->Flav().Mass()) return 0;
      return 1;
    }
  };

  typedef std::set<ATOOLS::Decay_Table *,Decay_Table_Sorting_Criterion> DecayTables;
  typedef DecayTables::iterator DecIt;


  class Decaying_Particle_Table_Sorting_Criterion {
  public :
    bool operator() (ATOOLS::Particle * part1,ATOOLS::Particle * part2) {
      if (part1->Flav().Width()>part2->Flav().Width()) return 0;
      return 1;
    }
  };

  typedef std::map<ATOOLS::Particle *,ATOOLS::Decay_Channel *,
		   Decaying_Particle_Table_Sorting_Criterion> DecayingParticleTable;
  typedef DecayingParticleTable::iterator DPTIt;


  class Hard_Decay_Handler {
  private :
    bool                     m_on;
    std::string              m_path,m_file;
    DecayTables              m_decaytables;
    DecayingParticleTable    m_table;
    Kinematics             * p_tools;
    Matrix_Element_Handler * p_mehandler;
    bool                     m_meflag;
    AMEGIC::Amegic         * p_amegic;
    bool                     m_amegicflag;

    void                    ReadInDecays();
    void                    EvaluateWidths(std::string,MODEL::Model_Base *);
    void                    SetWidths();
    ATOOLS::Decay_Channel * SpecifyHardDecay(ATOOLS::Particle *,double &);  
  public :
    Hard_Decay_Handler(std::string,std::string,std::string,MODEL::Model_Base *);
    ~Hard_Decay_Handler();

    double DefineSecondaryDecays(ATOOLS::Blob *,bool = true);
    bool   PerformDecay(ATOOLS::Blob *);

    void                     ResetTables();
    std::string              Name();
    Matrix_Element_Handler * GetMEHandler() { return p_mehandler; }
    AMEGIC::Amegic         * GetAmegic()    { return p_amegic; }
    bool                     On()           { return m_on; }
  };
}

#endif

