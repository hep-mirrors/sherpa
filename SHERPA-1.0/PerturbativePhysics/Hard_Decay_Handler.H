#ifndef Hard_Decay_Handler_h
#define Hard_Decay_Handler_h

#include "Amegic.H"
#include "Process_Base.H"
#include "Model_Base.H"
#include "Decay_Table.H"

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

  class Hard_Decay_Handler {
  private :
    std::string m_path,m_file;

    DecayTables m_decaytables;
    AMEGIC::Amegic * p_amegic;

    void ReadInDecays();
    void EvaluateWidths(std::string,MODEL::Model_Base *);
  public :
    Hard_Decay_Handler(std::string,std::string,std::string,MODEL::Model_Base *);
    ~Hard_Decay_Handler();
    std::string Name();
    AMEGIC::Amegic * GetAmegic() { return p_amegic; }
  };
}

#endif

