#ifndef Hadrons_Decays_H
#define Hadrons_Decays_H

#include "Decay_Table_Reader.H"
#include "Decay_Table.H"
#include <string>

namespace HADRONS {
  class Hadron_Decays {
  private:
    std::string m_path, m_file;

    void ReadInDecayTables();
    ATOOLS::Decay_Table * InitialiseOneDecayTable(std::vector<std::string>);
  public: 
    Hadron_Decays(std::string _path=std::string("./"),
		  std::string _file=std::string("HadronDecays.dat"));
  };
}


#endif
