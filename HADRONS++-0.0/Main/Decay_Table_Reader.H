#ifndef Decay_Table_Reader_H
#define Decay_Table_Reader_H

#include "Decay_Table.H"
#include <string>
#include <vector>

namespace HADRONS {
  class Decay_Table_Reader {
  private:
    std::vector<std::vector<std::string> > m_helpsvv;
    bool ExtractFlavours(std::vector<int> &,std::string);
  public: 
    Decay_Table_Reader(std::string path=std::string("./"),
		       std::string file=std::string(""));
    int  FillDecayTable(ATOOLS::Decay_Table * dt=NULL);
    void FillInMatrixElementsAndPS(ATOOLS::Decay_Table * dt=NULL);
  };
}


#endif
