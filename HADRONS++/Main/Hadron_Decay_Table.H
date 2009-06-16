#ifndef HADRONS_Main_Hadron_Decay_Table_H
#define HADRONS_Main_Hadron_Decay_Table_H

#include "ATOOLS/Phys/Decay_Table.H"
#include "HADRONS++/Main/Tools.H"

namespace HADRONS {

class Hadron_Decay_Channel;

class Hadron_Decay_Table : public ATOOLS::Decay_Table {
private:
  bool ExtractFlavours(std::vector<int> &,std::string);
  void ExtractBRInfo( std::string entry, double & sbr, double & sdbr, std::string & sorigin );
public:
  Hadron_Decay_Table(ATOOLS::Flavour decayer);
  ~Hadron_Decay_Table();
  inline Hadron_Decay_Channel* Select() {
    return (Hadron_Decay_Channel*) Decay_Table::Select();
  }
  inline Hadron_Decay_Channel* at(const size_t i) {
    return (Hadron_Decay_Channel*) Decay_Table::at(i);
  }
  void Read(std::string path, std::string file);
  void Write(std::ostream& ostr);
  void Initialise(GeneralModel& startmd);
  void UpdateWidth(Hadron_Decay_Channel * hdc,const double &width);
}; // end of class

} // end of namespace


#endif
