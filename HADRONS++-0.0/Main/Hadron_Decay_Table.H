#ifndef Hadron_Decay_Table_H
#define Hadron_Decay_Table_H

#include "Decay_Table.H"
#include "Tools.H"

namespace HADRONS {

class Hadron_Decay_Channel;

class Hadron_Decay_Table : public ATOOLS::Decay_Table {
private:
  bool ExtractFlavours(std::vector<int> &,std::string);
  void ExtractBRInfo( std::string entry, double & sbr, double & sdbr, std::string & sorigin );
public:
  Hadron_Decay_Table(ATOOLS::Flavour decayer);
  ~Hadron_Decay_Table();
  void Read(std::string path, std::string file);
  void Write(std::ostream& ostr);
  void Initialise(GeneralModel& startmd);
  void UpdateWidth(Hadron_Decay_Channel * hdc,const double &width);
}; // end of class

} // end of namespace


#endif
