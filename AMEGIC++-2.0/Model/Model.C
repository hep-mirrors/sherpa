#include "Model.H"
#include "Message.H"

using namespace AMEGIC;

void    Model::Init_Vertex() {
  AORGTOOLS::msg.Debugging()<<"Initialize new vertices !"<<std::endl;
  if (v!=NULL) delete v;
  v = new Vertex;
} 

void Model::Init() {
  AORGTOOLS::msg.Error()<<"Error: No Model::Init !"<<std::endl;
}

double Model::SinTW(){
  AORGTOOLS::msg.Out()<<"No SinTW included in this Model"<<std::endl;
  return 0.;
}
double Model::CosTW(){
  AORGTOOLS::msg.Out()<<"No CosTW included in this Model"<<std::endl;
  return 0.;
}
double Model::Aqed(){
  AORGTOOLS::msg.Out()<<"No Aqed included in this Model"<<std::endl;
  return 0.;
}
double Model::Aqed(double){
  AORGTOOLS::msg.Out()<<"No running Aqed included in this Model"<<std::endl;
  return 0.;
}
double Model::Aqcd(){
  AORGTOOLS::msg.Out()<<"No Aqcd included in this Model"<<std::endl;
  return 0.;
}
double Model::Aqcd(double){
  AORGTOOLS::msg.Out()<<"No running Aqcd included in this Model"<<std::endl;
  return 0.;
}
double Model::TanB(){
  AORGTOOLS::msg.Out()<<"No tan(beta) included in this Model"<<std::endl;
  return 0.;
}
double Model::CosA(){
  AORGTOOLS::msg.Out()<<"No cos(alpha) included in this Model"<<std::endl;
  return 0.;
}
double Model::SinA(){
  AORGTOOLS::msg.Out()<<"No sin(alpha) included in this Model"<<std::endl;
  return 0.;
}
Spectrum_EW*  Model::Get_Spectrum_EW(){
  AORGTOOLS::msg.Out()<<"No Spectrum_EW included in this Model"<<std::endl;
  return 0;
}
Spectrum_Higgs*  Model::Get_Spectrum_Higgs(){
  AORGTOOLS::msg.Out()<<"No Spectrum_Higgs included in this Model"<<std::endl;
  return 0;
}
Spectrum_sUpquarks*  Model::Get_Spectrum_sUpquarks(){
  AORGTOOLS::msg.Out()<<"No Spectrum_sUpquarks included in this Model"<<std::endl;
  return 0;
}
Spectrum_sDownquarks*  Model::Get_Spectrum_sDownquarks(){
  AORGTOOLS::msg.Out()<<"No Spectrum_sDownquarks included in this Model"<<std::endl;
  return 0;}
Spectrum_sLeptons*  Model::Get_Spectrum_sLeptons(){
  AORGTOOLS::msg.Out()<<"No Spectrum_sLeptons included in this Model"<<std::endl;
  return 0;
}
Spectrum_Neutralinos*  Model::Get_Spectrum_Neutralinos(){
  AORGTOOLS::msg.Out()<<"No Spectrum_Neutralinos included in this Model"<<std::endl;
  return 0;
}
Spectrum_Charginos*  Model::Get_Spectrum_Charginos(){
  AORGTOOLS::msg.Out()<<"No Spectrum_Charginos included in this Model"<<std::endl;
  return 0;
}

