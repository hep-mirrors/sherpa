#include "Photon_Definition.H"
#include "Primitive_Detector.H"
#include "Primitive_ElMag_Calorimeter.H"
#include "MyStrStream.H"
#include "Particle_Qualifier.H"
#include <iomanip>

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Photon_Definition_Getter,"Photon",
 	       Object_Definition_Base,Argument_Matrix);

void Photon_Definition_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"blabla\n"
     <<std::setw(width+4)<<" "<<"}";
}


Object_Definition_Base * 
Photon_Definition_Getter::operator()(const Argument_Matrix &parameters) const {
  return NULL;
}


Photon_Definition::Photon_Definition(const std::string order) :
  Object_Definition_Base("photons",kf::photon,order) 
{}

void Photon_Definition::FillMCTruthList(Particle_List * finalstate,
					Primitive_Detector * det) {
  p_data->ResetPList();
  for (std::deque<Particle*>::iterator pit=finalstate->begin();
       pit!=finalstate->end();) {
    if ((*pit)->Flav().Kfcode()==m_code) {
      p_data->AddPToPList((*pit));
      pit = finalstate->erase(pit);
    }
    else ++pit;
  }
}

void Photon_Definition::FillSimpleDetectorList(Particle_List * finalstate,
					       Primitive_Detector * det) {
  Primitive_ElMag_Calorimeter * ecal = 
    dynamic_cast<Primitive_ElMag_Calorimeter *>(det->GetElement("ECal"));
  if (ecal==NULL) { abort(); }
  long int neta(-1),nphi(-1);
  ecal->GetNumbersOfCells(neta,nphi);
  for (int etapos=0;etapos<neta;etapos++) {
    for (int phipos=0;phipos<nphi;phipos++) {
      if (ecal->GetFlav(etapos,phipos)==22) {
	p_data->AddPToPList(new Particle(0,Flavour(kf::photon),
					 ecal->ReconstructMasslessFourMom(etapos,phipos),'X'));
      }
    }
  }
  p_data->SortPList();
}
