#include "Primitive_Analysis.H"


using namespace ATOOLS;

Primitive_Analysis::Primitive_Analysis(std::string _name)
{
  name = std::string("Analysis : ") + _name;
}

Primitive_Analysis::Primitive_Analysis() 
{
  name = std::string("noname");
}

Primitive_Analysis::~Primitive_Analysis()
{
  for(int i=observables.size();i>0;i--) {
    if (observables[i-1]) delete observables[i-1];
  }
}

void Primitive_Analysis::AddObservable(Primitive_Observable_Base * obs) {
  observables.push_back(obs);
}


void Primitive_Analysis::SetBlobType(std::string _btype) {
  for (int i=0;i<observables.size();i++) observables[i]->SetBlobType(_btype);
}


void Primitive_Analysis::DoAnalysis(double _weight) {
  for (int i=0;i<observables.size();i++) observables[i]->Evaluate(_weight);
}

void Primitive_Analysis::DoAnalysis(const Particle_List & pl, double _weight) {
  for (int i=0;i<observables.size();i++) observables[i]->Evaluate(pl,_weight);
}

void Primitive_Analysis::DoAnalysis(const Blob_List & bl, double _weight) {
  for (int i=0;i<observables.size();i++) observables[i]->Evaluate(bl,_weight);
}

void Primitive_Analysis::FinishAnalysis(std::string resdir,int tables) {
  std::cout<<"In Primitive_Analysis::FinishAnalysis : "<<resdir<<std::endl;
  int  mode_dir = 448;
  mkdir(resdir.c_str(),mode_dir); 
  for (int i=0;i<observables.size();i++) {
    observables[i]->EndEvaluation();
    observables[i]->Output((resdir).c_str());
  }
}





