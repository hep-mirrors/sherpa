#include "Primitive_Analysis.H"


using namespace APHYTOOLS;

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

void Primitive_Analysis::DoAnalysis(double value) {
  for (int i=0;i<observables.size();i++) {
    observables[i]->Evaluate(value);
  }
}

void Primitive_Analysis::DoAnalysis(const APHYTOOLS::Parton_List & pl, double value) {
  for (int i=0;i<observables.size();i++) {
    observables[i]->Evaluate(pl,value);
  }
}

void Primitive_Analysis::DoAnalysis(const APHYTOOLS::Blob_List & bl, double value) {
  for (int i=0;i<observables.size();i++) {
    observables[i]->Evaluate(bl,value);
  }
}

void Primitive_Analysis::FinishAnalysis(std::string resdir,int tables) {
  int  mode_dir = 448;
  mkdir(resdir.c_str(),mode_dir); 
  for (int i=0;i<observables.size();i++) {
    observables[i]->EndEvaluation();
    observables[i]->Output((resdir).c_str());
  }
}





