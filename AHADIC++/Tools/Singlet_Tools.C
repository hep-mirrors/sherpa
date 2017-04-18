#include "AHADIC++/Tools/Singlet_Tools.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Singlet::Singlet() {
  //msg_Out()<<METHOD<<"("<<this<<")\n";
}

Singlet::~Singlet() {
  //msg_Out()<<METHOD<<"("<<this<<", size = "<<size()<<").\n";
  while (size()>0) {
    delete front();
    pop_front();
  }
  clear();
}

bool Singlet::ReorderCriterion(Proto_Particle * first) {
  // checks if first particle is given particle (which is assumed to be
  // one at the beginning of the singlet, or, if no particle is given, if
  // the first particle is a gluon, anti-quark or diquark.
  if (first) return (first!=front());
  return (front()->Flavour().IsGluon() ||
	  (front()->Flavour().IsQuark() && front()->Flavour().IsAnti()) ||
	  (front()->Flavour().IsDiQuark() && !front()->Flavour().IsAnti()));
}

void Singlet::Reorder(Proto_Particle * first) {
  // rotating the singlet until a quark or anti-diquark is fronting it
  //msg_Out()<<METHOD<<"("<<first<<", "<<front()<<").\n";
  while (ReorderCriterion(first)) {
    push_back(front());
    pop_front();
  }
}

bool Singlet::Combine(Proto_Particle * part1,Proto_Particle * part2) {
  // quark-gluon pair not heavy enough, just add them.
  if (!part1->Flavour().IsGluon() && part2->Flavour().IsGluon()) {
    part1->SetMomentum(part1->Momentum()+part2->Momentum());
    Erase(part2);
    return true;
  }
  // gluon-quark pair not heavy enough, just add them.
  else if (part1->Flavour().IsGluon() && !part2->Flavour().IsGluon()) {
    part2->SetMomentum(part1->Momentum()+part2->Momentum());
    Erase(part1);
    return true;
  }
  // gluons only - add them and kill one
  else if (part1->Flavour().IsGluon() && part2->Flavour().IsGluon()) {
    part2->SetMomentum(part1->Momentum()+part2->Momentum());
    Erase(part1);
    return true;
  }
  return false;
}

void Singlet::Erase(Proto_Particle * ref) {
  for (list<Proto_Particle *>::iterator pit=begin();
       pit!=end();pit++) {
    if ((*pit)==ref) {
      delete ref;
      erase(pit);
      break;
    }
  }
}

double Singlet::Mass2() const {
  return Momentum().Abs2();
}

Vec4D Singlet::Momentum() const {
  Vec4D mom(0.,0.,0.,0.);
  for (list<Proto_Particle *>::const_iterator pliter=begin();
       pliter!=end();pliter++) mom += (*pliter)->Momentum();
  return mom;
}

std::ostream& AHADIC::operator<<(std::ostream & str,const Singlet & sing) {
  str<<"****** Singlet ("<<sing.size()<<", mass = "<<sqrt(sing.Mass2())
     <<") **********:\n";
  for (list<Proto_Particle *>::const_iterator pliter=sing.begin();
       pliter!=sing.end();pliter++)
    str<<"*** "<<(*pliter)<<" : "
       <<(*pliter)->Flavour()
       <<" (lead = "<<(*pliter)->IsLeading()<<", "
       <<"beam = "<<(*pliter)->IsBeam()<<")"
       <<": "<<(*pliter)->Momentum()<<"\n";
  str<<"***********************************************************\n";
  return str;
}
  
Singlet_Tools::Singlet_Tools() {}

void Singlet_Tools::Init() {
  p_constituents = hadpars->GetConstituents();
  m_minQmass     = p_constituents->MinMass();

}

bool Singlet_Tools::CheckMass(Proto_Particle * part1,Proto_Particle * part2) {
  double factor = ((part1->Flavour().IsGluon()?2.:1.)*
		   (part2->Flavour().IsGluon()?2.:1.));
  m_mass = sqrt((part1->Momentum()+part2->Momentum()).Abs2());
  return (m_mass > (p_constituents->Mass(part1->Flavour())+
		    p_constituents->Mass(part2->Flavour())+
		    factor*m_minQmass));
}






