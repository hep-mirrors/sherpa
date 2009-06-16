#include "ATOOLS/Phys/Particle_List.H"

#include "ATOOLS/Math/Poincare.H"

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &s,const Particle_List &pl) 
{
  s<<"Particle List with "<<pl.size()<<" elements"<<std::endl;
  for (Particle_List::const_iterator pit=pl.begin(); pit!=pl.end(); ++pit) {
    if (*pit!=NULL) s<<**pit<<"\n";
    else s<<"NULL pointer\n";
  }
  return s;
}

Particle_List::Particle_List():
  m_destructor(NULL) {}

Particle_List::Particle_List(const bool destruct):
  m_destructor(destruct?this:NULL) {}

void Particle_List::Clear()
{
  while (!empty()) {
    delete back();
    pop_back();
  }
}

void Particle_List::Boost(Poincare *const boost) const
{
  for (const_iterator pit=begin();pit!=end();++pit) {
    Vec4D mom((*pit)->Momentum());
    boost->Boost(mom);
    (*pit)->SetMomentum(mom);
  }
}

void Particle_List::BoostBack(Poincare *const boost) const
{
  for (const_iterator pit=begin();pit!=end();++pit) {
    Vec4D mom((*pit)->Momentum());
    boost->BoostBack(mom);
    (*pit)->SetMomentum(mom);
  }
}

void Particle_List::Rotate(Poincare *const rot) const
{
  for (const_iterator pit=begin();pit!=end();++pit) {
    Vec4D mom((*pit)->Momentum());
    rot->Rotate(mom);
    (*pit)->SetMomentum(mom);
  }
}

void Particle_List::RotateBack(Poincare *const rot) const
{
  for (const_iterator pit=begin();pit!=end();++pit) {
    Vec4D mom((*pit)->Momentum());
    rot->RotateBack(mom);
    (*pit)->SetMomentum(mom);
  }
}

void Particle_List::Flip() const
{
  for (const_iterator pit=begin();pit!=end();++pit) {
    Vec4D mom((*pit)->Momentum());
    (*pit)->SetMomentum(Vec4D(mom[0],(-1.)*Vec3D(mom)));
  }
}
