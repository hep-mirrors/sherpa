#include "Kinematics.H"

using namespace SHERPA;
using namespace ATOOLS;

void Kinematics::ShuffleMomenta(std::vector<Vec4D> _moms,std::vector<double> _newmasses2)
{
  int number = _moms.size();

  double xmt     = 0.;
  double oldenergies2[number];
  double energies[number];
  Vec4D cms      = Vec4D(0.,0.,0.,0.);
  for (short int i=0;i<number;i++) {
    xmt            += sqrt(_newmasses2[i]);
    cms            += _moms[i];
    oldenergies2[i] = sqr(_moms[i][0]);
  }
  double ET  = sqrt(cms.Abs2()); 
  double x   = sqrt(1.-sqr(xmt/ET));
  double acc = ET*1.e-14;

  double f0,g0,x2;
  for (int j=0;j<10;j++) {
    f0 = -ET;g0 = 0.;x2 = x*x;
    for (short int i=0;i<number;i++) {
      energies[i] = sqrt(sqr(_newmasses2[i])+x2*oldenergies2[i]);
      f0         += energies[i];
      g0         += oldenergies2[i]/energies[i];
    }
    if (dabs(f0)<acc) break; 
    x -= f0/(x*g0);  
  }
  
  // Construct Momenta
  for (short int i=0;i<number;i++) _moms[i] = Vec4D(energies[i],x*Vec3D(_moms[i]));
}


void Kinematics::ShuffleMomenta(std::vector<Particle *> _parts)
{
  int number = _parts.size();
  for (int i=0;i<number;i++) {
    if (_parts[i]->Flav().IsStable()) _parts[i]->SetFinalMass(_parts[i]->Flav().Mass()); 
  }
  if (number==2) {
    double energy   = _parts[0]->Momentum()[0]+_parts[1]->Momentum()[0];
    double energy0  = (sqr(energy)+sqr(_parts[0]->FinalMass())-sqr(_parts[1]->FinalMass()))/(2.*energy);
    double energy1  = (sqr(energy)+sqr(_parts[1]->FinalMass())-sqr(_parts[0]->FinalMass()))/(2.*energy);
    Vec3D direction = Vec3D(_parts[0]->Momentum())/(Vec3D(_parts[0]->Momentum()).Abs());
    Vec3D p0        = direction*sqrt(sqr(energy0)-sqr(_parts[0]->FinalMass()));
    Vec3D p1        = (-1.)*direction*sqrt(sqr(energy1)-sqr(_parts[1]->FinalMass()));
    _parts[0]->SetMomentum(Vec4D(energy0,p0));
    _parts[1]->SetMomentum(Vec4D(energy1,p1));
    return; 
  }

  double xmt     = 0.;
  double oldenergies2[number];
  double energies[number];
  Vec4D cms      = Vec4D(0.,0.,0.,0.);
  for (short int i=0;i<number;i++) {
    xmt            += _parts[i]->FinalMass();
    cms            += _parts[i]->Momentum();
    oldenergies2[i] = sqr(Vec3D(_parts[i]->Momentum()).Abs());
  }
  double ET  = sqrt(cms.Abs2()); 
  double x   = sqrt(1.-sqr(xmt/ET));
  double acc = ET*1.e-14;

  double f0,g0,x2;
  for (int j=0;j<10;j++) {
    f0 = -ET;g0 = 0.;x2 = x*x;
    for (short int i=0;i<number;i++) {
      energies[i] = sqrt(sqr(_parts[i]->FinalMass())+x2*oldenergies2[i]);
      f0         += energies[i];
      g0         += oldenergies2[i]/energies[i];
    }
    if (dabs(f0)<acc) break; 
    x -= f0/(x*g0);  
  }
  
  // Construct Momenta
  for (short int i=0;i<number;i++) _parts[i]->SetMomentum(Vec4D(energies[i],x*Vec3D(_parts[i]->Momentum())));
}
