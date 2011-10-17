#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace std;

bool Momenta_Stretcher::MassThem(const int n0,const int n,Vec4D * momenta,const double * masses,
				 const double rel)
{
  if ((n-n0)==2) {
    Vec4D cms         = momenta[n0]+momenta[n-1];
    Poincare boost(cms);
    for (int i=n0;i<n;i++) boost.Boost(momenta[i]);
    double energy     = momenta[n0][0]+momenta[n-1][0];
    if (masses[n0]+masses[n-1]<energy) {
      double m12      = sqr(masses[n0]);
      double m22      = sqr(masses[n-1]);
      double energy0  = (sqr(energy)+m12-m22)/(2.*energy);
      double energy1  = (sqr(energy)-m12+m22)/(2.*energy);
      Vec3D direction = Vec3D(momenta[n0])/(Vec3D(momenta[n0]).Abs());
      Vec3D p0        = direction*sqrt(sqr(energy0)-m12);
      Vec3D p1        = (-1.)*p0;
      momenta[n0]      = Vec4D(energy0,p0);
      momenta[n-1]      = Vec4D(energy1,p1);
      for (int i=n0;i<n;i++) boost.BoostBack(momenta[i]);      
      return true; 
    }
    else {
      msg_Error()<<"==================================================="<<std::endl
		 <<"Warning in "<<METHOD<<" :"<<std::endl
		 <<"   Too little energy: "<<masses[n0]<<" + "<<masses[n-1]
		 <<" > "<<energy<<"."<<std::endl;
      for (int i=n0;i<n;i++) boost.BoostBack(momenta[i]);
      return false; 
    }
  }
  else {
    double xmt         = 0.;
    double * oldens2   = new double[n];
    double * ens       = new double[n];
    Vec4D cms          = Vec4D(0.,0.,0.,0.);
    for (short int k=n0;k<n;k++) {
      xmt       += masses[k];
      cms       += momenta[k];
      oldens2[k] = sqr(momenta[k][0]);
    }
    if (cms[0]>xmt) {
      double ET  = sqrt(cms.Abs2()); 
      double x   = sqrt(1.-sqr(xmt/ET));
      double acc = dabs(rel)*ET;
      
      double f0,g0,x2;
      for (int i=0;i<10;i++) {
        f0 = -ET;g0 = 0.;x2 = x*x;
        for (short int k=n0;k<n;k++) {
          ens[k] = sqrt(sqr(masses[k])+x2*oldens2[k]);
          f0    += ens[k];
          g0    += oldens2[k]/ens[k];
        }
        if (dabs(f0)<acc) break; 
        x -= f0/(x*g0);  
      }
      for (short int k=n0;k<n;k++) {
        momenta[k] = Vec4D(ens[k],x*Vec3D(momenta[k]));
      }
      delete [] oldens2;
      delete [] ens;
      return true;
    }
    delete [] oldens2;
    delete [] ens;
    msg_Error()<<"==================================================="<<std::endl
	       <<"Warning in "<<METHOD<<" :                             "<<std::endl
	       <<"   Not enough energy ("<<cms<<") for the "<<(n-n0)
	       <<" masses ("<<xmt<<"); return false"<<std::endl
	       <<"   Masses & momenta:"<<std::endl;
    for (int i=n0;i<n;i++) msg_Error()<<"  "<<masses[i]<<" : "<<momenta[i]<<std::endl;
    return false;
    abort();
  }
  return false;
}

bool Momenta_Stretcher::MassThem(const int n0,vector<Vec4D>& momenta,vector<double> masses,
				 const double rel)
{
  int n=0;
  if(momenta.size()==masses.size()) n = momenta.size();
  else {
    return false;
  }
  if ((n-n0)==2) {
    Vec4D cms         = momenta[n0]+momenta[n-1];
    Poincare boost(cms);
    for (int i=n0;i<n;i++) boost.Boost(momenta[i]);
    double energy     = momenta[n0][0]+momenta[n-1][0];
    if (masses[n0]+masses[n-1]<energy) {
      double m12      = sqr(masses[n0]);
      double m22      = sqr(masses[n-1]);
      double energy0  = (sqr(energy)+m12-m22)/(2.*energy);
      double energy1  = (sqr(energy)-m12+m22)/(2.*energy);
      Vec3D direction = Vec3D(momenta[n0])/(Vec3D(momenta[n0]).Abs());
      Vec3D p0        = direction*sqrt(sqr(energy0)-m12);
      Vec3D p1        = (-1.)*p0;
      momenta[n0]      = Vec4D(energy0,p0);
      momenta[n-1]      = Vec4D(energy1,p1);
      for (int i=n0;i<n;i++) boost.BoostBack(momenta[i]);
      return true;
    }
    else {
      msg_Error()<<"==================================================="<<std::endl
		 <<"Warning in "<<METHOD<<" :"<<std::endl
		 <<"   Too little energy: "<<masses[n0]<<" + "<<masses[n-1]
		 <<" > "<<energy<<"."<<std::endl;
      for (int i=n0;i<n;i++) boost.BoostBack(momenta[i]);
      return false;
    }
  }
  else {
    double xmt         = 0.;
    double * oldens2   = new double[n];
    double * ens       = new double[n];
    Vec4D cms          = Vec4D(0.,0.,0.,0.);
    for (short int k=n0;k<n;k++) {
      xmt       += masses[k];
      cms       += momenta[k];
      oldens2[k] = sqr(momenta[k][0]);
    }
    if (cms[0]>xmt) {
      double ET  = sqrt(cms.Abs2());
      double x   = sqrt(1.-sqr(xmt/ET));
      double acc = dabs(rel)*ET;
      
      double f0,g0,x2;
      for (int i=0;i<10;i++) {
        f0 = -ET;g0 = 0.;x2 = x*x;
        for (short int k=n0;k<n;k++) {
          ens[k] = sqrt(sqr(masses[k])+x2*oldens2[k]);
          f0    += ens[k];
          g0    += oldens2[k]/ens[k];
        }
        if (dabs(f0)<acc) break;
        x -= f0/(x*g0);
      }
      for (short int k=n0;k<n;k++) {
        momenta[k] = Vec4D(ens[k],x*Vec3D(momenta[k]));
      }
      delete [] oldens2;
      delete [] ens;
      return true;
    }
    delete [] oldens2;
    delete [] ens;
    msg_Error()<<"==================================================="<<std::endl
	       <<"Warning in "<<METHOD<<" : "<<endl
	       <<"   Not enough energy ("<<cms<<") for the "
	       <<(n-n0)<<" masses ("<<xmt<<"); return false"<<endl
	       <<"   Masses & momenta:"<<endl;
    for (int i=n0;i<n;i++) msg_Error()<<masses[i]<<" : "<<momenta[i]<<std::endl;
    return false;
    abort();
  }
  return false;
}

bool Momenta_Stretcher::ZeroThem(const int n0, const int n, Vec4D * momenta,
				 const double rel)
{
  if ((n-n0)==2) {
    double energy   = momenta[n0][0]+momenta[n-1][0];
    Vec3D direction = Vec3D(momenta[n0])/(Vec3D(momenta[n0]).Abs());
    momenta[n0]      = energy/2.*Vec4D(1.,direction);
    momenta[n-1]      = energy/2.*Vec4D(1.,-1.*direction);
    return true;
  }
  else {
    double xmt         = 0.;
    double * oldps2    = new double[n];
    double * ens       = new double[n];
    Vec4D cms          = Vec4D(0.,0.,0.,0.);
    for (short int i=n0;i<n;i++) {
      xmt      += sqrt(Max(0.,momenta[i].Abs2()));
      oldps2[i] = sqr(Vec3D(momenta[i]).Abs());
      cms       += momenta[i];
    }
    double ET  = sqrt(cms.Abs2());
    double x   = 1./sqrt(1.-sqr(xmt/ET));
    double acc = dabs(rel)*ET;
    xmt        = 0.;

    double f0,g0,x2;
    for (int i=0;i<10;i++) {
      f0 = -ET;g0 = 0.;x2 = x*x;
      for (short int i=n0;i<n;i++) {
        ens[i] = sqrt(x2*oldps2[i]);
        f0    += ens[i];
        g0    += oldps2[i]/ens[i];
      }
      if (dabs(f0)<acc) break;
      x -= f0/(x*g0);
    }
    for (short int k=n0;k<n;k++) {
      momenta[k] = Vec4D(ens[k],x*Vec3D(momenta[k]));
    }
    delete [] oldps2;
    delete [] ens;
    return true;
  }
  return false;
}


bool Momenta_Stretcher::ZeroThem(const int n0,vector<Vec4D>& momenta,
				 const double rel)
{
  int n = momenta.size();
  
  if ((n-n0)==2) {
    double energy   = momenta[n0][0]+momenta[n-1][0];
    Vec3D direction = Vec3D(momenta[n0])/(Vec3D(momenta[n0]).Abs());
    momenta[n0]      = energy/2.*Vec4D(1.,direction);
    momenta[n-1]      = energy/2.*Vec4D(1.,-1.*direction);
    return true;
  }
  else {
    double xmt         = 0.;
    double * oldps2    = new double[n];
    double * ens       = new double[n];
    Vec4D cms          = Vec4D(0.,0.,0.,0.);
    for (short int i=n0;i<n;i++) {
      xmt      += sqrt(Max(0.,momenta[i].Abs2()));
      oldps2[i] = sqr(Vec3D(momenta[i]).Abs());
      cms       += momenta[i];
    }
    double ET  = sqrt(cms.Abs2());
    double x   = 1./sqrt(1.-sqr(xmt/ET));
    double acc = dabs(rel)*ET;
    xmt        = 0.;
    
    double f0,g0,x2;
    for (int i=0;i<10;i++) {
      f0 = -ET;g0 = 0.;x2 = x*x;
      for (short int i=n0;i<n;i++) {
        ens[i] = sqrt(x2*oldps2[i]);
        f0    += ens[i];
        g0    += oldps2[i]/ens[i];
      }
      if (dabs(f0)<acc) break;
      x -= f0/(x*g0);
    }
    for (short int k=n0;k<n;k++) {
      momenta[k] = Vec4D(ens[k],x*Vec3D(momenta[k]));
    }
    delete [] oldps2;
    delete [] ens;
    return true;
  }
  return false;
}

bool Momenta_Stretcher::StretchBlob(Blob* blob)
{
  if(blob->GetOutParticles().size()<2) return true;
  std::vector<double> masses;
  Particle_Vector outparts = blob->GetOutParticles();
  vector<Vec4D> momenta;
  Vec4D total(0.0,0.0,0.0,0.0);
  for (size_t i=0;i <outparts.size();i++) {
    if (outparts[i]->DecayBlob()&&outparts[i]->DecayBlob()->NOutP()>0) continue;
    masses.push_back(outparts[i]->FinalMass());
    momenta.push_back(outparts[i]->Momentum());
    total+=outparts[i]->Momentum();
// =======
//   //msg_Out()<<"Check the "<<outparts.size()<<" momenta of blob in "<<METHOD<<":"<<std::endl;
//   for(Particle_Vector::iterator pit=outparts.begin();pit!=outparts.end();pit++) {
//     if( use_finalmasses ) masses.push_back( (*pit)->FinalMass() );
//     momenta.push_back( (*pit)->Momentum() );
// //     msg_Out()<<"  "<<(*pit)->Flav()<<" "<<(*pit)->FinalMass()<<" "<<(*pit)->Momentum()<<std::endl;
// >>>>>>> .merge-right.r13247
  }
  Poincare cms(total);
  for (size_t i=0; i<momenta.size(); ++i) cms.Boost(momenta[i]);
  if(!ZeroThem(0,momenta)) return false;
  if(!MassThem(0,momenta,masses)) {
    msg_Error()<<"Error in "<<METHOD<<"(Blob *)."<<std::endl;
    return false;
  }
  size_t j=0;
  for(size_t i=0;i<outparts.size();i++) {
    if (outparts[i]->DecayBlob()&&outparts[i]->DecayBlob()->NOutP()>0) continue;
    cms.BoostBack(momenta[j]);
    outparts[i]->SetMomentum(momenta[j]);
    ++j;
  }
  return true;
}

bool Momenta_Stretcher::StretchMomenta( const Particle_Vector& outparts, std::vector<Vec4D>& moms)
{
  if(outparts.size() != moms.size()) return false;
  if(outparts.size()==1 && abs(outparts[0]->FinalMass()-moms[0].Mass())<Accu() ) return true;

  Vec4D cms;
  vector<double> masses;
  for(size_t i=0; i<moms.size(); i++) {
    cms += moms[i];
    masses.push_back(outparts[i]->FinalMass());
  }
  Poincare boost(cms);
  for(size_t i=0; i<moms.size(); i++) {
    moms[i] = boost*moms[i];
  }
  if(!ZeroThem(0,moms)) return false;
  if(! MassThem(0,moms,masses)) {
    msg_Error()<<"Error in "<<METHOD<<"(const Particle_Vector&, moms)."<<std::endl;
    return false;
  }
  boost.Invert();
  for(size_t i=0; i<moms.size(); i++) {
    moms[i] = boost*moms[i];
  }
  return true;
}

bool Momenta_Stretcher::StretchMomenta( const Particle_Vector& outparts,
                                        std::vector<double>& masses )
{
  if(outparts.size() != masses.size()) return false;
  if(outparts.size()==1 && abs(outparts[0]->FinalMass()-masses[0])<Accu() ) return true;

  Vec4D cms;
  vector<Vec4D> moms;
  for(size_t i=0; i<masses.size(); i++) {
    moms.push_back(outparts[i]->Momentum());
    cms += moms[i];
  }
  Poincare boost(cms);
  for(size_t i=0; i<masses.size(); i++) {
    boost.Boost(moms[i]);
  }
  if(!ZeroThem(0,moms)) return false;
  if(!MassThem(0,moms,masses)) {
    msg_Error()<<"Warning in "<<METHOD<<"(const Particle_Vector&, masses)."<<std::endl;
    return false;
  }
  for(size_t i=0; i<moms.size(); i++) {
    boost.BoostBack(moms[i]);
    outparts[i]->SetMomentum(moms[i]);
    outparts[i]->SetFinalMass(masses[i]);
  }
  return true;
}
