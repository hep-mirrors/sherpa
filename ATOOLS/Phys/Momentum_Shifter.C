#include "ATOOLS/Phys/Momentum_Shifter.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const ms::error_code code)
{
  switch (code) {
  case ms::shift_error:     return ostr<<"shift vector error";
  case ms::direction_error: return ostr<<"shift direction error";
  case ms::sperp_error:     return ostr<<"s_\\perp error";
  case ms::momenta_error:   return ostr<<"momenta error";
  case ms::boost_error_1:   return ostr<<"boost error (1)";
  case ms::boost_error_2:   return ostr<<"boost error (2)";
  case ms::scale_error:     return ostr<<"scale error";
  case ms::setup_error:     return ostr<<"setup error";
  case ms::no_error:        return ostr<<"no error";
  }
  return ostr;
}

Momentum_Shifter::Momentum_Shifter(Particle *const initial1,
				   Particle *const initial2)
{
  p_initial[0]=initial1;
  p_initial[1]=initial2;
  Reset();
}

double Momentum_Shifter::Lambda2(double sp,double sp1,double sp2) 
{ 
  double lambda2=(sp-sp1-sp2)*(sp-sp1-sp2)-4.0*sp1*sp2;
  if (!(lambda2>0.)) {
    msg_Tracking()<<"Momentum_Shifter::Lambda2("<<sp<<","<<sp1<<","<<sp2<<"): "
			  <<"\\Lambda^2(s,s_1,s_2) < 0."<<std::endl;
  }
  return lambda2;
}

bool Momentum_Shifter::CalculateShift()
{
  if (!m_setshift) {
    m_shift=Vec4D();
  }
  else {
    for (short unsigned int i=0;i<4;++i) if (IsZero(m_shift[i])) m_shift[i]=0.;
  }
  return true;
}

bool Momentum_Shifter::DetermineDirection()
{
  if (m_setshift && !m_setdirection) {
    double abs=Vec3D(m_shift).Abs();
    if (abs==0.0) {
      msg_Tracking()<<"Momentum_Shifter::DetermineDirection(): "
		    <<"Shift has vanishing 3-momentum. Abort."<<std::endl;
      return false;
    }
    m_direction=Vec4D(0.0,Vec3D(m_shift));
    m_direction=1./abs*m_direction;
    double max=0.0, sign=1.0;
    for (short unsigned int i=0;i<4;++i) {
      if (IsZero(m_direction[i])) m_direction[i]=0.;
      if (dabs(m_direction[i])>max) {
	max=dabs(m_direction[i]);
	sign=Sign(m_direction[i]);
      }
    }
    m_direction=sign*m_direction;
  }
  else if (!m_setdirection) {
    m_direction=Vec4D(0.,0.,0.,1.);
  }
  return true;
}

bool Momentum_Shifter::CalculateSPerp()
{
  m_pold[0]=Vec4D();
  for (short unsigned int i=1;i<3;++i) {
    m_pold[i]=p_initial[i-1]->Momentum();
    m_pperp[i]=Vec4D(0.0,Vec3D(m_pold[i]+(m_pold[i]*m_direction)*m_direction));
    m_pold[0]=m_pold[0]+m_pold[i];
  }
  m_pperp[0]=Vec4D(0.0,Vec3D(m_pold[0]+(m_pold[0]*m_direction)*m_direction));
  m_pnew[0]=m_pold[0]+m_shift;
  for (short unsigned int i=1;i<3;++i) {
    if (!m_setsp[i]) {
      m_sp[i]=(m_pold[i]-m_pperp[i]).Abs2();
      if (m_sp[i]<0.0) {
	msg_Tracking()<<"Momentum_Shifter::CalculateSPerp(): "
		      <<"s_{\\perp "<<i<<"} < 0. Abort."<<std::endl;
	return false;
      }
    }
  }
  m_sp[0]=(m_pnew[0]-m_pperp[0]).Abs2();
  if (m_sp[0]<0.0) {
    msg_Tracking()<<"Momentum_Shifter::CalculateSPerp(): "
		  <<"s_\\perp < 0. Abort."<<std::endl;
    return false;
  }
  return true;
}

bool Momentum_Shifter::ConstructMomenta()
{
  double E1=0.0, E2=0.0, plong1=0.0, plong2=0.0;
  double lambda2=Lambda2(m_sp[0],m_sp[1],m_sp[2]);
  if (lambda2<0.) {
    msg_Tracking()<<"Momentum_Shifter::ConstructMomenta(..): "
			  <<"\\Lambda^2("<<m_sp[0]<<","<<m_sp[1]<<","
			  <<m_sp[2]<<") < 0. Cannot shift momenta."<<std::endl;
    return false;
  }
  lambda2=sqrt(lambda2);
  double plong=-1.0*m_pnew[0]*m_direction;
  double yto=(m_pnew[0][0]+plong)/(m_pnew[0][0]-plong);
  for (double sign=1.0;sign>=-1.0;sign-=2.0) {
    E1=0.5/m_sp[0]*((m_sp[0]+m_sp[1]-m_sp[2])*m_pnew[0][0]+sign*lambda2*plong);
    E2=m_pnew[0][0]-E1;
    if (E1<0.0 || E2<0.0) continue;
    plong1=-Sign(m_pold[1]*m_direction)*sqrt(E1*E1-m_sp[1]);
    plong2=-Sign(m_pold[2]*m_direction)*sqrt(E2*E2-m_sp[2]);
    double spn1=sqr(m_pnew[0][0])-sqr(plong1+plong2);
    double spn2=sqr(m_pnew[0][0])-sqr(-plong1+plong2);
    if (dabs(spn1-m_sp[0])>dabs(spn2-m_sp[0])) plong1*=-1.0; 
    double ytn=(m_pnew[0][0]+plong1+plong2)/(m_pnew[0][0]-plong1-plong2);
    if (dabs(ytn-yto)>dabs(1.0/ytn-yto)) { 
      plong1*=-1.0; 
      plong2*=-1.0; 
    }
    if (dabs(plong1)>dabs(plong2)) { 
      if (Sign(plong1)==-Sign(m_pold[1]*m_direction)) break; 
    }
    else { 
      if (Sign(plong2)==-Sign(m_pold[2]*m_direction)) break; 
    }
  }
  if (E1<0.0 || E2<0.0) {
    m_pnew[1]=m_pold[1];
    m_pnew[2]=m_pold[2];
    return false;
  }
  m_pnew[1]=Vec4D(E1,plong1*Vec3D(m_direction))+m_pperp[1];
  m_pnew[2]=Vec4D(E2,plong2*Vec3D(m_direction))+m_pperp[2];
  return true;
}

bool Momentum_Shifter::Boost(Particle *const particle,const size_t catcher)
{
  if (m_boosted.find(particle)!=m_boosted.end()) return true;
  if (catcher>=m_maxdepth) {
    msg_Tracking()<<"Momentum_Shifter::Boost(..): "
		  <<"Nesting of event structure is deeper than "<<m_maxdepth
		  <<" levels.\n   Cannot adjust momenta."<<std::endl;
    return false;
  }
  if (particle->DecayBlob()!=NULL) {
    Blob *cur=particle->DecayBlob();
    for (size_t i=0;i<(size_t)cur->NOutP();++i) {
      if (!Boost(cur->OutParticle(i),catcher+1)) return false;
    }
  }
  Vec4D p=particle->Momentum();
  m_oldcms.Boost(p);
  m_rotate.Rotate(p);
  m_newcms.BoostBack(p);
  particle->SetMomentum(p);
  m_boosted.insert(particle);
  return true;
}

bool Momentum_Shifter::Boost(Particle *const particle)
{
  if (!m_initboost) return false;
  return Boost(particle,0);
}

bool Momentum_Shifter::BoostBack(Particle *const particle,const size_t catcher)
{
  if (m_boosted.find(particle)!=m_boosted.end()) return true;
  if (catcher>=m_maxdepth) {
    msg_Tracking()<<"Momentum_Shifter::Boost(..): "
		  <<"Nesting of event structure is deeper than "<<m_maxdepth
		  <<" levels.\n   Cannot adjust momenta."<<std::endl;
    return false;
  }
  if (particle->DecayBlob()!=NULL) {
    Blob *cur=particle->DecayBlob();
    for (size_t i=0;i<(size_t)cur->NOutP();++i) {
      if (!BoostBack(cur->OutParticle(i),catcher+1)) return false;
    }
  }
  Vec4D p=particle->Momentum();
  m_newcms.Boost(p);
  m_rotate.RotateBack(p);
  m_oldcms.BoostBack(p);
  particle->SetMomentum(p);
  m_boosted.insert(particle);
  return true;
}

bool Momentum_Shifter::BoostBack(Particle *const particle)
{
  if (!m_initboost) return false;
  return BoostBack(particle,0);
}

ms::error_code Momentum_Shifter::Boost()
{
  if (!CalculateShift()) return ms::shift_error;
  if (!DetermineDirection()) return ms::direction_error;
  if (!CalculateSPerp()) return ms::sperp_error;
  if (!ConstructMomenta()) return ms::momenta_error;
  m_oldcms=Poincare(m_pold[0]);
  m_newcms=Poincare(m_pnew[0]);
  m_newcms.Boost(m_pnew[1]); 
  m_newcms.Boost(m_pnew[2]);
  if (m_pnew[2]*m_direction<0.0) m_rotate=Poincare(Vec4D::ZVEC,m_pnew[2]);
  else m_rotate=Poincare(Vec4D::ZVEC,m_pnew[1]);
  for (short unsigned int i=1;i<3;++i) {
    m_rotate.RotateBack(m_pnew[i]);
    m_oldcms.BoostBack(m_pnew[i]);
  }
  m_boosted.clear();
  if (!Boost(p_initial[0],0)) return ms::boost_error_1;
  if (!Boost(p_initial[1],0)) {
    m_boosted.clear();
    BoostBack(p_initial[0],0);
    return ms::boost_error_2;
  }
  m_initboost=true;
  return ms::no_error;
}

ms::error_code Momentum_Shifter::BoostBack()
{
  if (!m_initboost) return ms::setup_error;
  m_boosted.clear();
  if (!BoostBack(p_initial[0],0)) return ms::boost_error_1;
  if (!BoostBack(p_initial[1],0)) {
    m_boosted.clear();
    Boost(p_initial[0],0);
    return ms::boost_error_2;
  }
  return ms::no_error;
}

ms::error_code Momentum_Shifter::Scale()
{
  if (!CalculateShift()) return ms::shift_error;
  if (!DetermineDirection()) return ms::direction_error;
  if (!CalculateSPerp()) return ms::sperp_error;
  if (!ConstructMomenta()) return ms::momenta_error;
  p_initial[0]->SetMomentum(m_pnew[1]);
  p_initial[1]->SetMomentum(m_pnew[2]);
  m_initscale=true;
  return ms::no_error;
}

void Momentum_Shifter::Reset()
{
  m_maxdepth=100;
  m_initboost=false;
  m_initscale=false;
  m_setshift=false;
  m_setdirection=false;
  for (short unsigned int i=0;i<3;++i) m_setsp[i]=false; 
}
