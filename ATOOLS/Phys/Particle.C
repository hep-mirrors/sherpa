#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Smart_Pointer.C"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

#ifdef PROFILE__all
#include "prof.hh"
#else
#ifdef PROFILE__Particle
#include "prof.hh"
#else
#define PROFILE_HERE
#endif
#endif

namespace ATOOLS {
  int Particle::s_totalnumber=0;
  long unsigned int Particle::s_currentnumber=0;
}

using namespace ATOOLS;

namespace ATOOLS { template class SP(Part_List); }

bool ATOOLS::Particle::operator==(Particle part)
{
  if ((part.m_status==m_status)&&         
      (part.m_info==m_info)&&
      (part.m_fl==m_fl)&&         
      (part.m_momentum==m_momentum)&&
      (part.m_dec_time==m_dec_time)) {
    return true;
  }
  return false;
}

std::ostream& ATOOLS::operator<<(std::ostream& str, const Particle &part) {
  int io;
  switch (part.Status()) {
  case part_status::undefined : // null entry
    return str<<"--- empty entry ---"<<std::endl;
  case part_status::active :  // active (final state) particle
  case part_status::decayed : // decayed particle 
  case part_status::fragmented : // or fragmented particle
    io=str.precision(4);
    str<<std::setiosflags(std::ios::left);
    str<<"["<<part.Info()<<"] "<<part.Status()<<" "
       <<std::setw(16)<<part.Flav()<<" "<<std::setiosflags(std::ios::right)
       <<std::setw(6)<<part.Number()<<" (";
    if (part.ProductionBlob()) str<<std::setw(4)<<part.ProductionBlob()->Id();
    else str<<"    ";
    if (part.DecayBlob()) str<<" -> "<<std::setw(4)<<part.DecayBlob()->Id();
    else str<<" ->     ";
    str<<")"<<std::resetiosflags(std::ios::right);
    str<<std::resetiosflags(std::ios::scientific)<<std::resetiosflags(std::ios::left);
    break;
  case part_status::documentation : // documentation line
    io=str.precision(4);
    str<<std::setiosflags(std::ios::left);
    str<<"============================================================"<<std::endl
       <<"  "<<std::setw(3)<<part.Info()<<"  "<<std::setw(3)<<part.Status()<<std::setw(1)<<" "
       <<std::setw(22)<<part.Flav()<<std::setw(1)<<" "
       <<std::setw(10)<<part.Number()<<std::endl
       <<"============================================================"
       <<std::resetiosflags(std::ios::scientific)<<std::resetiosflags(std::ios::left);
    str.precision(io);
    return str;		  
  default : // user defined or reserved
    return str<<"--- unrecognized status:"<<part.Status()<<" ---"<<std::endl;
  }
  Vec4D p=part.Momentum();
  str<<std::setiosflags(std::ios::scientific)
     <<" [("<<std::setw(11)<<p[0]<<','<<std::setw(11)<<p[1]<<','
     <<std::setw(11)<<p[2]<<','<<std::setw(11)<<p[3]<<"), p^2="
     <<std::setw(11)<<p.Abs2()<<", m="<<std::setw(11)<<part.m_finalmass<<"]"
     <<" ("<<std::setw(3)<<part.GetFlow(1)<<","<<std::setw(3)<<part.GetFlow(2)<<")"
     <<std::resetiosflags(std::ios::scientific)<<std::resetiosflags(std::ios::left);
  str.precision(io);
  return str;
}

namespace ATOOLS {
  std::ostream & operator<<(std::ostream & s, const Part_List & pl) {
    s<<"Particle List with "<<pl.size()<<" elements"<<std::endl;
    for (Part_Const_Iterator pit=pl.begin(); pit!=pl.end(); ++pit) {
      s<<**pit<<std::endl;
    }
    return s;
  }
}


Particle::~Particle() 
{
  delete p_flow; 
  --s_totalnumber;
}

Particle::Particle():
  m_number(-1), m_jetnumber(-1), m_status(part_status::undefined), 
  m_info('X'), 
  m_fl(Flavour(kf_none)), m_momentum(Vec4D(0,0,0,0)), 
  p_flow(new Flow(this)),
  p_startblob(NULL),p_endblob(NULL), p_originalpart(this),
  m_dec_time(0.), m_finalmass(0.)
{
  ++s_totalnumber;
}

Particle::Particle(const Particle &in): 
  m_number(in.m_number), m_jetnumber(in.m_jetnumber), m_status(in.m_status), 
  m_info(in.m_info), 
  m_fl(in.m_fl), m_momentum(in.m_momentum), 
  p_flow(new Flow(this)),
  p_startblob(NULL),p_endblob(NULL), p_originalpart(in.p_originalpart),
  m_dec_time(in.m_dec_time), m_finalmass(in.m_finalmass)
{
  ++s_totalnumber;
  p_flow->SetCode(1,in.GetFlow(1));
  p_flow->SetCode(2,in.GetFlow(2));
}

Particle& Particle::operator=(const Particle &in)
{
  if (this!=&in) {
    m_number    = in.m_number;
    m_info      = in.m_info;
    m_status    = in.m_status;
    m_fl        = in.m_fl;
    m_momentum  = in.m_momentum;
    m_dec_time  = in.m_dec_time;
    m_finalmass = in.m_finalmass;
    p_startblob = NULL;
    p_endblob   = NULL;
    p_flow->SetCode(1,in.GetFlow(1));
    p_flow->SetCode(2,in.GetFlow(2));
  }
  return *this;
}


Particle::Particle(int number, Flavour fl, Vec4D p, char a) :
  m_number(number), m_jetnumber(-1), m_status(part_status::active),
  m_info(a), 
  m_fl(fl), m_momentum(p),
  p_flow(new Flow(this)),
  p_startblob(NULL),p_endblob(NULL), p_originalpart(this),
  m_dec_time(0.), m_finalmass(fl.Mass())
{
  ++s_totalnumber;
}


void Particle::Copy(Particle * in)  {
  m_number    = in->m_number;
  m_info      = in->m_info;
  m_status    = in->m_status;
  m_fl        = in->m_fl;
  m_momentum  = in->m_momentum;
  m_dec_time  = in->m_dec_time;
  m_finalmass = in->m_finalmass;
  p_startblob = in->p_startblob;
  p_endblob   = in->p_endblob;
  p_originalpart = in->p_originalpart,
  p_flow->SetCode(1,in->GetFlow(1));
  p_flow->SetCode(2,in->GetFlow(2));
}

double Particle::ProperTime() 
{
  double q2    = m_momentum.Abs2();
  double m2    = sqr(m_fl.Mass());
  double tau2  = 1.e96;
  if (!( (q2-m2 < rpa.gen.Accu()) && 
         (m_fl.Width() < rpa.gen.Accu()) )) { // off-shell or big width
    if (m2>rpa.gen.Accu()) { 
      tau2 = q2/(sqr(q2-m2)+sqr(q2*m_fl.Width())/m2);
    }
    else {
      if (dabs(q2)>rpa.gen.Accu()) tau2 = 1/dabs(q2);
    }
  }
  else {
    if (m_fl.Strong()) tau2 = 1./sqr(0.2); 
    else if (!m_fl.IsStable()) tau2 = 1./sqr(m_fl.Width());
  }
  return rpa.hBar() * sqrt(tau2);
}

double Particle::LifeTime() {
  double t   = -ProperTime()*log(1.-ran.Get());  
  if (t>1.e6) t = 1.e6;
  double gamma = 1./rpa.gen.Accu();
  if (m_fl.Mass()>rpa.gen.Accu()) gamma = E()/m_fl.Mass();
  else {
    double q2    = dabs(m_momentum.Abs2());
    if (q2>rpa.gen.Accu()) gamma = E()/sqrt(q2);
  }
  return gamma * t;      
}

Vec3D Particle::Distance(double _lifetime) {
  Vec3D v = Vec3D(m_momentum)/E()*rpa.c();
  if (_lifetime<0.) _lifetime = LifeTime();
  return v*_lifetime;
}

void Particle::SetProductionBlob(Blob *blob)  
{ 
  if (p_startblob!=NULL && blob!=NULL) {
    if (p_startblob->Id()>-1) 
      msg_Out()<<"WARNING in Particle::SetProductionBlob("<<blob<<"):"<<std::endl
		       <<"   blob->Id() = "<<blob->Id()<<std::endl
		       <<"   Particle already has a production blob!"<<std::endl
		       <<"   "<<*this<<std::endl;
  }
  p_startblob=blob; 
}

// Numbers etc.
int  Particle::Number() const                   { return m_number; }
int  Particle::JetNumber() const                { return m_jetnumber; }
void Particle::SetJetNumber(const int n)        { m_jetnumber = n; }

// Status etc.
part_status::code Particle::Status() const               { return m_status; }
void Particle::SetStatus(const part_status::code status) { m_status = status; }
char              Particle::Info() const                 { return m_info;}
void              Particle::SetInfo(char info)           { m_info = info; }

// Momentum, energy, and lifetime
const Vec4D& Particle::Momentum() const              { return m_momentum; }
double       Particle::E() const                     { return m_momentum[0];}
double       Particle::FinalMass() const             { return m_finalmass; }
void         Particle::SetMomentum(const Vec4D& vc4) { m_momentum = vc4; }
double       Particle::Time() const                  { return m_dec_time; }
void         Particle::SetTime(const double t)       { m_dec_time = t; }
void         Particle::SetTime()                     { m_dec_time = m_fl.DiceLifeTime(); }

// Production and decay vertices
Vec4D Particle::XProd() const
{ if (p_startblob) return p_startblob->Position(); return Vec4D(); }
Blob *       Particle::ProductionBlob() const {return p_startblob;}
Vec4D Particle::XDec() const
{ if (p_endblob) return p_endblob->Position(); return Vec4D(); }
Blob *       Particle::DecayBlob() const      {return p_endblob;}
Particle *   Particle::OriginalPart() const   {return p_originalpart;}

// Flavour and flow
Flavour        Particle::Flav() const                   { return m_fl; }
const Flavour& Particle::RefFlav() const                { return m_fl; }
void           Particle::SetFlav(const Flavour& fl)     { m_fl   = fl; }
Flow         * Particle::GetFlow() const                { return p_flow; }
unsigned int   Particle::GetFlow(const unsigned int index) const {
  return p_flow->Code(index);
}
void           Particle::SetFlow(Flow * _flow)          { p_flow = _flow; }
void           Particle::SetFlow(const int index, const int code) {
  if ((!m_fl.IsDiQuark()) && (!m_fl.Strong())) return;
  p_flow->SetCode(index,code);
}


void Particle::SetDecayBlob(Blob *blob)       
{ 
  p_endblob=blob; 
}

void Particle::SetOriginalPart(Particle *part)
{ 
  p_originalpart=part; 
}

void Particle::SetNumber(const int n)           
{ 
  if (n<0) m_number = -n;
  else {
    if (m_number<=0) m_number=++s_currentnumber;
  }
}

void   Particle::SetFinalMass(const double _lower,const double _upper) {
  if (_lower==-1. &&_upper==-1.) {
    m_finalmass = m_fl.HadMass();
    return;
  }
  if (_upper<0.) {
    m_finalmass = _lower;
    return;
  }
  double mass2  = m_fl.Mass()*m_fl.Mass();
  double mw     = m_fl.Mass()*m_fl.Width();
  double low2   = _lower*_lower;
  double up2    = _upper*_upper;
  double range  = up2-low2;
  double yup    = (up2-mass2)/mw;
  double ylow   = (low2-mass2)/mw;
  double ymin   = atan(ylow);
  double yrange = atan(range/(mw*(1.+ylow*yup)));
  if (ylow*yup<-1.) {
    if (yup>0) yrange = yrange + M_PI;
    if (yup<0) yrange = yrange - M_PI;
  }     
  m_finalmass = sqrt(mass2+mw*tan(ran.Get()*yrange + ymin));
}


















