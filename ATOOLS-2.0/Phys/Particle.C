#include "Particle.H"
#include "Blob.H"
#include "Random.H"
#include "Run_Parameter.H"
#include <iomanip>

using namespace ATOOLS;

std::ostream& ATOOLS::operator<<(std::ostream& str,Particle* part) {
  str<<std::setprecision(4)<<std::setiosflags(std::ios::left);
  switch (part->Status()) {
  case 0 : // null entry
    return str<<"--- empty entry ---"<<std::endl;
  case 1 : // active (final state) particle
  case 2 : // decayed or fragmented
    str<<"  "<<std::setw(3)<<part->Info()<<"  "<<std::setw(3)<<part->Status()<<std::setw(1)<<" "
       <<std::setw(22)<<part->Flav()<<std::setw(1)<<" "
       <<std::setw(10)<<part->Number()<<std::setw(1)<<" (";
    if (part->ProductionBlob()) str<<std::setw(5)<<part->ProductionBlob()->Id();
                           else str<<"     ";
    if (part->DecayBlob()) str<<" -> "<<std::setw(5)<<part->DecayBlob()->Id();
                      else str<<" ->      ";
    str<<std::setw(1)<<") ";
    break;
  case 3 : // documentation line
    return str<<"============================================================"<<std::endl;
  default : // user defined or reserved
    return str<<"--- unrecognized status:"<<part->Status()<<" ---"<<std::endl;
  }
  str<<std::setiosflags(std::ios::scientific)
     <<" ["<<part->Momentum()<<", "<<part->Momentum().Abs2()<<"]"
     <<" ("<<part->GetFlow(1)<<","<<part->GetFlow(2)<<")";
  return str;
}

Particle::~Particle() {
  if (p_flow) { delete p_flow; p_flow = NULL; } 
}

Particle::Particle() {
  m_number    = -1;
  m_info      = ' ';
  m_status    = 1;
  m_fl        = Flavour(kf::none);
  m_momentum  = Vec4D(0,0,0,0); 
  m_dec_time  = 0.;
  p_startblob = NULL;
  p_endblob   = NULL;
  p_flow      = new Flow(this);
}

Particle::Particle(const Particle * in)  {
  m_number    = in->m_number;
  m_info      = in->m_info;
  m_status    = in->Status();
  m_fl        = in->m_fl;
  m_momentum  = in->m_momentum;
  m_dec_time  = in->m_dec_time;
  p_startblob = in->p_startblob;
  p_endblob   = in->p_endblob;
  p_flow      = new Flow(this);
  p_flow->SetCode(1,in->GetFlow(1));
  p_flow->SetCode(2,in->GetFlow(2));
}

Particle::Particle(const Particle & in)  
{
  p_flow=0;
  *this = in;
}

Particle& Particle::operator=(const Particle & in)
{
  if (this!=&in) {
    m_number    = in.m_number;
    m_info      = in.m_info;
    m_status    = in.Status();
    m_fl        = in.m_fl;
    m_momentum  = in.m_momentum;
    m_dec_time  = in.m_dec_time;
    p_startblob = in.p_startblob;
    p_endblob   = in.p_endblob;
    if (p_flow) delete p_flow;
    p_flow      = new Flow(this);
    p_flow->SetCode(1,in.GetFlow(1));
    p_flow->SetCode(2,in.GetFlow(2));

  }
  return *this;
}


Particle::Particle(int number,Flavour fl,Vec4D p)  {
  m_number    = number;
  m_status    = 1;
  m_info      = ' ';
  m_fl        = fl;
  m_momentum  = p;
  m_dec_time  = 0.;
  p_startblob = 0;
  p_endblob   = 0;
  p_flow      = new Flow(this);
}

void Particle::Copy(Particle * in)  {
  m_number    = in->m_number;
  m_info      = in->m_info;
  m_status    = in->Status();
  m_fl        = in->m_fl;
  m_momentum  = in->m_momentum;
  m_dec_time  = in->m_dec_time;
  p_startblob = in->p_startblob;
  p_endblob   = in->p_endblob;
  if (!p_flow)  p_flow  = new Flow(this);
  p_flow->SetCode(1,in->GetFlow(1));
  p_flow->SetCode(2,in->GetFlow(2));
}

double Particle::ProperTime() {
  double q2    = m_momentum.Abs2();
  double m2    = sqr(m_fl.Mass());
  double tau2  = 1.e96;
  if (!( (q2-m2 < rpa.gen.Accu()) && (m_fl.Width() < rpa.gen.Accu()))) { 
    if (m2>rpa.gen.Accu()) { 
      tau2 = q2/(sqr(q2-m2)+sqr(q2*m_fl.Width())/m2);
    }
    else {
      if (dabs(q2)>rpa.gen.Accu()) tau2 = 1/dabs(q2);
    }
  }
  else {
    if (m_fl.Strong()) tau2 = 1./sqr(0.2); 
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

Vec3D Particle::Distance() {
  Vec3D v = Vec3D(m_momentum)/E()*rpa.c();
  return v*LifeTime();
}

// Numbers etc.
int    Particle::Number()   const                 { return m_number; }
void   Particle::SetNumber(const int n)           { m_number    = n; }
int    Particle::JetNumber()   const              { return m_jetnumber; }
void   Particle::SetJetNumber(const int n)        { m_jetnumber = n; }

  // Status etc.
int    Particle::Status() const                   { return m_status; }
void   Particle::SetStatus( int status )          { m_status = status; }
char   Particle::Info() const                     { return m_info;}
void   Particle::SetInfo(char info)               { m_info = info; }

  // Momentum, energy, and lifetime
Vec4D  Particle::Momentum() const                 { return m_momentum; }
double Particle::E()                              { return m_momentum[0];}
void   Particle::SetMomentum(const Vec4D & vec4 ) { m_momentum = vec4; } 
double Particle::Time() const                     { return m_dec_time; }
void   Particle::SetTime(const int t)             { m_dec_time = t; }
void   Particle::SetTime()                        { m_dec_time = LifeTime(); }

// Production and decay vertices
Vec4D  Particle::XProd()                          { return p_startblob->Position(); }
Blob * Particle::ProductionBlob()                 { return p_startblob; }
void   Particle::SetProductionBlob(Blob * _blob)  { p_startblob = _blob; }
Vec4D  Particle::XDec()                           { return p_endblob->Position(); }
Blob * Particle::DecayBlob()                      { return p_endblob; }
void   Particle::SetDecayBlob(Blob * _blob)       { p_endblob = _blob; }

// Flavour and flow
Flavour   Particle::Flav() const                     { return m_fl; }
void      Particle::SetFlav(Flavour & fl) { m_fl      = fl; }
Flow    * Particle::GetFlow() const                  { return p_flow; }
int       Particle::GetFlow(const int index) const   { return p_flow->Code(index); }
void      Particle::SetFlow(Flow * _flow)            { p_flow    = _flow; }
void      Particle::SetFlow(const int index, const int code) {
  if ((!m_fl.IsDiQuark()) && (!m_fl.Strong())) return;
  p_flow->SetCode(index,code);
}



void Particle2MPI(const Particle * p , MPI_Particle & mpi_p) {
  mpi_p.id  =p->Number();
  mpi_p.m_fl=int(p->Flav());
  for (int i=0; i<4; ++i)  
  mpi_p.m_mom[i]=p->Momentum()[i];
  mpi_p.m_flow[0]=p->GetFlow(1);
  mpi_p.m_flow[1]=p->GetFlow(2);
}
  
Particle * MPI2Particle(const MPI_Particle & mpi_p ) {
  Particle * p ;
  if (mpi_p.m_fl>0) p= new Particle(mpi_p.id, Flavour((kf::code)mpi_p.m_fl),
			  Vec4D(mpi_p.m_mom[0],mpi_p.m_mom[1],mpi_p.m_mom[2],mpi_p.m_mom[3]));
  else p= new Particle(mpi_p.id, Flavour((kf::code)(-mpi_p.m_fl)).Bar(),
			  Vec4D(mpi_p.m_mom[0],mpi_p.m_mom[1],mpi_p.m_mom[2],mpi_p.m_mom[3]));
  if (mpi_p.m_flow[0]) p->SetFlow(1,mpi_p.m_flow[0]);
  if (mpi_p.m_flow[1]) p->SetFlow(2,mpi_p.m_flow[1]);
  return p;
}
