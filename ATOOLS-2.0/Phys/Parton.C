#include "Parton.H"
#include "Blob.H"
#include "Random.H"
#include "Run_Parameter.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;

std::ostream& APHYTOOLS::operator<<(std::ostream& str,Parton* part) {
  switch (part->Status()) {
  case 0 : // null entry
    return str<<"--- empty entry ---"<<std::endl;
  case 1 : // active (final state) particle
    str<<part->Info();
    str.width(8);
    str<<part->Flav();
    str<<' ';
    str<<part->Number();
    str<<' ';
    if (part->ProductionBlob()) str<<" produced in "<<part->ProductionBlob()->Id()<<std::endl<<" ";
    if (part->DecayBlob())      str<<" decayed in "<<part->DecayBlob()->Id()<<std::endl<<" ";
    break;
  case 2 : // decayed or fragmented
    str<<part->Info();
    str<<'!';
    str.width(7);
    str<<part->Flav();
    str<<'!';
    str<<part->Number();
    str<<' ';
    if (part->ProductionBlob()) str<<" produced in "<<part->ProductionBlob()->Id()<<std::endl<<" ";
    if (part->DecayBlob())      str<<" decayed in "<<part->DecayBlob()->Id()<<std::endl<<" ";
    break;
  case 3 : // documentation line
    return str<<"============================================================"<<std::endl;
  default : // user defined or reserved
    return str<<"--- unrecognized status:"<<part->Status()<<" ---"<<std::endl;
  }
  str<<"      "<<part->Momentum()<<" "<<part->Momentum().Abs2()
     <<" colours : ("<<part->GetFlow(1)<<","<<part->GetFlow(2)<<")"<<std::endl;
  return str;
}

Parton::~Parton() {
  if (p_flow) { delete p_flow; p_flow = NULL; } 
}

Parton::Parton() {
  m_number    = -1;
  m_info      = ' ';
  m_status    = 1;
  m_fl        = Flavour(kf::none);
  m_momentum  = Vec4D(0,0,0,0); 
  m_dec_time  = 0.;
  p_startblob = NULL;
  p_endblob   = NULL;
  cout<<"Blobs done."<<endl;
  p_flow      = new Flow(this);
}

Parton::Parton(Parton * in)  {
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

Parton::Parton(int number,Flavour fl,Vec4D p)  {
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

void Parton::Copy(Parton* in)  {
  m_number    = in->m_number;
  m_info      = in->m_info;
  m_status    = in->Status();
  m_fl        = in->m_fl;
  m_momentum  = in->m_momentum;
  m_dec_time  = in->m_dec_time;
  p_startblob = in->p_startblob;
  p_endblob   = in->p_endblob;
  p_flow      = in->p_flow;
}

double Parton::ProperTime() {
  double q2    = m_momentum.Abs2();
  double m2    = sqr(m_fl.Mass());
  double tau2  = 1.e6;
  if (!( (q2-m2 < rpa.gen.Accu()) && (m_fl.Width() < rpa.gen.Accu()))) { 
    if (m2>AORGTOOLS::rpa.gen.Accu()) { 
      tau2 = q2/(sqr(q2-m2)+sqr(q2*m_fl.Width())/m2);
    }
    else {
      if (dabs(q2)>AORGTOOLS::rpa.gen.Accu()) tau2 = 1/dabs(q2);
    }
  }
  else {
    if (m_fl.Strong()) tau2 = 1./sqr(0.2); 
  }
  return AORGTOOLS::rpa.hBar() * sqrt(tau2);
}

double Parton::LifeTime() {
  double t   = -ProperTime()*log(1.-ran.Get());  
  if (t>1.e6) t = 1.e6;
  double gamma = 1./AORGTOOLS::rpa.gen.Accu();
  if (m_fl.Mass()>AORGTOOLS::rpa.gen.Accu()) gamma = E()/m_fl.Mass();
  return gamma * t;      
}

Vec3D Parton::Distance() {
  double gamma     = E()/m_fl.Mass();
  Vec3D v = Vec3D(m_momentum)/E()*AORGTOOLS::rpa.c();
  return v*LifeTime();
}

// Numbers etc.
int    Parton::Number()   const                 { return m_number; }
void   Parton::SetNumber(const int n)           { m_number    = n; }
int    Parton::JetNumber()   const              { return m_jetnumber; }
void   Parton::SetJetNumber(const int n)        { m_jetnumber = n; }

  // Status etc.
int    Parton::Status() const                   { return m_status; }
void   Parton::SetStatus( int status )          { m_status = status; }
char   Parton::Info() const                     { return m_info;}
void   Parton::SetInfo(char info)               { m_info = info; }

  // Momentum, energy, and lifetime
Vec4D  Parton::Momentum() const                 { return m_momentum; }
double Parton::E()                              { return m_momentum[0];}
void   Parton::SetMomentum(const Vec4D & vec4 ) { m_momentum = vec4; } 
double Parton::Time() const                     { return m_dec_time; }
void   Parton::SetTime(const int t)             { m_dec_time = t; }
void   Parton::SetTime()                        { m_dec_time = LifeTime(); }

// Production and decay vertices
Vec4D  Parton::XProd()                          { return p_startblob->Position(); }
Blob * Parton::ProductionBlob()                 { return p_startblob; }
void   Parton::SetProductionBlob(Blob * _blob)  { p_startblob = _blob; }
Vec4D  Parton::XDec()                           { return p_endblob->Position(); }
Blob * Parton::DecayBlob()                      { return p_endblob; }
void   Parton::SetDecayBlob(Blob * _blob)       { p_endblob = _blob; }

// Flavour and flow
Flavour   Parton::Flav() const                     { return m_fl; }
void      Parton::SetFlav(APHYTOOLS::Flavour & fl) { m_fl      = fl; }
Flow    * Parton::GetFlow() const                  { return p_flow; }
int       Parton::GetFlow(const int index) const   { return p_flow->Code(index); }
void      Parton::SetFlow(Flow * _flow)            { p_flow    = _flow; }
void      Parton::SetFlow(const int index, const int code) {
  cout<<"In Parton SetFlow("<<index<<","<<code<<")"<<endl;
  if (!m_fl.Strong()) return;
  if (code==-1) p_flow->SetUniqueCode(index);
  else p_flow->SetCode(index,code);
}



void Parton2MPI(const Parton * p , MPI_Parton & mpi_p) {
  mpi_p.id  =p->Number();
  mpi_p.m_fl=int(p->Flav());
  for (int i=0; i<4; ++i)  
  mpi_p.m_mom[i]=p->Momentum()[i];
  mpi_p.m_flow[0]=p->GetFlow(1);
  mpi_p.m_flow[1]=p->GetFlow(2);
}
  
Parton * MPI2Parton(const MPI_Parton & mpi_p ) {
  Parton * p ;
  if (mpi_p.m_fl>0) p= new Parton(mpi_p.id, Flavour((kf::code)mpi_p.m_fl),
			  Vec4D(mpi_p.m_mom[0],mpi_p.m_mom[1],mpi_p.m_mom[2],mpi_p.m_mom[3]));
  else p= new Parton(mpi_p.id, Flavour((kf::code)(-mpi_p.m_fl)).Bar(),
			  Vec4D(mpi_p.m_mom[0],mpi_p.m_mom[1],mpi_p.m_mom[2],mpi_p.m_mom[3]));
  if (mpi_p.m_flow[0]) p->SetFlow(1,mpi_p.m_flow[0]);
  if (mpi_p.m_flow[1]) p->SetFlow(2,mpi_p.m_flow[1]);
  return p;
}






