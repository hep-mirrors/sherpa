#include "Parton.H"
#include "Blob.H"
#include "Random.H"
#include "Run_Parameter.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;

std::ostream& APHYTOOLS::operator<<(std::ostream& str,Parton* part) {
  switch (part->status()) {
  case 0 : // null entry
    return str<<"--- empty entry ---"<<std::endl;
  case 1 : // active (final state) particle
    str<<part->info();
    str.width(8);
    str<<part->flav();
    str<<' ';
    str<<part->Get_Numb();
    str<<' ';
    if (part->prod()) str<<" produced in "<<part->prod()->Id()<<std::endl<<" ";
    if (part->dec())  str<<" decayed in "<<part->dec()->Id()<<std::endl<<" ";
    break;
  case 2 : // decayed or fragmented
    str<<part->info();
    str<<'!';
    str.width(7);
    str<<part->flav();
    str<<'!';
    str<<part->Get_Numb();
    str<<' ';
    if (part->prod()) str<<" produced in "<<part->prod()->Id()<<std::endl<<" ";
    if (part->dec())  str<<" decayed in "<<part->dec()->Id()<<std::endl<<" ";
    break;
  case 3 : // documentation line
    return str<<"============================================================"<<std::endl;
  default : // user defined or reserved
    return str<<"--- unrecognized status:"<<part->status()<<" ---"<<std::endl;
  }
  str<<"      "<<part->momentum()<<" "<<part->momentum().abs2()
     <<" colours : ("<<part->flow(1)<<","<<part->flow(2)<<")"<<std::endl;
  return str;
}

Parton::~Parton() {
  //  cout<<" delete "<<this<<endl;
}

Parton::Parton() {
  Number     = -1;
  m_info     = ' ';
  m_status   = 1;
  m_fl       = Flavour(kf::none);
  m_momentum = AMATOOLS::vec4d(0,0,0,0); 
  startblob  = 0;
  endblob    = 0;
  T_Dec      = 0.;
};

Parton::Parton(Parton * in)  {
  Number     = in->Number;
  m_info     = in->m_info;
  m_status   = in->status();
  m_fl       = in->m_fl;
  m_momentum = in->m_momentum;
  startblob  = in->startblob;
  endblob    = in->endblob;
  T_Dec      = in->T_Dec;
  m_flow     = in->m_flow;
}

Parton::Parton(int number,Flavour fl,AMATOOLS::vec4d p)  {
  Number     = number;
  m_status   = 1;
  m_info     = ' ';
  m_fl       = fl;
  m_momentum = p;
  startblob  = 0;
  endblob    = 0;
  T_Dec      = 0.;
}


void Parton::Copy(Parton* in)  {
  Number     = in->Number;
  m_info     = in->m_info;
  m_status   = in->status();
  m_fl       = in->m_fl;
  m_momentum = in->m_momentum;
  startblob  = in->startblob;
  endblob    = in->endblob;
  T_Dec      = in->T_Dec;
  m_flow     = in->m_flow;
}

double Parton::ProperTime() {
  double q2    = m_momentum.abs2();
  double m2    = AMATOOLS::sqr(m_fl.mass());
  double tau2  = 1.e6;
  if (!( (q2-m2 < rpa.gen.Accu()) && (m_fl.width() < rpa.gen.Accu()))) { 
    if (m2>AORGTOOLS::rpa.gen.Accu()) { 
      tau2 = q2/(AMATOOLS::sqr(q2-m2)+AMATOOLS::sqr(q2*m_fl.width())/m2);
    }
    else {
      if (dabs(q2)>AORGTOOLS::rpa.gen.Accu()) tau2 = 1/dabs(q2);
    }
  }
  else {
    if (m_fl.strong()) tau2 = 1./sqr(0.2); 
  }
  return AORGTOOLS::rpa.hbar() * sqrt(tau2);
};

double Parton::Lifetime() {
  double t   = -ProperTime()*log(1.-AMATOOLS::Ran.get());  
  if (t>1.e6) t = 1.e6;
  double gamma = 1./AORGTOOLS::rpa.gen.Accu();
  if (m_fl.mass()>AORGTOOLS::rpa.gen.Accu()) gamma = E()/m_fl.mass();
  return gamma * t;      
};

AMATOOLS::vec3d Parton::Distance() {
  double gamma     = E()/m_fl.mass();
  AMATOOLS::vec3d v = AMATOOLS::vec3d(m_momentum)/E()*AORGTOOLS::rpa.c();
  return v*Lifetime();
};

  // Numbers etc.
int    Parton::Get_Numb()   const                          { return Number; }
void   Parton::Set_Numb(const int N)                       { Number    = N; }
int    Parton::Get_JetN()   const                          { return JetNumber; }
void   Parton::Set_JetN(const int N)                       { JetNumber = N; }

  // Status etc.
int    Parton::status() const                              { return m_status; }
void   Parton::set_status( int status )                    { m_status = status; }
char   Parton::info() const                                { return m_info;}
void   Parton::set_info(char info)                         { m_info = info; }

  // Momentum, energy, and lifetime
AMATOOLS::vec4d Parton::momentum() const                   { return m_momentum; }
double Parton::E()                                         { return m_momentum[0];}
void   Parton::set_momentum(const AMATOOLS::vec4d & vec4 ) { m_momentum = vec4; } 
double Parton::time() const                                { return T_Dec; }
void   Parton::set_time(const int t)                       { T_Dec     = t; }
void   Parton::set_time()                                  { T_Dec     = Lifetime(); }

  // Production and decay vertices
AMATOOLS::vec4d Parton::xprod()                            { return startblob->Position(); }
Blob * Parton::prod()                                      { return startblob; }
void   Parton::set_prod(Blob * _blob)                      { startblob = _blob; }
AMATOOLS::vec4d Parton::xdec()                             { return endblob->Position(); }
Blob * Parton::dec()                                       { return endblob; }
void   Parton::set_dec(Blob * _blob)                       { endblob = _blob; }

  // Flavour and flow
APHYTOOLS::Flavour Parton::flav() const                    { return m_fl; }
void   Parton::set_flav(APHYTOOLS::Flavour & fl)           { m_fl      = fl; }
APHYTOOLS::Flow Parton::flow() const                       { return m_flow; }
int    Parton::flow(const int index) const                 { return m_flow.icode(index); }
void   Parton::set_flow(const APHYTOOLS::Flow & f )        { m_flow    = f; }
void   Parton::set_flow(const int index, const int code)   {
  if (!m_fl.strong()) return;
  if ( code==-1 ) m_flow.set_unique_icode(index);
  else m_flow.set_icode(index,code);
};



void APHYTOOLS::Parton2MPI(const Parton * p , MPI_Parton & mpi_p) {
  mpi_p.id  =p->Get_Numb();
  mpi_p.m_fl=int(p->flav());
  for (int i=0; i<4; ++i)  
    mpi_p.m_mom[i]=p->momentum()[i];
  mpi_p.m_flow[0]=p->flow(1);
  mpi_p.m_flow[1]=p->flow(2);
}
  
Parton * APHYTOOLS::MPI2Parton(const MPI_Parton & mpi_p ) {
  Parton * p ;
  if (mpi_p.m_fl>0) p= new Parton(mpi_p.id, Flavour((kf::code)mpi_p.m_fl),
			  vec4d(mpi_p.m_mom[0],mpi_p.m_mom[1],mpi_p.m_mom[2],mpi_p.m_mom[3]));
  else p= new Parton(mpi_p.id, Flavour((kf::code)(-mpi_p.m_fl)).bar(),
			  vec4d(mpi_p.m_mom[0],mpi_p.m_mom[1],mpi_p.m_mom[2],mpi_p.m_mom[3]));
  if (mpi_p.m_flow[0]) p->set_flow(1,mpi_p.m_flow[0]);
  if (mpi_p.m_flow[1]) p->set_flow(2,mpi_p.m_flow[1]);
  return p;
}






