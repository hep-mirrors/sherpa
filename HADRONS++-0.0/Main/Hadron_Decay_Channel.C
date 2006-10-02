#include "Hadron_Decay_Channel.H"
#include "HD_ME_Selector.H"
#include "Spin_Correlation_Tensor.H"
#include "Poincare.H"
#include "Random.H"
#include <stdio.h>

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Hadron_Decay_Channel::Hadron_Decay_Channel( Decay_Channel * _dc, string _path ) :
  Integrable_Base(1,_dc->NumberOfDecayProducts()),
  p_dc(_dc), p_me(NULL), m_metype(p_dc->METype()),
  m_path(_path), m_fulldecay(3), m_createbooklet(0),
  p_amplitudes(NULL)
{
  m_resultpath      = string("./");                             // where to write results
  m_resultfile      = m_histofile = string("");                 // filename
  m_name            = p_dc->ProcessName();                      // name of proces
  p_flavours        = new Flavour[m_nin+m_nout];                    
  p_flavours[0]     = p_dc->GetDecaying();                      // decaying particle
  int counter = ATOOLS::Particle::Counter();
  m_particles.push_back(new Particle(-1,p_flavours[0]));
  m_channelname     = string("");
  m_chnamenumbers   = string("0");
  m_channelname     = p_flavours[0].IDName() + string(" --> ");
  m_chnamenumbers.append( (p_flavours[0].IDName()).length()-1, ' ' );
  m_chnamenumbers.append(" --> ");
  char helpch[2];
  for (int i=0;i<p_dc->NumberOfDecayProducts();i++) {           // decay products
    p_flavours[i+1] = p_dc->GetDecayProduct(i);
    m_particles.push_back(new Particle(-1,p_flavours[i+1]));
    m_channelname  += p_flavours[i+1].IDName() + string(" ");
    sprintf( helpch, "%i", i+1 );
    m_chnamenumbers.append(string(helpch));
    m_chnamenumbers.append( (p_flavours[i+1].IDName()).length(), ' ' );
  }
  ATOOLS::Particle::ResetCounter(counter);
  p_amplitudes = new Amplitude_Tensor(m_particles);
  HD_ME_Selector mesel;                                         // ME selector
  p_me = mesel.GetME(m_nin,m_nout,p_flavours);                  // get the appropr. ME
  p_me->SetPath( m_path );                                      // set Decaydata path 
  msg.Tracking()<<"Matrix Element for "<<m_channelname<<" : "<<p_me->METype()<<"."<<endl;
}

Hadron_Decay_Channel::~Hadron_Decay_Channel()
{
  if (p_dc) { delete p_dc; p_dc=NULL; }
  if (p_ps) { delete p_ps; p_ps=NULL; }
  if (p_me) { delete p_me; p_me=NULL; }
  if(p_amplitudes) { delete p_amplitudes; p_amplitudes=NULL; }
  for(size_t i=0;i<m_particles.size();i++) {
    if(m_particles[i]) { delete m_particles[i]; }
  }
}


bool Hadron_Decay_Channel::InitialisePhaseSpace(vector<string> & PStype, GeneralModel startmd) 
  // PStype: decay products | BR | DC filename  <-- line of Decay file
{
  bool mustinit;
  struct GeneralModel locmd (startmd);
  bool rewriteH (false);
  if ( PStype.size() == 2 ) {                                   // in case no DC file given
    string fn("");                                              // filename of DC file
    fn += p_dc->GetDecaying().ShellName() + string("_");
    for ( int i=0; i<p_dc->NumberOfDecayProducts(); i++ ) {
      fn += p_dc->GetDecayProduct(i).ShellName();
    }
    fn += string(".dat");
    PStype.push_back( fn );                                     // generate DC filename
    rewriteH = true;                                            // rewrite decay table file
  }
  // check if dc file exists
  string dcfilename = PStype[2];
  ifstream dcf( (m_path+dcfilename).c_str() );
  bool read_dc = dcf;                               // read DC file if it exists
  p_ps = new HD_PS_Base(this,m_path,PStype,mustinit,locmd,read_dc); // new PS
  WriteModelOnScreen(locmd);                        // show model parameters on screen
  p_me->SetModelParameters( locmd );                // set parameters for ME
  if (mustinit) p_ps->Initialise();                 // if in need => ini
  return rewriteH;                                  // rewrite decay table file ?
}

void Hadron_Decay_Channel::WriteModelOnScreen( GeneralModel _locmd )
{
  msg.Debugging()
    <<"-----------------------------------------------------\n"
    <<"Modelparameters for channel: "<<m_channelname<<"\n"<<endl;
  GeneralModel::iterator md_it;
  for ( md_it = _locmd.begin(); md_it != _locmd.end(); ++md_it ) {
    msg.Debugging()<<"   "<<md_it->first<<":\t"<<md_it->second<<endl;
  }
  msg.Debugging()
    <<"-----------------------------------------------------"<<endl;
}


// differential with random PS points; just for weight
double Hadron_Decay_Channel::Differential()
{
  p_momenta[0] = Vec4D(p_flavours[0].PSMass(),0.,0.,0.);        // decay from rest
  p_ps->GeneratePoint(p_momenta);                               // generate a PS point
  p_ps->GenerateWeight(p_momenta);                              // calculate its weight factor
  double weight = p_ps->Weight();                               // get weight factor
  CalculateAmplitudes(p_momenta,p_amplitudes,false);
  double value=p_amplitudes->SumSquare();
  value /= (p_dc->GetDecaying().IntSpin()+1);
  return value*weight;
}


// differential with incoming (possibly offshell) momentum; for weight and momenta
double Hadron_Decay_Channel::Differential( Vec4D * mom, bool anti )
{
#ifdef DEBUG__Hadrons
  if( !IsZero(mom[0][1]) || !IsZero(mom[0][2]) || !IsZero(mom[0][3]) ) {
    PRINT_INFO("Error: given momentum is not in CMS: "<<mom[0]);
  }
#endif
  p_ps->GeneratePoint(mom);                               // generate a PS point
  p_ps->GenerateWeight(mom);                              // calculate its weight factor
  
  double weight = p_ps->Weight();                               // weight factor
  CalculateAmplitudes(mom,p_amplitudes, anti);
  double value=p_amplitudes->SumSquare();
  value /= (p_dc->GetDecaying().IntSpin()+1);
  return value*weight;
}

void Hadron_Decay_Channel::CalculateAmplitudes( Vec4D * moms, Amplitude_Tensor* amps, bool anti )
{
  p_me->SetAnti(anti);
  (*p_me)(  moms,                                               // phase space point
            amps,                                               // ampl. tensor and indices
            1);
}

namespace ATOOLS {
  template <> Blob_Data<HADRONS::Hadron_Decay_Channel*>::~Blob_Data() { }
  template class Blob_Data<HADRONS::Hadron_Decay_Channel*>;
  template HADRONS::Hadron_Decay_Channel* &Blob_Data_Base::Get<HADRONS::Hadron_Decay_Channel*>();
}

