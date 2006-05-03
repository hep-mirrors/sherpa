#include "Hadron_Decay_Channel.H"
#include "HD_ME_Selector.H"
#include "Spin_Correlation_Tensor.H"
#include "Poincare.H"
#include "Random.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Hadron_Decay_Channel::Hadron_Decay_Channel( Decay_Channel * _dc, string _path ) :
  p_dc(_dc), p_me(NULL), m_metype(p_dc->METype()),
  Integrable_Base(1,_dc->NumberOfDecayProducts()),
  m_path(_path), m_fulldecay(3), m_mass_smearing(true), m_createbooklet(0),
  p_indices(NULL), p_ampls(NULL)
{
  m_resultpath      = string("./");                             // where to write results
  m_resultfile      = m_histofile = string("");                 // filename
  m_name            = p_dc->ProcessName();                      // name of proces
  p_flavours        = new Flavour[m_nin+m_nout];                    
  p_flavours[0]     = p_dc->GetDecaying();                      // decaying particle
  m_channelname     = string("");
  m_chnamenumbers   = string("0");
  m_channelname     = p_flavours[0].IDName() + string(" --> ");
  m_chnamenumbers.append( (p_flavours[0].IDName()).length()-1, ' ' );
  m_chnamenumbers.append(" --> ");
  char helpch[2];
  for (int i=0;i<m_nout;i++) {           // decay products
    p_flavours[i+1] = p_dc->GetDecayProduct(i);         
    m_channelname  += p_flavours[i+1].IDName() + string(" ");
    sprintf( helpch, "%i%", i+1 );
    m_chnamenumbers.append(string(helpch));
    m_chnamenumbers.append( (p_flavours[i+1].IDName()).length(), ' ' );
  }
  p_indices = new vector<pair<int,int> >;                        // index bookkeeping
  p_ampls   = new vector<Complex>;                              // new amplitude tensor
  HD_ME_Selector mesel;                                         // ME selector
  p_me = mesel.GetME(m_nin,m_nout,p_flavours);                  // get the appropr. ME
  p_me->SetPath( m_path );                                      // set Decaydata path 
  msg.Out()<<"Matrix Element for "<<m_channelname<<" : "<<p_me->METype()<<"."<<endl;
  // check for identical particles
  Flavour refflav;
  double symfactor (1);         
  int l(0), lfac (1);              
  for( int i=0; i<m_nout; ++i ) {
    refflav = p_flavours[i+1];
    l = 0;
    lfac = 1;
    for( int j=0; j<m_nout; ++j ) {
      if( p_flavours[j+1]==refflav ) {
        l++;
        lfac *= l;
      }
    }
    symfactor *= lfac;
  }
  m_symmetry = 1./sqrt(symfactor);
}

Hadron_Decay_Channel::~Hadron_Decay_Channel()
{
  if (p_dc) { delete p_dc; p_dc=NULL; }
  if (p_ps) { delete p_ps; p_ps=NULL; }
  if (p_me) { delete p_me; p_me=NULL; }
  if (p_ampls) { delete p_ampls; p_ampls=NULL; }
  if (p_indices) { delete p_indices; p_indices=NULL; }
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
    rewriteH = true;                                            // rewrite hadron decay file
  }
  // check if dc file exists
  string dcfilename = PStype[2];
  ifstream dcf( (m_path+dcfilename).c_str() );
  bool read_dc = dcf;                               // read DC file if it exists
  p_ps = new HD_PS_Base(this,m_path,PStype,mustinit,locmd,read_dc); // new PS
  WriteModelOnScreen(locmd);                        // show model parameters on screen
  p_me->SetModelParameters( locmd );                // set parameters for ME
  if (mustinit) p_ps->Initialise();                 // if in need => ini
  return rewriteH;                                  // rewrite H file ?
}

void Hadron_Decay_Channel::WriteModelOnScreen( GeneralModel _locmd )
{
  if( msg.LevelIsDebugging() ) {
  msg.Out()
    <<"-----------------------------------------------------\n"
    <<"Modelparameters for channel: "<<m_channelname<<"\n"<<endl;
  GeneralModel::iterator md_it;
  for ( md_it = _locmd.begin(); md_it != _locmd.end(); ++md_it ) {
    msg.Out()<<"   "<<md_it->first<<":\t"<<md_it->second<<endl;
  }
  msg.Out()
    <<"-----------------------------------------------------"<<endl;
  }
}

// w/o spin correlation

double Hadron_Decay_Channel::Differential()
{
  p_momenta[0] = Vec4D(p_flavours[0].PSMass(),0.,0.,0.);   // decay from rest
  p_ps->GeneratePoint(p_momenta);                               // generate a PS point
  p_ps->GenerateWeight(p_momenta);                              // calculate its weight factor
  double weight = p_ps->Weight();                               // get weight factor
  // get amplitude tensor
  (*p_me)(  p_momenta,                                          // phase space point
            p_ampls, p_indices,                                 // ampl. tensor and indices       
            1 );                                            // spinor base
  double value (0.);
  if( p_ampls->size() ) {
    for( size_t i=0; i<p_ampls->size(); ++i ) {
      value += norm( (*p_ampls)[i] );
    }
    value /= (p_dc->GetDecaying().IntSpin()+1);
  }
  else value = 1.;                                              // isotropic
  return value*weight*m_symmetry;
}


// with spin correlation
double Hadron_Decay_Channel::Differential( Vec4D * mom, Spin_Density_Matrix * sigma )
{
  if( !mom ) return Differential();                             // if no momentum
  double ret;
  if( !sigma ) {                                                // no SDM
    ret = Differential();
    // boost into Lab frame
    Poincare lambda(mom[0]);
    lambda.Invert();
    for( int i=0; i<m_nout+1; ++i ) {
      p_momenta[i] = lambda*p_momenta[i];
      mom[i] = p_momenta[i];
    }
    return ret;
  }
  // spin correlations 
  // get PS point in rest frame
  p_momenta[0] = Vec4D(p_flavours[0].PSMass(),0.,0.,0.);        // decay from rest
  p_ps->GeneratePoint(p_momenta);                               // generate a PS point
  p_ps->GenerateWeight(p_momenta);                              // calculate its weight factor
  double weight = p_ps->Weight();                               // weight factor
  // boost into Lab system 
  Poincare lambda(mom[0]);
  lambda.Invert();
  for( int i=0; i<m_nout+1; ++i ) {
    p_momenta[i] = lambda*p_momenta[i];
    mom[i] = p_momenta[i];
  }
  // get amplitude tensor
  (*p_me)(  p_momenta,                                          // phase space point
            p_ampls, p_indices,                                 // ampl. tensor and indices
            Spin_Correlation_Tensor::Get_k0_n() );              // spinor base
  double value;
  if( p_ampls->size()==0 ) {                                    // no ampls <=> isotropic
    CreateTrivial(sigma);                                       // create trivial amplitude tensor
    value = 1.;
  }
  else {                                                        // not isotropic
    if( p_indices->size() ) {                                   // if there are indices
      Spin_Correlation_Tensor help_sct ( p_indices, p_ampls );  // temporary SCT
      help_sct.Contract(0,sigma);                               // contract over mother sigma
      value = real( help_sct.Trace() );                         // get T by taking trace of rest
    }
    else {                                                      // no spin correlation
      value = norm( (*p_ampls)[0] );                            // only one amplitude
    }
  }
  return value*weight*m_symmetry;
}

void Hadron_Decay_Channel::CreateTrivial( Spin_Density_Matrix * sigma )
{
  // create trivial amplitude tensor and its index bookkeeping
  p_indices->clear();
  p_ampls->clear();
  if( sigma ) {                                               // if spin > 0
    p_indices->push_back(pair<int,int>(0,sigma->Spin()));     // order of ind. does not matter
    for( int j=0; j<sigma->NrEntries(); ++j ) 
      p_ampls->push_back(1.);                                 // M=1 <=> isotropic
  }
  int spin;
  for( int i=0; i<m_nout; ++i ) {
    spin = p_flavours[i+1].IntSpin();                         // 2*spin of daughter
    if( spin ) {                                              // if spin > 0
      p_indices->push_back(pair<int,int>(i+1,spin));          // order of ind. does not matter
      for( int j=0; j<sqr(spin+1); ++j ) 
        p_ampls->push_back(1.);                               // M=1 <=> isotropic
    }
  }
}
