#include "Hadrons.H"
#include "Data_Reader.H"
#include "Flavour.H"
#include "Message.H"
#include "Random.H"
#include <time.h>
#include "Hadron_Decay_Channel.H"
#include <fstream>
#include "Run_Parameter.H"
#include "MyStrStream.H"
#include <utility>
#include <algorithm>
#include "XYZFuncs.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Hadrons::Hadrons( string _path, string _file, string _constfile ) : 
  m_path(_path), m_file(_file), m_constfile(_constfile), m_createbooklet(0)
{ 
  msg_Tracking()<<"In Hadrons: ("<<_path<<") "<<_file<<std::endl;
  msg_Tracking()<<"In Hadrons: ("<<_path<<") "<<_constfile<<std::endl;
  ReadInConstants();
  ReadInDecayTables();
//  CreateBookletNow();
}


void Hadrons::CreateBookletNow()
{
  string fn = string("DecayBooklet");
  string texfn = m_path + fn + string(".tex");
  ofstream f(texfn.c_str());
   
  // header
  f<<"\\documentclass{article}\n"
   <<"\\usepackage{latexsym,amssymb,amsmath,amsxtra}\n\n"
   <<"\\begin{document}\n\n"; 
  f<<"\\newcommand{\\m}{-}\n"; 
  f<<"\\newcommand{\\p}{+}\n"; 
  f<<"\\title{Decay Information}\n\\maketitle\n";

  // text 
  for ( map<kf::code,Decay_Table *>::iterator pos = p_decaymap->begin(); pos != p_decaymap->end(); ++pos) {
	kf::code      kfc (pos->first);
	Decay_Table * dt  (pos->second);
	f<<"\\section{Decaying Particle: "<<dt->Flav()<<"}\n";
	f<<"\\begin{center} \n\\begin{tabular}{ll}\n";
	f<<" number of decay channels:	& "<<dt->NumberOfDecayChannels()<<"\\\\ \n";
	f<<" total width:               & "<<dt->TotalWidth()<<"\\\\ \n";
	f<<"\\end{tabular} \n\\end{center}\n";
  }
   
  // end 
  f<<"\\end{document}\n"; 
}

bool Hadrons::FindDecay(const ATOOLS::Flavour & Decayer)
{
  if (p_decaymap->find(Decayer.Kfcode())==p_decaymap->end()) return false;
  p_table = (*p_decaymap)[Decayer.Kfcode()];
  return true;
}

Hadron_Decay_Channel * Hadrons::ChooseDecayChannel()
{
  const int nchan = p_table->NumberOfDecayChannels();
  Decay_Channel * dec_channel;
  if (p_table->Flav().Kfcode() != kf::K ) {
    double TotalWidth = p_table->TotalWidth();
    bool channel_chosen (0);
    int k (0);

    // dice decay channel acc. to BR
    if( nchan>1 ) {
      while( !channel_chosen ) {
        k = int( ran.Get() * nchan );					// dice decay channel
        double r = ran.Get();							// random number for rejection
        if( r < p_table->Width(k) / TotalWidth ) {
          dec_channel = p_table->GetDecayChannel(k);
          channel_chosen = 1;
        }
      }	
    }
    else dec_channel = p_table->GetDecayChannel(0);
  }
  else {                                            // treatment for K0's
    if (nchan!=2) {
      msg.Error()<<"ERROR in Hadrons::ChooseDecayChannel() : "<<endl
                 <<"     K0 must have two channels K(L) and K(S), but is does not."<<endl
                 <<"     Don't know what to do, will abort."<<endl;
      abort();
    }
    if (ran.Get() < 0.5) dec_channel = p_table->GetDecayChannel(0); // 50% K(L)
    else                 dec_channel = p_table->GetDecayChannel(1); // 50% K(S)
  }

  if( p_channelmap->find(dec_channel) == p_channelmap->end() ) {
	msg.Error()<<"Error in Hadrons::ChooseDecayChannel() \n"
	  		   <<"      Couldn't find appropriate channel pointer for "<<dec_channel->GetDecaying()
			   <<" decay. \n"
			   <<"      Don't know what to do, will abort."<<endl;
    abort();
  }
  return (*p_channelmap)[dec_channel];
}

// implementation with spin correlation
// hep-ph/0110108

void Hadrons::ChooseDecayKinematics( 
    Vec4D                   * _p, 
    Hadron_Decay_Channel    * _hdc, 
    Spin_Density_Matrix     * sigma )
{
  double value(0.);
  double over_factor (1.);                                  // overestimate maximum factor
  if( sigma ) over_factor = double((*_hdc->Flavours())[0].IntSpin() + 1);
  const double max = _hdc->GetPS()->Maximum() * over_factor;    // note: no flux in max.
//  int trials(0);											// number of trials
  do {
    value = _hdc->Differential(_p,sigma);					// current val. of T
//    trials++;
  } while( ran.Get() > value/max );
//  ofstream f("trials.out",ios::app|ios::out );
//  f<<trials<<endl;
}

Spin_Density_Matrix Hadrons::PerformDecay( 
    Particle            * part, 
    Blob_List           * blob_list, 
    Particle_List       * part_list,
    Spin_Density_Matrix * sigma )
{
  msg_Tracking()<<"Hadrons::PerformDecay() for "<<part->Flav()<<" with spin correlation @ "<<sigma<<endl;
  msg.Debugging()<<"Momentum: "<<part->Momentum()<<endl;
  if( part->Flav().Kfcode() != p_table->Flav().Kfcode() ) {
	msg.Error()<<"ERROR in Hadrons::PerformDecay() : Particle in selected decay table is not \n"
	           <<"        identical to the decayer.\n"
			   <<"        >>>   "<<p_table->Flav()<<" <-> "<<part->Flav()<<"   <<<"<<endl
			   <<"        Don't know which to take, will abort."<<endl;
	abort();
  }

  bool       anti = part->Flav().IsAnti();			    // antiparticle ?
  int mother_spin = part->Flav().IntSpin();             // 2*spin

  // choose decay channel acc. to BR
  Hadron_Decay_Channel * hdc = ChooseDecayChannel();
  FlavourSet       daughters = hdc->DecayChannel()->GetDecayProducts();
  const int                n = hdc->DecayChannel()->NumberOfDecayProducts()+1;
  if( (anti) && !(hdc->FullDecay() & 1) ) {             // if NO_ANTI
    if( sigma ) {
      Spin_Density_Matrix unitmatrix = Spin_Density_Matrix(mother_spin);
      unitmatrix.SetUnitMatrix();
      return unitmatrix;                                // return unit matrix
    }
    return Spin_Density_Matrix();                       // return empty decay matrix
  }
   
  // choose a kinematics that corresponds to the ME kinematic distribution
  Vec4D mom[n];
  mom[0] = part->Momentum();
  if( n>2 ) ChooseDecayKinematics( mom, hdc, sigma ); 
  else      {
    mom[1] = mom[0];
  }

  // Create blob
  Blob * blob;											// decay blob
  blob = new Blob();
  blob->SetStatus(1);
  blob->SetType( btp::Hadron_Decay );
  blob->SetTypeSpec("Sherpa");
  blob->SetId();
  double        time = part->LifeTime();
  Vec3D      spatial = part->Distance( time );
  Vec4D     position = Vec4D( time*rpa.c(), spatial );
  blob->SetPosition( part->XProd() + position );
  msg.Debugging()<<"created new blob: #"<<blob->Id()<<" with status 1"<<endl;
  blob->AddToInParticles( part );
  if( part->Info() == 'P' ) part->SetInfo('p');
  if( part->Info() == 'D' ) part->SetInfo('d');
  blob_list->push_back( blob );
   
  Particle * particle;									// daughter part.
  Vec4D momentum;										// daughter mom.
  Flavour flav;											// daughter flav.
  FlSetConstIter dit;									// Iterator
  int i(0);
  vector<pair<Flavour,int> > daughters_pair;
  for( dit = daughters.begin(); dit != daughters.end(); dit++ ) {
    i++;
    daughters_pair.push_back( pair<Flavour,int>( (*dit), i ) );
  }

  // create new spin correlation tensor and shuffle randomly daughters
  Spin_Correlation_Tensor * SCT (NULL);
  if( Spin_Correlation_Tensor::Mode()!=scmode::None && hdc->GetIndexList()->size() ) {
    SCT = 
      new Spin_Correlation_Tensor( hdc->GetIndexList(), hdc->GetAmplitudeTensor() );
    random_shuffle( daughters_pair.begin(), daughters_pair.end() );
  }

  // treat every daughter
  vector<pair<Flavour,int> >::iterator dpit;
  for( dpit = daughters_pair.begin(); dpit != daughters_pair.end(); dpit++ ) {
    i = dpit->second;
	momentum = mom[i];
	flav = dpit->first;
	if( anti ) flav = flav.Bar();
	msg.Debugging()<<"Daughter: "<<flav<<"  "<<momentum<<endl;
	particle = new Particle( -1, flav, momentum );
	if( part_list ) particle->SetNumber( part_list->size() );
	else particle->SetNumber( 0 );
	particle->SetStatus(1);
	particle->SetInfo('D');
	if( part_list ) part_list->push_back( particle ); 
	blob->AddToOutParticles( particle );
	// check if daughter can be treated as well
    int spin = particle->Flav().IntSpin();    // 2*spin
	if ( ( hdc->FullDecay() & 2 ) &&
         ( !(particle->Flav().IsStable()) ) &&
         ( FindDecay(particle->RefFlav()) )  
       )  
    {                          // it is not stable
      if( spin && SCT ) {
        Spin_Density_Matrix * new_sigma    (NULL);
        Spin_Density_Matrix * decay_matrix (NULL);
        if( sigma )             // if mother has SDM
          new_sigma    = new Spin_Density_Matrix( SCT->GetSigma(i,sigma) );
        else                   
          new_sigma    = new Spin_Density_Matrix( SCT->GetSigma(i) );
        new_sigma->Normalise();
        decay_matrix   = new Spin_Density_Matrix();
        (*decay_matrix) = PerformDecay( particle, blob_list, part_list, new_sigma );
        decay_matrix->Normalise();
        // contract over decay matrix to reduce the size of the SCT
        SCT->Contract( i, decay_matrix );
        delete new_sigma;
        delete decay_matrix;
      }
      else PerformDecay( particle, blob_list, part_list, NULL );
    }
    else {                      // it is stable
      if( spin && SCT ) SCT->Contract( i, NULL );   // contract with unitmatrix
    }
  }
  // finished with all daughters
  if( SCT ) {
    if( sigma ) {
      Spin_Density_Matrix decm = SCT->GetSigma(0);
      delete SCT;
      return decm;                              // return decay matrix for this particle
    }
    delete SCT;
    return Spin_Density_Matrix(0);              // empty Decay Matrix
  }
  return Spin_Density_Matrix(0);              // empty Decay Matrix
}

void Hadrons::ReadInConstants()
{
  m_md0.clear();                            // clear model 
  Data_Reader reader = Data_Reader(string("|"),string(";"),string("!"));
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_constfile);

  vector<vector<string> > constants;
  reader.SetMatrixType(mtc::transposed);
  if(!reader.MatrixFromFile(constants)) {
    msg.Out()<<"Warning! The file "<<m_path<<m_constfile<<" does not exist"<<endl
             <<"     or has some syntax error."<<endl;
    msg.Out()<<"     Will ignore it and hope for the best."<<endl;     
    return;
  }

  for (size_t i=0;i<constants.size();++i) {
    if( constants[i][1] == "=" ) {              // <name> = <value>
      m_md0[constants[i][0]] = ToType<double> (
          reader.Interpreter()->Interprete(constants[i][2]) );
    }
  }
}
 
void Hadrons::ReadInDecayTables()
{
  Data_Reader reader = Data_Reader(string("->"),string(";"),string("!"));
  reader.SetAddCommandLine(false);
  reader.AddComment("#");
  reader.AddComment("//");
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_file);

  vector<vector<string> > Decayers;
  reader.SetMatrixType(mtc::transposed);
  if(!reader.MatrixFromFile(Decayers)) {
    msg.Error()<<"ERROR in Hadrons::ReadInDecayTables() :\n"
	       <<"   Read in failure "<<m_path<<m_file<<", will abort."<<endl;
    abort();
  }

  p_decaymap = new map<ATOOLS::kf::code,ATOOLS::Decay_Table *>;
  p_channelmap = new map< Decay_Channel*, Hadron_Decay_Channel* >;
  Decay_Table * dt;
  Flavour fl;
  for (size_t i=0;i<Decayers.size();++i) {
    fl = Flavour(kf::code(atoi((Decayers[i][0]).c_str())));
    if (p_decaymap->find(fl.Kfcode())!=p_decaymap->end()) {
      msg.Error()<<"ERROR in Hadrons::ReadInDecayTables() :"<<endl
		 <<"   Flavour "<<fl
		 <<" already in map. Don't know what to do, will abort."<<endl;
      abort();
    }
    dt = InitialiseOneDecayTable(Decayers[i]);
    (*p_decaymap)[fl.Kfcode()] = dt;
  }
}

// interprete one line of HadronDecays.dat
// line: kfcode -> filepath/  filename
Decay_Table * Hadrons::InitialiseOneDecayTable(vector<string> line)
{
  Decay_Table * dt              = new Decay_Table(Flavour(kf::code(atoi((line[0]).c_str()))));
  string lcpath (line[1]);      // path of decay files
  Decay_Table_Reader * dtreader = new Decay_Table_Reader(m_path+lcpath,line[2]);
  if (dtreader->FillDecayTable(dt)>0) {     // if at least one channel defined
    msg.Out()<<om::blue<<"Found "<<dt->NumberOfDecayChannels()<<" decay channels for "<<dt->Flav()<<om::reset<<endl;
    dtreader->FillInMatrixElementsAndPS(dt,p_channelmap,m_md0);
    if( msg.LevelIsInfo() ) {
      msg.Info()<<"Initialised a new decay table : "<<endl;
      msg.Info()<<"(Using information given in .dat files, such as BR, width, ...)"<<endl;
      dt->Output();
      msg.Info()<<"Calculated decay widths using implemented models :"<<endl
        <<"(only for information; they are NOT used for the simulation)"<<endl;
      for (int i=0;i<dt->NumberOfDecayChannels();i++) {			
        msg.Info()
          <<(*p_channelmap)[dt->GetDecayChannel(i)]->ChannelName()<<" : "
          <<(*p_channelmap)[dt->GetDecayChannel(i)]->GetPS()->Result()<<" ("
          <<(*p_channelmap)[dt->GetDecayChannel(i)]->GetPS()->RelError()<<" %)"
          <<" GeV"<<endl;
      }
      msg.Info()<<endl;
    }
  }
  else { 
    msg.Error()<<"WARNING in Hadrons::InitialiseOneDecayTable : "<<endl
	       <<"   No decay channels found for "<<dt->Flav()<<" in file "<<line[1]<<endl
	       <<"   Will continue and hope for the best."<<endl;
    delete dt; 
    dt = NULL;
  }
  delete dtreader;
  return dt;
}
