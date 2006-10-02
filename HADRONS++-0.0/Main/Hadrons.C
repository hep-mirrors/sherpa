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
#include "Decay_Table.H"
#include "Vector.H"
#include "Particle.H"
#include "Particle_List.H"
#include "Spin_Density_Matrix.H"
#include "Spin_Correlation_Tensor.H"
#include "Blob_List.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Hadrons::Hadrons( string _path, string _file, string _constfile ) : 
  m_path(_path), m_file(_file), m_constfile(_constfile), m_sc(false), m_createbooklet(0)
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

Hadron_Decay_Channel * Hadrons::ChooseDecayChannel(Blob* blob)
{
  Blob_Data_Base* data = (*blob)["hdc"];
  if(data) return data->Get<Hadron_Decay_Channel*>();
  
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
  blob->AddData("hdc",new Blob_Data<Hadron_Decay_Channel*>((*p_channelmap)[dec_channel]));
  return (*p_channelmap)[dec_channel];
}

// implementation with spin correlation
// hep-ph/0110108

void Hadrons::ChooseDecayKinematics( 
    vector<Vec4D>           & moms,
    Hadron_Decay_Channel    * hdc,
    bool                      anti )
{
  double value(0.);
  const double max = hdc->GetPS()->Maximum();    // note: no flux in max.
  do {
    value = hdc->Differential(&moms.front(),anti);		// current val. of T
  } while( ran.Get() > value/max );
}


Return_Value::code Hadrons::PerformDecay( Blob* blob, const Vec4D& labmom )
{
  Particle*  inpart    = blob->InParticle(0);
#ifdef DEBUG__Hadrons
  if(blob->NOutP()>0) {
    msg.Error()<<METHOD<<" Blob ["<<blob->Id()<<"] in event "<<rpa.gen.NumberOfDicedEvents()
      <<" already has outparticles: "<<endl
      <<(*blob)<<endl;
  }
#endif
  if (!blob->Has(blob_status::needs_hadrondecays) ||
      blob->NInP()!=1 ||
      blob->InParticle(0)->Status()!=part_status::active)
  {
    msg.Error()<<METHOD<<" Blob or particle have wrong status."<<endl;
    return Return_Value::Error;
  }
  
  // choose decay channel or use the one chosen before for this blob
  if( inpart->Flav().IsStable() || !FindDecay(inpart->Flav()) ) return Return_Value::Nothing;
  Hadron_Decay_Channel* hdc = ChooseDecayChannel(blob);
  if( (inpart->Flav().IsAnti()) && !(hdc->FullDecay()&1) ) return Return_Value::Nothing;
  if(inpart->FinalMass()<hdc->DecayChannel()->MinimalMass()) return Return_Value::Retry_Method;
  
  FlavourSet daughters = hdc->DecayChannel()->GetDecayProducts();
  const int  n         = hdc->DecayChannel()->NumberOfDecayProducts()+1;

  // add daughters to blob
  Vec4D momentum; Flavour flav; Particle* particle=NULL;
  for( FlavourSet::iterator dpit = daughters.begin(); dpit != daughters.end(); dpit++ ) {
    flav = (*dpit);
    if( inpart->Flav().IsAnti() ) flav = flav.Bar();
    particle = new Particle( 0, flav );
    particle->SetStatus(part_status::active);
    particle->SetInfo('D');
    blob->AddToOutParticles( particle );
  }
  Particle_Vector outparticles = blob->GetOutParticles();

  // choose a kinematics that corresponds to the ME kinematic distribution
  vector<Vec4D> moms(n);
  moms[0] = inpart->Momentum();
  if( n<3 ) moms[1] = moms[0];
  else      {
    if(m_sc) { // weight correction if spin correlations are enabled
      // fixme : if all particles have to be retried for mass smearing
      Particle_Vector particles;
      particles.push_back(inpart);
      Particle_Vector::iterator pit;
      for(pit=outparticles.begin();pit!=outparticles.end();pit++) {
        particles.push_back(*pit);
      }

      Poincare labboost(labmom);
      labboost.Invert();
      vector<Vec4D> boosted_moms(n);
      boosted_moms[0]=labboost*moms[0];

      Blob* motherblob = inpart->ProductionBlob();
      Blob_Data_Base* scdata = (*motherblob)["amps"];
      if(scdata) {
        Amplitude_Tensor* amps = new Amplitude_Tensor(particles);
        Amplitude_Tensor* motheramps = scdata->Get<Amplitude_Tensor*>();
        Amplitude_Tensor* contractedamps = NULL;
        double motherampssumsquare = motheramps->SumSquare();

        double weight;
        int trials=0;
        do {
          ChooseDecayKinematics( moms, hdc, inpart->Flav().IsAnti() );
          // Calculate amplitude in lab frame for correlation
          for( int i=1; i<n; i++ ) boosted_moms[i]=labboost*moms[i];
          hdc->CalculateAmplitudes( &boosted_moms[0], amps, inpart->Flav().IsAnti() );

          double denominator = amps->SumSquare()*motherampssumsquare;

          if(contractedamps) { delete contractedamps; contractedamps = NULL; }
          contractedamps = new Amplitude_Tensor(Contraction(inpart,motheramps,amps));
          double numerator = contractedamps->SumSquare();
          trials++;
          weight = numerator/denominator;
          if(weight>1.0) {
            PRINT_INFO("Error: weight="<<weight<<" in "<<rpa.gen.NumberOfDicedEvents());
            PRINT_INFO(*blob);
          }
          if(trials>1000) PRINT_INFO("trials="<<trials);
        } while( ran.Get() > weight );

        blob->AddData("amps",new Blob_Data<Amplitude_Tensor*>(contractedamps));
        motheramps->Recreate(contractedamps); // save the contraction with the mother
        delete amps; amps=NULL;
      }
      else { // if first particle in SC chain
        ChooseDecayKinematics( moms, hdc, inpart->Flav().IsAnti() );
        for( int i=1; i<n; i++ ) boosted_moms[i]=labboost*moms[i];
        Amplitude_Tensor* amps = new Amplitude_Tensor(particles);
        hdc->CalculateAmplitudes( &boosted_moms[0], amps, inpart->Flav().IsAnti() );
        blob->AddData("amps",new Blob_Data<Amplitude_Tensor*>(amps));
      }
    }
    else {
      ChooseDecayKinematics( moms, hdc, inpart->Flav().IsAnti());
    }
  }

  for( size_t i=0; i<outparticles.size(); i++ ) {
    outparticles[i]->SetMomentum(moms[i+1]);
  }

  if( inpart->Info() == 'P' ) inpart->SetInfo('p');
  if( inpart->Info() == 'D' ) inpart->SetInfo('d');
  return Return_Value::Success;
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
    msg.Error()<<"Warning! The file "<<m_path<<m_constfile<<" does not exist"<<endl
             <<"     or has some syntax error."<<endl;
    msg.Error()<<"     Will ignore it and hope for the best."<<endl;
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
    msg.Tracking()<<om::blue<<"Found "<<dt->NumberOfDecayChannels()<<" decay channels for "<<dt->Flav()<<om::reset<<endl;
    dtreader->FillInMatrixElementsAndPS(dt,p_channelmap,m_md0);
    if( msg.LevelIsInfo() ) {
      msg.Tracking()<<"Initialised a new decay table : "<<endl;
      msg.Tracking()<<"(Using information given in .dat files, such as BR, width, ...)"<<endl;
      if(msg.LevelIsTracking()) dt->Output();
      msg.Tracking()<<"Calculated decay widths using implemented models :"<<endl
        <<"(only for information; they are NOT used for the simulation)"<<endl;
      for (int i=0;i<dt->NumberOfDecayChannels();i++) {			
        msg.Tracking()
          <<(*p_channelmap)[dt->GetDecayChannel(i)]->ChannelName()<<" : "
          <<(*p_channelmap)[dt->GetDecayChannel(i)]->GetPS()->Result()<<" ("
          <<(*p_channelmap)[dt->GetDecayChannel(i)]->GetPS()->RelError()<<" %)"
          <<" GeV"<<endl;
      }
      msg.Tracking()<<endl;
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
