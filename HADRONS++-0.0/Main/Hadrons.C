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
#include "Channel_Elements.H"
#include "Momenta_Stretcher.H"
#include "Poincare.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Hadrons::Hadrons( string _path, string _file, string _constfile ) : 
  m_path(_path), m_file(_file), m_constfile(_constfile), m_createbooklet(0)
{ 
  msg_Tracking()<<"In Hadrons: ("<<_path<<") "<<_file<<std::endl;
  msg_Tracking()<<"In Hadrons: ("<<_path<<") "<<_constfile<<std::endl;
  ReadInConstants();
  bool createbooklet (false);
  ReadInDecayTables(createbooklet);
  if( createbooklet ) CreateBookletNow();
}


void Hadrons::CreateBookletNow()
{
  string fn = string("hadrons");
  string texfn = m_path + fn + string(".tex");
  ofstream f(texfn.c_str());
   
  // header
  f<<"\\documentclass[a4paper]{article}\n"
   <<"\\usepackage{latexsym,amssymb,amsmath,amsxtra,longtable,fullpage}\n\n"
   <<"\\begin{document}\n"<<endl; 
  f<<"\\newcommand{\\m}{-}"<<endl; 
  f<<"\\newcommand{\\p}{+}"<<endl; 
  f<<"\\title{Decay Channels of the {\\tt HADRONS++} Module}\n\\maketitle"<<endl;
  f<<"\\tableofcontents"<<endl<<endl;

  // text 
  Decay_Channel * dc (NULL);
  FlavourSet outs;
  char helpstr[10];
  for ( map<kf::code,Decay_Table *>::iterator pos = p_decaymap->begin(); pos != p_decaymap->end(); ++pos) {
	kf::code      kfc (pos->first);
	Decay_Table * dt  (pos->second);
	f<<"\\section{Decaying Particle: $"<<dt->Flav().TexName()<<"$}"<<endl;
	f<<"\\begin{tabular}{ll}"<<endl;
	f<<" number of decay channels:	& "<<dt->NumberOfDecayChannels()<<"\\\\ "<<endl;
	f<<" total width:               & "<<dt->TotalWidth()<<" GeV \\\\ "<<endl;
    f<<" experimental width:        & "<<dt->Flav().Width()<<" GeV \\\\ "<<endl;
	f<<"\\end{tabular}"<<endl;
    f<<"\\begin{longtable}[l]{lll}"<<endl;
    f<<"\\multicolumn{3}{c}{\\bf Exclusive Decays}\\\\"<<endl;
    f<<"\\hline"<<endl;
    f<<"Decay Channel & Branching Ratio & Origin \\\\"<<endl;
    f<<"\\hline\n\\hline"<<endl;
    for( size_t i=0; i<dt->NumberOfDecayChannels(); ++i ) {
      dc = dt->GetDecayChannel(i);
      outs = dc->GetDecayProducts();
      f<<"$"<<dc->GetDecaying().TexName()<<"$ $\\to$ ";
      for (FlSetConstIter fl=outs.begin();fl!=outs.end();++fl) f<<"$"<<fl->TexName()<<"$ ";
      f<<" & ";
      sprintf( helpstr, "%.4f", dc->Width()/dt->TotalWidth()*100. );
      f<<helpstr;
      if( dc->DeltaWidth() > 0. ) {
        sprintf( helpstr, "%.4f", dc->DeltaWidth()/dt->TotalWidth()*100. );
        f<<" $\\pm$ "<<helpstr;
      }
      f<<" \\% & \\verb;"<<dc->Origin()<<";\\\\"<<endl;
    }
    f<<"\\hline"<<endl;
    map<string,pair<double,double> > brmap;
    GetInclusives(dt,brmap);
    f<<"\\multicolumn{3}{c}{\\hfill}\\\\"<<endl;
    f<<"\\multicolumn{3}{c}{\\bf Inclusive Decays with $BR>10^{-6}$}\\\\"<<endl;
    f<<"\\hline"<<endl;
    f<<"Decay Channel & Branching Ratio \\\\"<<endl;
    f<<"\\hline\n\\hline"<<endl;
    for( BRMapString::iterator brit=brmap.begin(); brit != brmap.end(); ++brit ) {
      double br  = (brit->second.first+brit->second.second)/2.;
      double dbr = (brit->second.first-brit->second.second)/2.;
      if( br>1.e-6 ) {
        f<<"$"<<dc->GetDecaying().TexName()<<"$ $\\to$ $"<<brit->first<<"$ & ";
        sprintf( helpstr, "%.4f", br*100. );
        f<<helpstr;
        if( dbr > 0. ) {
          sprintf( helpstr, "%.4f", dbr*100. );
          f<<" $\\pm$ "<<helpstr;
        }
        f<<" \\% \\\\"<<endl;
      }
    }
    f<<"\\hline"<<endl;
    f<<"\\end{longtable}"<<endl;
  }
   
  // end 
  f<<"\\end{document}"<<endl; 
  msg.Out()<<"Created HADRONS++ booklet:"<<endl;
  msg.Out()<<"Run "<<om::bold<<"latex "<<fn<<om::reset<<" in "<<m_path<<" for compilation."<<endl;
  abort();
}

std::vector<BRPairFlavourSet> Hadrons::GetInclusives( 
    Decay_Table * dt,                               // decay table to be considered
    BRMapString & brmap,                            // map with flavourset <-> (upper, lower)
    FlavourSet flset,                               // flavour set
    DoublePair br,                                  // upper, lower br
    bool eoi )                                      // end of iteration
{
  Decay_Channel * dc;
  FlavourSet outs;
  vector<BRPairFlavourSet> new_flset, flvec;
  // for each decay channel
  for( size_t I=0; I<dt->NumberOfDecayChannels(); ++I ) {
    dc = dt->GetDecayChannel(I);                            // decay channel
    outs = dc->GetDecayProducts();                          // FlavourSet of products
    // for each decay product
    new_flset.clear();
    double dbr = (dc->DeltaWidth()>0.)? dc->DeltaWidth()/dt->TotalWidth() : 0.;
    new_flset.push_back(BRPairFlavourSet(
          flset,
          DoublePair(
            br.first*(dc->Width()/dt->TotalWidth()+dbr),
            br.second*(dc->Width()/dt->TotalWidth()-dbr))
          ));
    for (FlSetConstIter fl=outs.begin();fl!=outs.end();++fl) {
      if (p_decaymap->find(fl->Kfcode())==p_decaymap->end() ||
          fl->IsStable()) {                                  // if daughter is stable
        for( size_t i=0; i<new_flset.size(); ++i ) {
          new_flset[i].first.insert(*fl);
        }
      }
      else {                                                 // if daughter has DT 
        vector< vector<BRPairFlavourSet> > dauchans; 
        for( size_t i=0; i<new_flset.size(); ++i ) {
          dauchans.push_back( GetInclusives( (*p_decaymap)[fl->Kfcode()], brmap, new_flset[i].first, new_flset[i].second, 0 ) );
        }
        new_flset.clear();
        for( size_t i=0; i<dauchans.size(); ++i ) 
          for( size_t j=0; j<dauchans[i].size(); ++j )
            new_flset.push_back( dauchans[i][j] );
      }
    }
    // end of iteration?
    if( eoi ) {
      // add to map
      for( size_t i=0; i<new_flset.size(); ++i ) {
        string channel = string("");
        for (FlSetConstIter fl=new_flset[i].first.begin();fl!=new_flset[i].first.end();++fl) 
          channel += fl->TexName()+string(" ");
        brmap[channel].first += new_flset[i].second.first;          // upper limit
        brmap[channel].second += new_flset[i].second.second;        // lower limit
      }
    }
    else {
      for( size_t i=0; i<new_flset.size(); ++i ) 
        flvec.push_back(new_flset[i]);
    }
  }
  if( flvec.size() ) {
  }
  return flvec;
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
  int trials(0);                                            // number of trials
  do {
    value = _hdc->Differential(_p,sigma);                   // current val. of T
    if( value/max>1.+1.e-4 ) {
      msg.Error()<<om::red<<"Warning in Hadrons::ChooseDecayKinematics for "<<_hdc->ChannelName()<<" with sigma @ "<<sigma<<endl;
      msg.Error()<<"  "<<value<<" > "<<max<<"   factor 1.+"<<value/max-1.<<"     @ trial "<<trials+1<<om::reset<<endl;
      msg.Debugging()<<om::green<<endl;
    }
    trials++;
  } while( ran.Get() > value/max );
  msg.Debugging()<<"          Decay kinematics before mass smearing:"<<endl;
  for(int i=0;i<_hdc->NOut()+1;i++) {
    msg.Debugging()<<"            p["<<i<<"]="<<_p[i]<<": ["<<_p[i].Mass()<<"]"<<endl;
  }
}


ATOOLS::Blob* Hadrons::CreateDecayBlobSkeleton(
    Particle            * particle, 
    Blob_List           * blob_list, 
    Particle_List       * part_list )
{
  if( particle->Flav().IsStable() || !FindDecay(particle->RefFlav()) ) {
    msg.Debugging()<<"            is stable."<<endl;
    return NULL; // don't create decayskeleton
  }
  
  // choose decay channel for daughter acc. to BR
  Hadron_Decay_Channel * hdc = ChooseDecayChannel();
  
  if( particle->Flav().IsAnti() && !(hdc->FullDecay() & 1) ) { // if NO_ANTI
    msg.Debugging()<<"            anti particle does not decay."<<endl;
    return NULL;
  }
  
  // Create daughters decay blob and add incoming particle
  Blob * blob;
  blob = new Blob();
  blob->SetStatus(3);
  blob->SetType( btp::Hadron_Decay );
  blob->SetTypeSpec("Sherpa");
  blob->SetId();
  double        time = particle->LifeTime();
  Vec3D      spatial = particle->Distance( time );
  Vec4D     position = Vec4D( time*rpa.c(), spatial );
  blob->SetPosition( particle->XProd() + position );
  blob->AddToInParticles( particle );
  if( particle->Info() == 'P' ) particle->SetInfo('p');
  if( particle->Info() == 'D' ) particle->SetInfo('d');
  msg.Debugging()<<"            created decayblob skeleton ["<<blob->Id()<<"]."<<endl;
    // fill it with outgoing particles (no momenta yet)
  FlavourSet daughters = hdc->DecayChannel()->GetDecayProducts();
  FlavourSet::iterator dit;
  for(dit = daughters.begin(); dit != daughters.end(); dit++ ) {
    Flavour daughter_flav=(*dit);
    if( particle->Flav().IsAnti() ) daughter_flav = daughter_flav.Bar();
    Particle* daughter_particle = new Particle( -1, daughter_flav );
    msg.Debugging()<<"            found daughter "<<daughter_flav<<": adding it to blob ["<<blob->Id()<<"]"<<endl;
    if( part_list ) {
      daughter_particle->SetNumber( part_list->size() );
      part_list->push_back( particle );
    }
    else daughter_particle->SetNumber( 0 );
    daughter_particle->SetStatus(part_status::active);
    daughter_particle->SetInfo('D');
    blob->AddToOutParticles( daughter_particle );
  }

  blob_list->push_back( blob );
  m_blob_hdc_map.insert( pair<int,Hadron_Decay_Channel*>(blob->Id(),hdc) );
  return blob;
}



Spin_Density_Matrix Hadrons::PerformDecay( 
    Blob                * blob, 
    Blob_List           * blob_list, 
    Particle_List       * part_list,
    Spin_Density_Matrix * sigma )
{
  msg.Debugging()<<"      Performing decay for blob ["<<blob->Id()
		 <<"] in 2 steps. 1st creating skeleton for daughter decay blobs, "
		 <<"then dicing kinematics for the blob itself."<<endl;
  Particle* part = blob->InParticle(0);
  msg_Tracking()<<"Hadrons::PerformDecay() for "<<part->Flav()<<" with spin correlation @ "<<sigma<<endl;
  msg.Debugging()<<"        Momentum of inparticle: "<<part->Momentum()<<endl;

  // find decay channel for current decay blob (it was mapped during the skeleton creation)
  Hadron_Decay_Channel * hdc;
  map<int,Hadron_Decay_Channel*>::iterator iter = m_blob_hdc_map.find(blob->Id());
  if (iter != m_blob_hdc_map.end()) {
    hdc = iter->second;
  }
  else {
    msg.Error()<<"Hadron decay channel for blob ["<<blob->Id()<<"] not found"
        <<"altough it should exist already. Aborting ..."<<std::endl;
    abort();
  }
  
  bool       anti = part->Flav().IsAnti();			    // antiparticle ?
  int mother_spin = part->Flav().IntSpin();             // 2*spin
  
  // create decay blob skeleton for daughters' decay (already needed because of mass smearing)
  Particle_Vector daughters = blob->GetOutParticles();
  Particle_Vector::iterator dit;
  if(hdc->FullDecay() & 2) {  // only if NO_FULLDECAY is not specified
    msg.Debugging()<<"        Creating skeleton for daughters of blob ["<<blob->Id()<<"]"<<endl;
    for( dit = daughters.begin(); dit != daughters.end(); dit++ ) {
      msg.Debugging()<<"          Daughter of blob ["<<blob->Id()<<"]: "<<(*dit)->Flav()<<": "<<(*dit)<<endl;
      Blob* daughter_blob=CreateDecayBlobSkeleton((*dit),blob_list,part_list);
    }
  }
  
  msg.Debugging()<<"        Dicing kinematics for ["<<blob->Id()<<"] with hadron decay channel "<<hdc->ChannelName()<<endl;
  

  if( (anti) && !(hdc->FullDecay() & 1) ) {             // if NO_ANTI
    if( sigma ) {
      Spin_Density_Matrix unitmatrix = Spin_Density_Matrix(mother_spin);
      unitmatrix.SetUnitMatrix();
      return unitmatrix;                                // return unit matrix
    }
    return Spin_Density_Matrix();                       // return empty decay matrix
  }
//  cout<<om::red<<hdc->ChannelName()<<om::reset<<endl;
//  if( sigma ) sigma->Print();

  // actually decay the current blob (create its kinematics)   
  // choose a kinematics that corresponds to the ME kinematic distribution
  if( !FindDecay(part->RefFlav()) ) {
    msg.Error()<<"Hadrons::PerformDecay cannot decay particle "<<part->RefFlav()<<std::endl;
    abort();
  }
  
  const int n = blob->NOutP()+1;
  Vec4D mom[n];
  mom[0] = part->Momentum();
  
  if(n>2) {
    ChooseDecayKinematics( mom, hdc, sigma );
    if(hdc->MassSmearing()) {  // if mass smearing is turned on
      // sort daughter particle indices by width
      int by_width[n];
      if(n>2) {
        for(int i=0;i<n;i++) {
          by_width[i]=i;
        }
        // sort "by_width"-array by the width of daughters (insertion sort)
        for(int i=2; i<n; ++i) {
          int j;
          int value = by_width[i];
          for (j=i-1; j>=1 && daughters[by_width[j]-1]->Flav().Width() > 
                              daughters[i-1]->Flav().Width(); --j) {
            by_width[j+1] = by_width[j];
          }
          by_width[j+1] = value;
        }
      }
      double masses[n];
      double min_masses[n];
      double max_masses[n];
      double widths[n];
      masses[0]=part->Momentum().Mass();
      msg.Debugging()<<"          masses before smearing:"<<endl;
      msg.Debugging()<<"            masses[0]="<<masses[0]<<" for "<<part->Flav()<<endl;
      // get minimum mass as min_daughter = sum ( min_daughter_daughters )
      for(int ii=1;ii<n;ii++) {
        int i=by_width[ii];
//       for(int i=1;i<n;i++) {
        widths[i]=daughters[i-1]->Flav().Width();
        masses[i]=daughters[i-1]->Flav().Mass(); // temporary storage of peak mass
        min_masses[i]=0.0;
        if(daughters[i-1]->DecayBlob()) {
          Particle_Vector daughter_daughters = daughters[i-1]->DecayBlob()->GetOutParticles();
          for(int j=0;j<daughter_daughters.size();j++) {
            min_masses[i]+=daughter_daughters[j]->Flav().Mass();
          }
        }
      }
      // get maximum masses and smear masses
      for(int ii=1;ii<n;ii++) {
        int i=by_width[ii];
        max_masses[i]=masses[0];
        for(int kk=1;kk<i;kk++) {
          int k=by_width[kk];
          max_masses[i]-=masses[k];
        }
        for(int jj=i+1;jj<n;jj++) {
          int j=by_width[jj];
          max_masses[i]-=min_masses[j];
        }
        msg.Debugging()<<"            masses["<<i<<"]="<<masses[i]
                       <<" for "<<daughters[i-1]->Flav()<<endl;
        // BreitWigner with mass, width, min_mass, max_mass
        if( masses[i] > 1e-5 && widths[i] > 1e-5 ) {
          double myrandom = ran.Get();
          msg.Debugging()<<"            "<<daughters[i-1]->Flav()
                         <<": smearing with mass="<<masses[i]
                         <<" width="<<widths[i]
                         <<" min="<<min_masses[i]
                         <<" max="<<max_masses[i]
                         <<" random="<<myrandom<<endl;
          masses[i] = sqrt( PHASIC::CE.MassivePropMomenta(masses[i],widths[i],
                            1,sqr(min_masses[i]),sqr(max_masses[i]),myrandom) );
        }
      }
      msg.Debugging()<<"          diced target masses:"<<endl;
      for(int i=1;i<n;i++) {
        msg.Debugging()<<"            masses["<<i<<"]="<<masses[i]
                       <<" for "<<daughters[i-1]->Flav()<<endl;
      }
      Vec4D cms = Vec4D(0.,0.,0.,0.);
      for(int i=1;i<hdc->NOut()+1;i++) {
        cms += mom[i];
      }
      Poincare boost(cms);
      for(int i=1;i<hdc->NOut()+1;i++) {
        boost.Boost(mom[i]);
      }
      Momenta_Stretcher stretch;
      bool okay = stretch.ZeroThem(1,hdc->NOut()+1, mom);
      okay = okay && stretch.MassThem(1,hdc->NOut()+1, mom, masses);
      if(!okay) {
        msg.Error()<<METHOD<<": Momenta_Stretcher delivered an error for stretching "
                   <<hdc->ChannelName()<<std::endl;
        msg.Error()<<METHOD<<":  "<<"min_masses["<<0<<"]="<<min_masses[0]
                   <<" masses["<<0<<"]="<<masses[0]<<" max_masses["<<0<<"]="
                   <<max_masses[0]<<std::endl;
        for(int ii=1;ii<n;ii++) {
          int i=by_width[ii];
          msg.Error()<<METHOD<<":  "<<"min_masses["<<i<<"]="<<min_masses[i]
                     <<" masses["<<i<<"]="<<masses[i]<<" max_masses["<<i<<"]="
                     <<max_masses[i]<<std::endl;
        }
      }
      for(int i=1;i<hdc->NOut()+1;i++) {
        boost.BoostBack(mom[i]);
      }
      msg.Debugging()<<"          Decay kinematics after mass smearing:"<<endl;
      for(int i=0;i<hdc->NOut()+1;i++) {
        msg.Debugging()<<"            p["<<i<<"]="<<mom[i]
                       <<": ["<<mom[i].Mass()<<"]"<<endl;
      }
    }
  }
  else { // if n=2
    mom[1] = mom[0];
  }  
  /*------------------------------*/
  

  blob->SetStatus(1);
  int i(0);
  vector<pair<Particle*,int> > daughters_pair;
  for( dit = daughters.begin(); dit != daughters.end(); dit++ ) {
    i++;
    (*dit)->SetMomentum(mom[i]);
    daughters_pair.push_back( pair<Particle*,int>( (*dit), i ) );
  }

  msg.Debugging()<<"        Finished with kinematics for ["<<blob->Id()<<"]. Going on recursively if necessary..."<<endl;
  
  // create new spin correlation tensor and shuffle randomly daughters
  Spin_Correlation_Tensor * SCT (NULL);
  if( Spin_Correlation_Tensor::Mode()!=scmode::None && hdc->GetIndexList()->size() ) {
    SCT = 
      new Spin_Correlation_Tensor( hdc->GetIndexList(), hdc->GetAmplitudeTensor() );
    random_shuffle( daughters_pair.begin(), daughters_pair.end() );
  }

  // treat every daughter
  vector<pair<Particle*,int> >::iterator dpit;
  for( dpit = daughters_pair.begin(); dpit != daughters_pair.end(); dpit++ ) {
	// check if daughter can be treated as well
    Particle* daughter_particle = dpit->first;
    Blob* daughter_decay_blob = daughter_particle->DecayBlob();
    int spin = daughter_particle->Flav().IntSpin();    // 2*spin
	if ( ( hdc->FullDecay() & 2 ) &&
         ( !(daughter_particle->Flav().IsStable()) ) &&
         ( FindDecay(daughter_particle->RefFlav()) )  
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
        (*decay_matrix) = PerformDecay( daughter_decay_blob, blob_list, part_list, new_sigma );
        decay_matrix->Normalise();
        // contract over decay matrix to reduce the size of the SCT
        SCT->Contract( i, decay_matrix );
        delete new_sigma;
        delete decay_matrix;
      }
      else PerformDecay( daughter_decay_blob, blob_list, part_list, NULL );
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
    if( SCT->GetDepth()!= 0 ) {
      PRINT_INFO("not eth. was contracted!");
      cout<<"SCT : "<<(*SCT)<<" with depth "<<SCT->GetDepth()<<endl;
      PRINT_INFO(om::red<<"will abort."<<om::reset);
      abort();
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
 
void Hadrons::ReadInDecayTables( bool & createbooklet )
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
    if( Decayers[i][0] != string("CREATE_BOOKLET") ) {
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
    else createbooklet = true;
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
    msg.Info()<<om::blue<<"Found "<<dt->NumberOfDecayChannels()<<" decay channels for "<<dt->Flav()<<om::reset<<endl;
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
