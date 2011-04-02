#include "HADRONS++/Main/Hadrons.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include <time.h>
#include "HADRONS++/Main/Hadron_Decay_Map.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include <fstream>
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <utility>
#include <algorithm>
#include "ATOOLS/Phys/Decay_Table.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Matrix.H"
#include "ATOOLS/Math/Poincare.H"
#include "HADRONS++/PS_Library/HD_PS_Base.H"
#include "HADRONS++/Main/Mixing_Handler.H"
#include "METOOLS/Main/Spin_Structure.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Hadrons::Hadrons( string path, string file, string constfile ) :
  m_spincorrelations(false)
{
  msg_Tracking()<<"In Hadrons: ("<<path<<") "<<file<<std::endl;
  msg_Tracking()<<"In Hadrons: ("<<path<<") "<<constfile<<std::endl;
  p_decaymap = new Hadron_Decay_Map(path, file, constfile);
  p_decaymap->ReadInConstants();
  p_decaymap->ReadHadronAliases(path, "HadronAliases.dat");
  p_decaymap->SetHadronProperties();
  p_decaymap->Read(path, file);
  p_decaymap->Read(path, "AliasDecays.dat");
  p_decaymap->Initialise();
  p_decaymap->ReadFixedTables();
  p_mixinghandler = new Mixing_Handler(p_decaymap);
  p_color_unitmatrix = new CMatrix(1);
  (*p_color_unitmatrix)[0][0] = Complex(1.0,0.0);
}

Hadrons::~Hadrons()
{
  delete p_decaymap; p_decaymap=NULL;
  delete p_mixinghandler; p_mixinghandler=NULL;
  delete p_color_unitmatrix; p_color_unitmatrix=NULL;
}

bool Hadrons::CreateDecayBlob(Blob* blob)
{
  Hadron_Decay_Table* table = p_decaymap->FindDecay(blob->InParticle(0)->Flav());
  if(table==NULL) return false;
  ChooseDecayChannel(blob,table);
  return true;
}

bool Hadrons::FillDecayBlob(Blob* blob, const Vec4D& labmom)
{
  DEBUG_FUNC("blob->Id()="<<blob->Id());
  Hadron_Decay_Table* table = p_decaymap->FindDecay(blob->InParticle(0)->Flav());
  if(table==NULL) return false;
  Hadron_Decay_Channel* hdc = ChooseDecayChannel(blob,table);
  if(!hdc) return false;
  
  Particle*       inpart = blob->InParticle(0);
  FlavourSet daughters = hdc->GetDecayProducts();
  Flavour flav; Particle* particle=NULL;
  for( FlavourSet::iterator dpit = daughters.begin(); dpit != daughters.end(); dpit++ ) {
    flav = (*dpit);
    if( inpart->Flav().IsAnti() ) flav = flav.Bar();
    particle = new Particle( 0, flav );
    particle->SetFinalMass();
    particle->SetStatus(part_status::active);
    particle->SetNumber();
    particle->SetInfo('D');
    blob->AddToOutParticles( particle );
  }
  
  ChooseDecayKinematics( blob, labmom, hdc);
  if(!SetColorFlow( blob )) {
    msg_Error()<<METHOD<<" wasn't able to set the color flow for"<<endl<<*blob<<endl;
  }
  inpart->SetStatus(part_status::decayed);
  DEBUG_VAR(*blob);
  return true;
}

void Hadrons::ChooseDecayKinematics(
                                    Blob* blob,
                                    const Vec4D& labmom,
                                    Hadron_Decay_Channel* hdc
                                   )
{
  Particle*       inpart = blob->InParticle(0);
  bool            anti   = inpart->Flav().IsAnti();
  Particle_Vector outparticles = blob->GetOutParticles();
  size_t          n = outparticles.size()+1;
  
  vector<Vec4D> moms(n);
  moms[0] = inpart->Momentum();
  
  if(n<3) moms[1] = moms[0];
  else {
    if(m_spincorrelations) { // weight correction if spin correlations are enabled
      Particle_Vector particles;
      particles.push_back(inpart);
      particles.insert(particles.end(),outparticles.begin(),outparticles.end());
      
      Poincare labboost(labmom);
      labboost.Invert();
      p_spblob->SetCMS();
      Poincare spcmsboost = Poincare( p_spblob->CMS() );
      vector<Vec4D> boosted_moms(n);
      boosted_moms[0]=spcmsboost*(labboost*moms[0]);

      Particle* contracting_part=inpart;
      Amplitude_Tensor* motheramps = GetMotherAmplitudes(inpart,contracting_part,labmom);
      if(motheramps) {
        Amplitude_Tensor* amps = new Amplitude_Tensor(particles);
        amps->SetColorMatrix(p_color_unitmatrix);
        double motherampssumsquare = motheramps->SumSquare();
        
        double weight;
        int trials=0;
        do {
          GenerateUncorrelatedKinematics( moms, hdc, anti );
          
          // Calculate amplitude in signal process cms frame for correlation
          for( size_t i=1; i<n; i++ ) boosted_moms[i]=spcmsboost*(labboost*moms[i]);
          hdc->CalculateAmplitudes( &boosted_moms[0], amps, anti );
          
          double denominator = amps->SumSquare()*motherampssumsquare;
          double numerator = motheramps->SoftContract(contracting_part,
              amps,
              inpart);
          weight = numerator/denominator;
          if(weight>(1.0+Accu())) {
            msg_Error()<<METHOD<<" Error: weight="<<weight<<" in "
                <<rpa.gen.NumberOfGeneratedEvents()<<endl
                <<(*blob)<<endl;
          }
          trials++;
          if(trials%1000==0) {
            msg_Error()<<METHOD<<" Warning: spin correlation trials="<<trials<<endl
                <<"    possibly wrong amplitudes?"<<endl
                <<(*blob)<<endl;
          }
        } while( ran.Get() > weight );

        motheramps->Contract(contracting_part,amps,inpart);
        delete amps; amps=NULL;
      }
      else { // if first particle in SC chain
        GenerateUncorrelatedKinematics( moms, hdc, anti );
        
        for( size_t i=1; i<n; i++ ) boosted_moms[i]=spcmsboost*(labboost*moms[i]);
        Amplitude_Tensor* amps = new Amplitude_Tensor(particles);
        amps->SetColorMatrix(p_color_unitmatrix);
        hdc->CalculateAmplitudes( &boosted_moms[0], amps, anti );
        blob->AddData("amps",new Blob_Data<Amplitude_Tensor*>(amps));
      }
    }
    else {
      GenerateUncorrelatedKinematics( moms, hdc, anti );
    }
  }
  
  for( size_t i=1; i<n; i++ ) {
    outparticles[i-1]->SetMomentum(moms[i]);
  }
}

Hadron_Decay_Channel * Hadrons::ChooseDecayChannel(Blob* blob, Hadron_Decay_Table* table)
{
  Blob_Data_Base* data = (*blob)["hdc"];
  if(data) {
    if(blob->Has(blob_status::internal_flag)) {
      bool partonic_finalstate(false);
      Hadron_Decay_Channel* hdc;
      do {
        hdc = table->Select();
        FlavourSet flavs=hdc->GetDecayProducts();
        for (FlSetIter fl=flavs.begin(); fl!=flavs.end(); fl++) {
          if(fl->Strong()) {
            partonic_finalstate=true;
            break;
          }
        }
      } while (!partonic_finalstate);
      DEBUG_INFO("retrying with "<<hdc->Name());
      blob->UnsetStatus(blob_status::internal_flag);
      blob->AddData("hdc",new Blob_Data<Hadron_Decay_Channel*>(hdc));
      return hdc;
    }
    else return data->Get<Hadron_Decay_Channel*>();
  }

  // set CP asymmetries due to interference between with mixing and without
  bool setasymmetries = p_mixinghandler->SetCPAsymmetries(blob->InParticle(0), table);

  // generate decay channel acc. to BR
  // fixme : treat K0 special
  Hadron_Decay_Channel* dec_channel = table->Select();

  if(setasymmetries) p_mixinghandler->ResetCPAsymmetries(table);

  blob->AddData("hdc",new Blob_Data<Hadron_Decay_Channel*>(dec_channel));
  return dec_channel;
}

void Hadrons::GenerateUncorrelatedKinematics(
                                          std::vector<ATOOLS::Vec4D> & moms,
                                          Hadron_Decay_Channel       * hdc,
                                          bool                         anti)
{
  static std::string mname(METHOD);
  rvalue.IncCall(mname);
  if(moms.size()==2) {
    moms[1]=moms[0];
    return;
  }
  double value(0.);
  const double max = hdc->GetPS()->Maximum();    // note: no flux in max.
  int trials(0);
  do {
    if(trials>10000) {
      msg_Error()<<METHOD<<"("<<hdc->Name()<<"): "
                 <<"Rejected hadron decay kinematics 10000 times. "
                 <<"This indicates a wrong maximum. "
                 <<"Will accept kinematics."
                 <<endl;
      rvalue.IncRetryMethod(mname);
      break;
    }
    value = hdc->Differential(&moms.front(),anti);
    if(value/max>1.05 && max>1e-30) {
      msg_Tracking()<<METHOD<<"("<<hdc->Name()<<") warning:"<<endl
                    <<"  d\\Gamma(x)="<<value<<" > max(d\\Gamma)="<<max
                    <<std::endl;
      hdc->GetPS()->SetMaximum(value);
      rvalue.IncRetryMethod(mname);
      break;
    }
    trials++;
  } while( ran.Get() > value/max );
}


void Hadrons::CleanUp()
{
  p_decaymap->ResetCounters();
}

void Hadrons::SetSpinCorrelations(bool sc)
{
  m_spincorrelations = sc;
}

bool Hadrons::SetColorFlow(Blob* blob)
{
  int n_q(0), n_g(0);
  for(int i=0;i<blob->NOutP();i++) {
    if(blob->OutParticle(i)->Flav().IsQuark())      n_q++;
    else if(blob->OutParticle(i)->Flav().IsGluon()) n_g++;
  }
  if(n_q>0 || n_g>0) {
    blob->SetStatus(blob_status::needs_showers);
    Hadron_Decay_Channel* hdc = (*blob)["hdc"]->Get<Hadron_Decay_Channel*>();
    return hdc->SetColorFlow(blob->GetOutParticles(),n_q,n_g);
  }
  else return true;
}

Amplitude_Tensor* Hadrons::GetMotherAmplitudes(Particle* part,
                                               Particle*& contracting_part,
                                               const Vec4D& labmom)
{
  /** try to go up step by step along _decay_ blobs, to find a blob where an
  amplitude tensor containing this particle has been stored */
  Blob* motherblob = part->ProductionBlob();
  while(motherblob) {
    Blob_Data_Base* scdata = (*motherblob)["amps"];
    if(scdata) {
      Amplitude_Tensor* amps = scdata->Get<Amplitude_Tensor*>();
      if(amps->Contains(part)) return amps;
      else                     return NULL;
    }
    else {
      if(motherblob->Type()==btp::Hadron_Decay && motherblob->InParticle(0)->ProductionBlob())
        motherblob = motherblob->InParticle(0)->ProductionBlob();
      else break;
    }
    /** stop the search if we encounter a pythia blob  (no spin corrrelations there) */
    if(motherblob->TypeSpec()=="Pythia_v6.214") return NULL;
  }

  /** if above didn't succeed, let's look if the signal process blob has an
  amplitude tensor for this particle */
  Blob_Data_Base* scdata = (*p_spblob)["amps"];
  if(scdata) {
    Amplitude_Tensor* amps = scdata->Get<Amplitude_Tensor*>();
    /** If it directly contains our particle -> fine, return. */
    if(amps->Contains(part)) return amps;
    /** If it contains a duplicate of our particle, e.g. because there was a
    duplicating fragmentation blob or ME_PS_Interface inbetween: */
    else if(amps->Contains(part->OriginalPart())) {
      contracting_part = part->OriginalPart();
      return amps;
    }
    else return NULL; // didn't find a duplicate of our particle
  }
  else return NULL;
}
