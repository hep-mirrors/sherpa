#include "HADRONS++/Main/Mixing_Handler.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "HADRONS++/Main/Hadron_Decay_Table.H"
#include "HADRONS++/Main/Hadron_Decay_Channel.H"
#include "HADRONS++/Main/Hadron_Decay_Map.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

Mixing_Handler::Mixing_Handler()
{
}

Mixing_Handler::~Mixing_Handler()
{
}

double Mixing_Handler::DetermineMixingTime(Particle* decayer,
                                           bool checkforpartstatus) const
{
  double time = decayer->Time();

  Blob* motherblob = decayer->ProductionBlob();
  if(motherblob->Type()==btp::Hadron_Mixing) {
    decayer = motherblob->InParticle(0);
    motherblob = decayer->ProductionBlob();
  }
  Particle* sister = NULL;
  if(motherblob->Type()!=btp::Fragmentation ||
     (motherblob->NInP()==2 && motherblob->NOutP()==2 &&
      motherblob->InParticle(0)->Flav()==motherblob->OutParticle(0)->Flav() &&
      motherblob->InParticle(1)->Flav()==motherblob->OutParticle(1)->Flav())) {
    // check if particle was produced coherently
    Particle_Vector sisters = motherblob->GetOutParticles();
    for(Particle_Vector::const_iterator it=sisters.begin(); it!=sisters.end(); it++) {
      if((*it)!=decayer && decayer->Flav()==(*it)->Flav().Bar()) {
        sister = (*it);
        break;
      }
    }
  }
  // in coherent production, special rules apply:
  if(sister) {
    if((checkforpartstatus && sister->Status()==part_status::decayed) ||
       (!checkforpartstatus && sister->DecayBlob()))
    {
      time = time - sister->Time();
    }
    else time = 0.0;
  }
  return time;
}


ATOOLS::Blob* Mixing_Handler::PerformMixing(Particle* decayer) const
{
  // explicit mixing in event record
  Flavour flav = decayer->Flav();
  string tag = flav.IsAnti() ? flav.Bar().IDName() : flav.IDName();
  if(m_model("Mixing_"+tag,0.0)!=0.0 && decayer->Info()!=char('M')) {
    double t = DetermineMixingTime(decayer,true)/rpa->hBar();
    if(t==0.0) return false;
    double factor = decayer->Flav().QOverP2();
    if(decayer->Flav().IsAnti()) factor = 1.0/factor;
    double dG = decayer->Flav().DeltaGamma()*t/4.0;
    double dm = decayer->Flav().DeltaM()*t/2.0;
    Complex i(0.0,1.0);
    double prob_not_mix = sqr(abs(exp(i*dm)*exp(dG)+exp(-i*dm)*exp(-dG)));
    double prob_mix = factor*sqr(abs(exp(i*dm)*exp(dG)-exp(-i*dm)*exp(-dG)));
    if(prob_mix > ran->Get()*(prob_mix+prob_not_mix)) {
      if(decayer->DecayBlob()) {
        decayer->DecayBlob()->RemoveOwnedParticles();
        delete decayer->DecayBlob();
      }
      decayer->SetStatus(part_status::decayed);
      decayer->SetInfo('m');
      Particle* mixed_part = new Particle(0, decayer->Flav().Bar(),
                                          decayer->Momentum(), 'M');
      mixed_part->SetFinalMass(decayer->FinalMass());
      mixed_part->SetStatus(part_status::active);
      mixed_part->SetTime(decayer->Time());
      Blob* mixingblob = new Blob();
      mixingblob->SetType(btp::Hadron_Mixing);
      mixingblob->SetId();
      mixingblob->SetStatus(blob_status::inactive);
      mixingblob->SetTypeSpec("HADRONS");
      mixingblob->AddToInParticles(decayer);
      mixingblob->AddToOutParticles(mixed_part);
      mixingblob->SetPosition(decayer->ProductionBlob()->Position());
      mixingblob->SetStatus(blob_status::needs_hadrondecays);
      return mixingblob;
    }
  }
  return NULL;
}


Hadron_Decay_Channel* Mixing_Handler::Select(Particle* decayer,
                                             const Hadron_Decay_Table& ot) const
{
  Flavour flav = decayer->Flav();
  string tag = flav.IsAnti() ? flav.Bar().IDName() : flav.IDName();
  if(m_model("Interference_"+tag,0.0)!=0.0) {
    double lifetime = DetermineMixingTime(decayer,false);
    bool anti_at_t0 = decayer->Flav().IsAnti();
    if(decayer->ProductionBlob()->Type()==btp::Hadron_Mixing) anti_at_t0 = !anti_at_t0;
    if(lifetime!=0.0) {
      Hadron_Decay_Table table(ot);
      double cos_term = cos(flav.DeltaM()/rpa->hBar()*lifetime);
      double sin_term = sin(flav.DeltaM()/rpa->hBar()*lifetime);
      double GX, GR, asymmetry, a;
      for(size_t i=0; i<table.size(); i++) {
        Hadron_Decay_Channel* hdc = table.at(i);
        if(hdc->CPAsymmetryS()==0.0 && hdc->CPAsymmetryC()==0.0) continue;
        if(flav.DeltaGamma()==0.0) {
          asymmetry = hdc->CPAsymmetryS()*sin_term - hdc->CPAsymmetryC()*cos_term;
        }
        else {
          Complex lambda = hdc->CPAsymmetryLambda();
          double l2 = sqr(abs(lambda));
          double dG = flav.DeltaGamma()/rpa->hBar();
          asymmetry = (2.0*lambda.imag()/(1.0+l2)*sin_term - (1.0-l2)/(1.0+l2)*cos_term)/
              (cosh(dG*lifetime/2.0) - 2.0*lambda.real()/(1.0+l2)*sinh(dG*lifetime/2.0));
        }
        GX = hdc->Width(); // partial width of this DC
        GR = table.TotalWidth()-GX; // partial width of other DCs
        if(asymmetry>0.0)
          a = -1.0*GR/2.0/GX/asymmetry+sqrt(sqr(GR)/4.0/sqr(GX)/sqr(asymmetry)+(GR+GX)/GX);
        else if(asymmetry<0.0)
          a = -1.0*GR/2.0/GX/asymmetry-sqrt(sqr(GR)/4.0/sqr(GX)/sqr(asymmetry)+(GR+GX)/GX);
        else
          a = 0.0;
        if(anti_at_t0) table.UpdateWidth(hdc, (1.0+a)*GX);
        else           table.UpdateWidth(hdc, (1.0-a)*GX);
      }
      return table.Select();
    }
  }
  return ot.Select();
}
