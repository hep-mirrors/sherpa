#include "SHRiMPS/Event_Generation/Inelastic_Event_Generator.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Inelastic_Event_Generator::
Inelastic_Event_Generator(Beam_Remnant_Handler * beams) :
  m_luminosity(Parton_Luminosity(beams)),
  m_laddergenerator(Ladder_Generator(&m_luminosity,0)),
  p_beams(beams), 
  m_rescatterhandler(Rescatter_Handler(p_beams)), 
  m_output(1), p_ladder(NULL)
{
  FillCrossSections();
}

Inelastic_Event_Generator::~Inelastic_Event_Generator() {
  Reset();
}

void Inelastic_Event_Generator::FillCrossSections() {
  m_xsec = m_sigma.Calculate();
  msg_Info()<<METHOD<<" yields inelastic cross section "
	    <<"sigma = "<<m_xsec/1.e9<<" mbarn.\n";
  m_sigma.FillDifferentialGrids();
  m_luminosity.FillGrids(this);
}

void Inelastic_Event_Generator::Reset() {
  delete p_ladder;
  p_ladder = NULL;
}

bool Inelastic_Event_Generator::DressShowerBlob(Blob * blob) {
  msg_Error()<<METHOD<<" not implemented for blob "
	     <<"["<<blob->Id()<<", "<<blob->Type()<<"].\n";
  abort();
}


int Inelastic_Event_Generator::GenerateEvent(Blob_List * blobs,
					     const bool & isUE) {
  msg_Info()<<"*** -----------------------------------------------------\n"
	    <<"*** "<<METHOD<<"(done = "<<m_done<<").\n";
  Blob * blob(blobs->FindFirst(btp::Soft_Collision));
  if (blob && blob->Status()==blob_status::needs_minBias) {
    InitInelasticEvent();
    msg_Info()<<"*"<<METHOD<<"("<<m_Nprim<<" of "<<m_Nladders<<").\n";
  }
  if (m_done) return 0;
  if (m_Nprim<=m_Nladders) {
    switch (AddScatter(blobs,m_sigma.XSec())) {
    case 1:
      return 1;
    case 0:
      blobs->push_front(p_beams->GetSoftColourBlob());
      blobs->SetExternalWeight(m_sigma.XSec());
      m_done = true;
      return 1;
    case -1:
    default:
      break;
    }
  }
  return -1;
}

void Inelastic_Event_Generator::InitInelasticEvent() {
  m_Nprim = m_Ngen = 0;
  m_first = true;
  m_done  = false;
  for (size_t trials=0;trials<1000;trials++) {
    FixEikonalAndImpact();    
    m_Nladders = ran->Poissonian((*p_eikonal)(m_B));
    msg_Out()<<"*** "<<METHOD<<"(B = "<<m_B<<") -> select N_ladders with "
	     <<(*p_eikonal)(m_B)<<" --> "<<m_Nladders<<".\n";
    if (m_Nladders>1) break;
    //&&
    //  p_beams->InitialiseCollision(m_Nladders,m_B,p_eikonal)) break;
  }
  msg_Info()<<"*** "<<METHOD<<" selected B = "<<m_B
	    <<" and eikonal ["<<p_eikonal->FF1()->Number()
	    <<p_eikonal->FF2()->Number()<<"]\n"
	    <<"***  ---> "<<m_Nladders<<" ladders to be generated.\n";
  //m_laddergenerator.InitCollision(p_eikonal,m_B);
  //m_laddergenerator.SetNPrim(m_Nladders);
  //m_rescatterhandler.ResetCollision(p_eikonal,Smin(),m_B);
  //m_luminosity.SetEikonal(p_eikonal);
  exit(1);
}

void Inelastic_Event_Generator::FixEikonalAndImpact() {
  p_eikonal = m_sigma.SelectEikonal();
  m_B = m_sigma.SelectB();
}


int Inelastic_Event_Generator::
AddScatter(Blob_List * blobs,const double & xsec) 
{
  msg_Tracking()<<METHOD<<"("<<m_Nprim<<" from "<<m_Nladders<<"):\n";
  Particle * part1(NULL), * part2(NULL);
  if (!m_rescatterhandler.ConnectBlobs(blobs,p_beams->GetCompensatorBlob())) {
    m_connectblobs++;
    return -1;
  }
  if (m_Nprim>0 && m_Nprim<=m_Nladders) {
    m_rescatterhandler.UpdateCollision(blobs);
    if (m_rescatterhandler.SelectRescatter(part1,part2)) {
      m_Nsec++;
      //m_laddergenerator.SetLadderGeneration(1+m_Nsec);
      p_ladder = m_laddergenerator(part1,part2,true);
      if (!p_ladder) return -1;
      p_beams->SetInitials(part1,part2);
      if (!p_beams->UpdateColours(p_ladder,false)) {
	m_updatecols++;
	return -1;
      }
      if (!CreateBlob(blobs,1.)) {
	m_laddercols++;
	return -1;
      }
      m_rescatterhandler.Map(part1,p_ladder->GetIn1()->GetParticle());
      m_rescatterhandler.Map(part2,p_ladder->GetIn2()->GetParticle());
      const double b1(m_laddergenerator.B1()), b2(m_laddergenerator.B2());
      if (m_analyse) {
	m_histograms[string("B1_all")]->Insert(b1);
	m_histograms[string("B2_all")]->Insert(b2);
      }
      m_Ngen++;
      return 1;
    }
  }
  if (m_Nprim<m_Nladders) {
    if (!p_beams->NextIS(part1,part2)) return -1;
    m_Nsec = 0;
    //m_laddergenerator.SetLadderGeneration(1);
    p_ladder = m_laddergenerator(part1,part2,false,m_Nprim==0,m_weighted);
    if (!p_ladder) return -1;
    if (!p_beams->UpdateColours(p_ladder,m_Nprim+1==m_Nladders)) {
      m_updatecols++;
      return -1;
    }
    m_rescatterhandler.ResetRescatter(m_Nprim==0);
    if (!CreateBlob(blobs,m_weighted?xsec*p_ladder->Weight():1.)) {
      m_laddercols++;
      return -1;
    }
    const double b1(m_laddergenerator.B1()), b2(m_laddergenerator.B2());
    if (m_analyse) {
      m_histograms[string("B1_prim")]->Insert(b1);
      m_histograms[string("B1_all")]->Insert(b1);
      m_histograms[string("B2_prim")]->Insert(b2);
      m_histograms[string("B2_all")]->Insert(b2);
    }
    m_Nprim++; m_Ngen++;
    return 1;
  }
  else {
    msg_Tracking()
      <<"##################################################################\n"
      <<"   Out of event with "<<m_Nprim<<"/"<<m_Ngen<<" from "<<m_Nladders
      <<"\n"//<<(*blobs)<<"\n"
      <<"##################################################################\n";

    if (m_analyse) {
      m_histograms[string("N_ladder_prim")]->Insert(m_Nprim);
      m_histograms[string("N_ladder_true")]->Insert(m_Ngen);
      m_histograms[string("N_ladder_sec")]->Insert(m_Ngen-m_Nprim);
      m_histograms[string("N_ladder1_B")]->Insert(m_B,m_Nprim);
      m_histograms[string("N_ladder_all_B")]->Insert(m_B,m_Ngen);
    }
    return 0;
  }
  msg_Tracking()<<"   @@@ undefined, Nprim = "<<m_Nprim<<", return -1.\n";
  return -1;
}



bool Inelastic_Event_Generator::
CreateBlob(Blob_List * blobs,const double & xsec) {
  Vec4D pos(p_ladder->Position()*rpa->hBar()*rpa->c());
  Blob * blob(blobs->FindFirst(btp::Soft_Collision));
  //msg_Out()<<METHOD<<"("<<blobs->size()<<", blob = "<<blob
  //	   <<", status = "<<blob->Status()<<"): \n";
  //if (blobs->size()<=2) msg_Out()<<(*blobs)<<"\n";
  bool add(true);
  if (blob && blob->Status()==blob_status::needs_minBias) {
    if (blob->NInP()>0)  {
      msg_Error()<<"Error in "<<METHOD<<": blob has particles."<<endl
		 <<(*blob)<<endl;
      blob->DeleteInParticles();
    }
    if (blob->NOutP()>0) {
      msg_Error()<<"Error in "<<METHOD<<": blob has particles."<<endl
		 <<(*blob)<<endl;
      blob->DeleteOutParticles();
    }    
    blob->UnsetStatus(blob_status::needs_minBias);
    add = false;
  }
  else {
    blob = new Blob();
    blob->SetId();
    blobs->push_back(blob);
  }

  blob->SetType(btp::Hard_Collision);
  if (m_isUE) blob->SetTypeSpec("UnderlyingEvent");
         else blob->SetTypeSpec("MinBias");    
  blob->SetStatus(blob_status::needs_showers);
  blob->SetPosition(pos);
  Blob_Data_Base *winfo((*blob)["Weight"]);
  if (!winfo) blob->AddData("Weight",new ATOOLS::Blob_Data<double>(1.));
  Blob_Data_Base *wninfo((*blob)["Weight_Norm"]);
  if (!wninfo) blob->AddData("Weight_Norm",new ATOOLS::Blob_Data<double>(1.));
  Blob_Data_Base *tinfo((*blob)["Trials"]);
  if (!tinfo) blob->AddData("Trials",new ATOOLS::Blob_Data<double>(1.));

  Particle * part;
  for (LadderMap::iterator liter=p_ladder->GetEmissionsBegin();
       liter!=p_ladder->GetEmissionsEnd();liter++) {
    part = liter->second.GetParticle();
    blob->AddToOutParticles(part);
  }
  double shat((p_ladder->GetIn1()->GetParticle()->Momentum()
	      +p_ladder->GetIn2()->GetParticle()->Momentum()).Abs2());
  m_rescatterhandler.FillInitialStateIntoBlob(blob,p_ladder);
  blob->SetCMS();  

  if (blob->CheckMomentumConservation().Abs2()/shat>1.e-6 ||
      blob->CheckMomentumConservation()[0]/sqrt(shat)>1.e-3 ||
      blob->CheckMomentumConservation()[3]/sqrt(shat)>1.e-3) {
    msg_Error()<<"Problem in "<<METHOD<<":\n"
	       <<"   Scattering blob ("<<blob->Id()<<") seems fishy: "
	       <<blob->CheckMomentumConservation()<<".\n"
	       <<(*blob)<<"\n"<<(*p_ladder)<<"\n";
  }
  if (!blob->CheckColour()) {
    msg_Error()<<"Problem in "<<METHOD<<":\n"
	       <<"   Scattering blob ("<<blob->Id()<<") seems fishy: "
	       <<"Bad colour configuration.\n"
	       <<(*blob)<<"\n"<<(*p_ladder)<<"\n";
    return false;
  }
  //msg_Out()<<METHOD<<":\n"<<(*blob)<<"\n"
  //	   <<"--> Hand over to shower now.\n"
  //	   <<"===============================================\n";
  return true;
}


double Inelastic_Event_Generator::Smin() const {
  double smin(m_luminosity.Smin()*m_Nladders);
  if (!p_ladder) return smin;
  //   smin *= m_kt2fac;
  if (p_ladder->IsHardDiffractive() && p_ladder->Size()==2) smin *= m_difffac;
  return smin; 
}

bool Inelastic_Event_Generator::IsLastRescatter() const {
  if (!p_ladder) return false;
  return p_ladder->IsRescatter();
}
/////////////////////////////////////////////////////////////////////////

void Inelastic_Event_Generator::
TestNumberOfLadders(Omega_ik * eikonal,const double & B){
  // int	  nval(10000);
  // double  mcmean(0.),anamean(0.),a,c,value;
  // double  Y(p_sigma->Y());
  // double  kappa(eikonal->Kappa_i());
  // double  beta0(eikonal->FF1()->Beta0());
  // double  Lambda2(eikonal->Lambda2());
  // double  Delta(eikonal->Delta());
  // a = Lambda2/(8.*(1.+kappa));
  // c = sqr(beta0)*Lambda2*(1.+kappa)*exp(2.*Delta*Y)/(8.*M_PI);
  // anamean = c*exp(-a*sqr(B));
  // for(int i=0; i<nval; i++){
  //   value = double(ran->Poissonian((*eikonal)(B)));
  //   mcmean += value/nval;
  // }
  // msg_Tracking()<<"In "<< METHOD <<" mean number of ladders: "<<endl
  // 	   << "		"<<mcmean<<" (Monte Carlo); " 
  // 	   <<(*eikonal)(B)<<" (eikonal); "<<anamean<<" (analytic)"<<endl;
}


