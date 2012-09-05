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
Inelastic_Event_Generator(Sigma_Inelastic * sigma,
			  list<Omega_ik *> * eikonals,
			  Beam_Remnant_Handler * beams,
			  const int & test) :
  p_sigma(sigma), p_eikonals(eikonals), 
  m_luminosity(Parton_Luminosity(2.,rpa->gen.Ecms())),
  m_laddergenerator(Ladder_Generator(&m_luminosity,test)), p_beams(beams), 
  m_rescatterhandler(Rescatter_Handler(p_beams)), 
  m_Bmin(MBpars("bmin")), m_Bmax(MBpars("bmax")), m_Bsteps(400), 
  m_deltaB((m_Bmax-m_Bmin)/double(m_Bsteps)),
  m_first(true), m_done(false), 
  m_Nladders_fix(MBpars("NLaddersFix")),
  m_kt2fac(MBpars("kt2_factor")), m_difffac(MBpars("diff_factor")),
  m_test(test), m_output(1), m_analyse(m_output>0), p_ladder(NULL)
{ 
  msg_Out()<<METHOD<<" for "<<m_Nladders_fix<<" vs "
	   <<MBpars("NLaddersFix")<<".\n";
  if (m_analyse) {
    m_histograms[string("N_ladder_naive")] = new Histogram(0,0.0,25.0,25);
    m_histograms[string("N_ladder_start")] = new Histogram(0,0.0,25.0,25);
    m_histograms[string("N_ladder_prim")]  = new Histogram(0,0.0,25.0,25);
    m_histograms[string("N_ladder_sec")]   = new Histogram(0,0.0,25.0,25);
    m_histograms[string("N_ladder_true")]  = new Histogram(0,0.0,25.0,25);
    m_histograms[string("B_naive")]        = new Histogram(0,0.0,25.0,50);
    m_histograms[string("B_real")]         = new Histogram(0,0.0,25.0,50);
    m_histograms[string("N_ladder1_B")]    = new Histogram(0,0.0,25.0,25);
    m_histograms[string("N_ladder_all_B")] = new Histogram(0,0.0,25.0,25);
    m_histograms[string("B1_prim")]        = new Histogram(0,0.0,25.0,50);
    m_histograms[string("B1_all")]         = new Histogram(0,0.0,25.0,50);
    m_histograms[string("B2_prim")]        = new Histogram(0,0.0,25.0,50);
    m_histograms[string("B2_all")]         = new Histogram(0,0.0,25.0,50);
  }
  m_connectblobs = m_laddercols = m_updatecols = 0;
  FillGrids();
}

Inelastic_Event_Generator::~Inelastic_Event_Generator() 
{
  msg_Info()<<"In "<<METHOD<<"(out = "<<m_output<<")\n";
  if (m_output) {
    msg_Info()
      <<"Mean number of number of ladders: "
      <<"naive = "<<m_histograms[string("N_ladder_naive")]->Average()<<", "
      <<"start = "<<m_histograms[string("N_ladder_start")]->Average()<<", "
      <<"prim = "<<m_histograms[string("N_ladder_prim")]->Average()<<", "
      <<"true = "<<m_histograms[string("N_ladder_true")]->Average()<<".\n";
    msg_Info()<<"Errors: \n"
	      <<"   Not able to connect blobs "<<m_connectblobs<<";\n"
	      <<"   Wrong colours from ladder "<<m_laddercols<<";\n"
	      <<"   Not able to update colours in event "<<m_updatecols<<".\n";
  }
  Histogram * histo;
  string name;
  for (map<string,Histogram *>::iterator hit=m_histograms.begin();
       hit!=m_histograms.end();hit++) {
    histo = hit->second;
    name  = string("Ladder_Analysis/")+hit->first+string(".dat");
    histo->Finalize();
    histo->Output(name);
    delete histo;
  }
  m_histograms.clear();
}


void Inelastic_Event_Generator::FillGrids() {
  double sigma;
  for (list<Omega_ik *>::iterator eikiter=p_eikonals->begin();
       eikiter!=p_eikonals->end();eikiter++) {
    sigma = p_sigma->Calculate(0.,m_Bmax,(*eikiter));
    msg_Info()<<"Sigma_inel("<<(*eikiter)->FF1()->Number()<<", "
	      <<(*eikiter)->FF2()->Number()<<") = "
	      <<sigma/1.e9<<" mbarn."<<endl;
    (*eikiter)->SetSigmaInelastic(sigma);
    p_sigma->FillGrid(m_Bmin,m_Bmax,m_deltaB,sigma);
  }
  p_sigma->SetSigma();
  p_sigma->SetEikonal(NULL);

  for (size_t beam=0;beam<2;beam++) 
    m_luminosity.SetPDF(p_beams->GetPDF(beam),beam);
  m_luminosity.FillGrids(p_eikonals);
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


int Inelastic_Event_Generator::
InelasticEvent(Blob_List * blobs,const double & xsec,
	       const bool & isUE,const bool & weighted) {
  Blob * blob(blobs->FindFirst(btp::Soft_Collision));
  if (blob && blob->Status()==blob_status::needs_minBias) {
    InitInelasticEvent(isUE,weighted);

    msg_Tracking()<<"-----------------------------------------------------\n"
		  <<METHOD<<"(done = "<<m_done<<", "
		  <<m_Nprim<<" of "<<m_Nladders<<" generated).\n";
  }
  if (m_done) {
    return 0;
  }
  if (m_Nprim<=m_Nladders) {
    switch (AddScatter(blobs,xsec)) {
    case 1:
      return 1;
    case 0:
      blobs->push_front(p_beams->GetSoftColourBlob());
      blobs->SetExternalWeight(xsec);
      m_done = true;
      return 1;
    case -1:
    default:
      break;
    }
  }
  return -1;
}

void Inelastic_Event_Generator::
InitInelasticEvent(const bool & isUE,const bool & weighted) {
  m_isUE     = isUE;
  m_weighted = weighted;
  m_Nprim    = m_Ngen = 0;
  m_first    = true;
  m_done     = false;
  bool success(false);
  do {
    m_B      = p_sigma->FixEikonalAndImpact(p_eikonal);
    int trials(0);
    if (m_analyse) m_histograms[string("B_naive")]->Insert(m_B);
    do {
      if (m_Nladders_fix<=0) 
	m_Nladders = ran->Poissonian((*p_eikonal)(m_B));//+(m_isUE?0:1);
      else m_Nladders = m_Nladders_fix;
      msg_Debugging()<<"   check this: "<<m_B<<" --> "<<m_Nladders<<".\n";
      if (m_analyse) m_histograms[string("N_ladder_naive")]->Insert(m_Nladders);
      if (m_Nladders<1) continue;
      if (m_analyse) m_histograms[string("N_ladder_start")]->Insert(m_Nladders);
      if (p_beams->InitialiseCollision(m_Nladders,p_eikonal)) success = true;
      trials++;
    } while (trials<100 && !success);
  } while (!success);
  if (m_analyse) m_histograms[string("B_real")]->Insert(m_B);
  m_laddergenerator.InitCollision(p_eikonal,m_B);
  m_laddergenerator.SetNPrim(m_Nladders);
  m_rescatterhandler.ResetCollision(p_eikonal,Smin(),m_B);
  m_luminosity.SetEikonal(p_eikonal);
  msg_Tracking()
    <<"######################################################"<<endl
    <<"######################################################"<<endl
    <<"######################################################"<<endl
    <<"######################################################"<<endl
    <<"######################################################"<<endl
    <<METHOD<<" yields "<<m_Nladders
    <<" ladders for B = "<<m_B<<"."<<endl
    <<"######################################################"<<endl
    <<"######################################################"<<endl
    <<"######################################################"<<endl
    <<"######################################################"<<endl
    <<"######################################################"<<endl;
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
      m_histograms[string("B1_all")]->Insert(b1);
      m_histograms[string("B2_all")]->Insert(b2);
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
    m_histograms[string("B1_prim")]->Insert(b1);
    m_histograms[string("B1_all")]->Insert(b1);
    m_histograms[string("B2_prim")]->Insert(b2);
    m_histograms[string("B2_all")]->Insert(b2);
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
  
  Particle * part;
  for (LadderMap::iterator liter=p_ladder->GetEmissionsBegin();
       liter!=p_ladder->GetEmissionsEnd();liter++) {
    part = liter->second.GetParticle();
    blob->AddToOutParticles(part);
  }
  double shat((p_ladder->GetIn1()->GetParticle()->Momentum()
	      +p_ladder->GetIn2()->GetParticle()->Momentum()).Abs2());
  //p_ladder->GetEmissionsBegin()->second.GetParticle()->SetInfo('B');
  //p_ladder->GetEmissionsRBegin()->second.GetParticle()->SetInfo('B');

  m_rescatterhandler.FillInitialStateIntoBlob(blob,p_ladder);
  blob->SetCMS();  

  if (blob->CheckMomentumConservation().Abs2()/shat>1.e-6 ||
      blob->CheckMomentumConservation()[0]/sqrt(shat)>1.e-3 ||
      blob->CheckMomentumConservation()[3]/sqrt(shat)>1.e-3) {
/*    if (blob->CheckMomentumConservation()[0]>1.e-3) {
      msg_Error()<<METHOD<<" try to flip orientation of momenta.\n";
      Vec4D mom;
      for (int i=0;i<blob->NOutP();i++) {
	mom = blob->OutParticle(i)->Momentum();
	mom[3] *= -1.;
	blob->OutParticle(i)->SetMomentum(mom);
      }
    }*/
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
  msg_Tracking()<<METHOD<<":\n"<<(*blob)<<"\n"
		<<"===============================================\n";
  return true;
}


double Inelastic_Event_Generator::Smin() const {
  double smin(m_luminosity.Smin());
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
  int	  nval(10000);
  double  mcmean(0.),anamean(0.),a,c,value;
  double  Y(p_sigma->Y());
  double  kappa(eikonal->Kappa_i());
  double  beta0(eikonal->FF1()->Beta0());
  double  Lambda2(eikonal->Lambda2());
  double  Delta(eikonal->Delta());
  a = Lambda2/(8.*(1.+kappa));
  c = sqr(beta0)*Lambda2*(1.+kappa)*exp(2.*Delta*Y)/(8.*M_PI);
  anamean = c*exp(-a*sqr(B));
  for(int i=0; i<nval; i++){
    value = double(ran->Poissonian((*eikonal)(B)));
    mcmean += value/nval;
  }
  msg_Tracking()<<"In "<< METHOD <<" mean number of ladders: "<<endl
	   << "		"<<mcmean<<" (Monte Carlo); " 
	   <<(*eikonal)(B)<<" (eikonal); "<<anamean<<" (analytic)"<<endl;
}


