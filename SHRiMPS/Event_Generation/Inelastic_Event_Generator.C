#include "SHRiMPS/Main/Cluster_Algorithm.H"
#include "SHRiMPS/Event_Generation/Inelastic_Event_Generator.H"
#include "SHRiMPS/Cross_Sections/Sigma_Inelastic.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Inelastic_Event_Generator::Inelastic_Event_Generator()
{
  Initialise();
}

Inelastic_Event_Generator::~Inelastic_Event_Generator() {
  Reset();
  while (!m_Bgrids.empty()) {
    delete m_Bgrids.begin()->second;
    m_Bgrids.erase(m_Bgrids.begin());
  }
  m_Bgrids.clear();
}

void Inelastic_Event_Generator::Initialise() {
  m_sigma = 0.;
  Sigma_Inelastic sigma;
  list<Omega_ik *> * eikonals(MBpars.GetEikonals());
  for (list<Omega_ik *>::iterator eikonal=eikonals->begin();
       eikonal!=eikonals->end(); eikonal++) {
    m_Bgrids[(*eikonal)] = sigma.FillBGrid((*eikonal));
    m_sigma += m_xsecs[(*eikonal)] =
      m_Bgrids[(*eikonal)]->back()*rpa->Picobarn();
  }
  msg_Info()<<METHOD<<" yields effective inelastic cross section "
	    <<"sigma = "<<m_sigma/1.e9<<" mbarn.\n";
  Reset();
}

void Inelastic_Event_Generator::Reset() {
  m_init = m_done = false;
  m_first = true;
}

int Inelastic_Event_Generator::
GenerateEvent(Blob_List * blobs,const bool & isUE) {
  if (m_done || !InitInelasticEvent(blobs)) return 0;
  if (!AddScatter(blobs)) {
    m_done = true;
    return 0;
  }
  if (m_Ngen++ > m_Nladders) m_done = true;
  return 1;
}

bool Inelastic_Event_Generator::InitInelasticEvent(Blob_List * blobs) {
  if (m_init) return true;
  Blob * blob(blobs->FindFirst(btp::Soft_Collision));
  if (!(blob && blob->Status()==blob_status::needs_minBias) ||
      !SetUpEvent()) return false;
  m_init = true;
  m_laddergenerator.InitCollision(p_eikonal,m_B,m_Nladders);
  return true;
}

bool Inelastic_Event_Generator::SetUpEvent() {
  p_eikonal  = 0; m_B = -1;
  m_Nladders = m_Nprim = m_Ngen = 0;
  size_t trials=0;
  do {
    SelectEikonal();
    SelectB();
    m_Nladders = ran->Poissonian(2.*(*p_eikonal)(m_B));
  } while (trials++<1000 && m_Nladders<1);
  return (trials<1000);
}

bool Inelastic_Event_Generator::SelectEikonal() {
  p_eikonal = 0;
  while (p_eikonal==NULL) {
    double disc = ran->Get()*m_sigma;
    for (std::map<Omega_ik *,double>::iterator eikiter=m_xsecs.begin();
	 eikiter!=m_xsecs.end();eikiter++) {
      disc-=eikiter->second;
      if (disc<=1.e-12) {
	p_eikonal = eikiter->first;
	break;
      }
    }
  }
  return (p_eikonal!=0);
}

bool Inelastic_Event_Generator::SelectB() {
  if (p_eikonal==0) {
    msg_Error()<<"Error in "<<METHOD<<": no eikonal selected.\n";
    return false;
  }
  std::vector<double> * grid = m_Bgrids[p_eikonal];  
  double deltaB(p_eikonal->DeltaB());
  m_B = -1.;
  do {
    double random = ran->Get()*(*grid)[grid->size()-1];
    size_t bin(0);
    while (bin<grid->size()-1 && (random-(*grid)[bin]>=0)) bin++;
    if (bin>=grid->size()) continue;
    double inthigh((*grid)[bin]), intlow((*grid)[bin-1]);
    double Bhigh(bin*deltaB), Blow((bin-1)*deltaB);
    m_B  = (Blow*(random-intlow)+Bhigh*(inthigh-random))/(inthigh-intlow);
  } while (m_B<0.);
  return (m_B>=0.);
}

int Inelastic_Event_Generator::AddScatter(Blob_List * blobs) {
  Blob * blob(blobs->FindFirst(btp::Soft_Collision));
  if (!(blob && blob->Status()==blob_status::needs_minBias)) return -1;
  if (m_Ngen>0) blob = CreateBlob();
  SetBlobType(blob);
  if (!m_laddergenerator.MakePrimaryLadder(blob,m_Ngen==0)) {
    delete blob;
    return 0;
  }
  if (m_Ngen!=0) blobs->push_back(blob);
  p_cluster->SetTMax(4.*m_laddergenerator.TMax());
  p_cluster->SetMinKT2(64.);
  //msg_Out()<<METHOD<<"["<<p_cluster<<"]: tmax = "
  //	   <<m_laddergenerator.TMax()<<".\n";
  if (m_Ngen==0) blobs->SetExternalWeight(m_sigma);
  return 1;
}

Blob * Inelastic_Event_Generator::CreateBlob() {
  Blob * blob = new Blob();
  blob->SetId();
  return blob;
}

void Inelastic_Event_Generator::SetBlobType(Blob * blob) {
  blob->SetType(btp::Hard_Collision);
  blob->SetTypeSpec("MinBias");    
  blob->SetStatus(blob_status::needs_showers|blob_status::needs_beams);
}
  
void Inelastic_Event_Generator::Test(const std::string & dirname) {
  TestSelectB(dirname);
  TestNumberOfLadders(dirname);
  m_laddergenerator.Test(dirname);
}

void Inelastic_Event_Generator::TestSelectB(const std::string & dirname) {
  Histogram histo_in(0,0.,MBpars.GetEikonalParameters().bmax,100);
  Histogram histo_out(0,0.,MBpars.GetEikonalParameters().bmax,100);
  SelectEikonal();
  Sigma_Inelastic sigma;
  sigma.SetEikonal(p_eikonal);
  double deltaB(p_eikonal->DeltaB()), B(deltaB/2.), val;
  do {
    val = 2.*M_PI*B*sigma.GetValue(B);
    histo_in.Insert(B,val);
    B  += deltaB;
  } while (B<MBpars.GetEikonalParameters().bmax);
  for (int i=0;i<1000000;i++) {
    SelectB();
    histo_out.Insert(m_B);
  }
  histo_in.Finalize();
  histo_in.Output(dirname+"/B_grid.dat");
  histo_out.Finalize();
  histo_out.Output(dirname+"/B_distribution.dat");
}

void Inelastic_Event_Generator::
TestNumberOfLadders(const std::string & dirname) {
  for (size_t i=1;i<7;i++) {
    double arg(2.*(*p_eikonal)(double(i)));
    msg_Info()<<METHOD<<"; eik("<<i<<") = "<<arg<<".\n";
    Histogram histo(0,0.,100.0,100);
    for (int trials=0;trials<1000000;trials++)
      histo.Insert(double(ran->Poissonian(arg))+0.5);
    ostringstream converter;
    converter << i;
    string istr(converter.str());
    string name(string("NLadders_B_")+istr+string("_GeV.dat"));
    histo.Finalize();
    histo.Output(dirname+"/"+name);
  }
}

  
// bool Inelastic_Event_Generator::DressShowerBlob(Blob * blob) {
//   msg_Error()<<METHOD<<" not implemented for blob "
// 	     <<"["<<blob->Id()<<", "<<blob->Type()<<"].\n";
//   Abort();
// }

// void Inelastic_Event_Generator::FixEikonalAndImpact() {
//   p_eikonal = m_sigma.SelectEikonal();
//   m_B = m_sigma.SelectB();
// }


// int Inelastic_Event_Generator::
// AddScatter(Blob_List * blobs,const double & xsec) 
// {
//   msg_Tracking()<<METHOD<<"("<<m_Nprim<<" from "<<m_Nladders<<"):\n";
//   Particle * part1(NULL), * part2(NULL);
//   if (!m_rescatterhandler.ConnectBlobs(blobs,p_beams->GetCompensatorBlob())) {
//     m_connectblobs++;
//     return -1;
//   }
//   if (m_Nprim>0 && m_Nprim<=m_Nladders) {
//     m_rescatterhandler.UpdateCollision(blobs);
//     if (m_rescatterhandler.SelectRescatter(part1,part2)) {
//       m_Nsec++;
//       //m_laddergenerator.SetLadderGeneration(1+m_Nsec);
//       p_ladder = m_laddergenerator(part1,part2,true);
//       if (!p_ladder) return -1;
//       p_beams->SetInitials(part1,part2);
//       if (!p_beams->UpdateColours(p_ladder,false)) {
// 	m_updatecols++;
// 	return -1;
//       }
//       if (!CreateBlob(blobs,1.)) {
// 	m_laddercols++;
// 	return -1;
//       }
//       m_rescatterhandler.Map(part1,p_ladder->GetIn1()->GetParticle());
//       m_rescatterhandler.Map(part2,p_ladder->GetIn2()->GetParticle());
//       const double b1(m_laddergenerator.B1()), b2(m_laddergenerator.B2());
//       if (m_analyse) {
// 	m_histograms[string("B1_all")]->Insert(b1);
// 	m_histograms[string("B2_all")]->Insert(b2);
//       }
//       m_Ngen++;
//       return 1;
//     }
//   }
//   if (m_Nprim<m_Nladders) {
//     if (!p_beams->NextIS(part1,part2)) return -1;
//     m_Nsec = 0;
//     //m_laddergenerator.SetLadderGeneration(1);
//     p_ladder = m_laddergenerator(part1,part2,false,m_Nprim==0,m_weighted);
//     if (!p_ladder) return -1;
//     if (!p_beams->UpdateColours(p_ladder,m_Nprim+1==m_Nladders)) {
//       m_updatecols++;
//       return -1;
//     }
//     m_rescatterhandler.ResetRescatter(m_Nprim==0);
//     if (!CreateBlob(blobs,m_weighted?xsec*p_ladder->Weight():1.)) {
//       m_laddercols++;
//       return -1;
//     }
//     const double b1(m_laddergenerator.B1()), b2(m_laddergenerator.B2());
//     if (m_analyse) {
//       m_histograms[string("B1_prim")]->Insert(b1);
//       m_histograms[string("B1_all")]->Insert(b1);
//       m_histograms[string("B2_prim")]->Insert(b2);
//       m_histograms[string("B2_all")]->Insert(b2);
//     }
//     m_Nprim++; m_Ngen++;
//     return 1;
//   }
//   else {
//     msg_Tracking()
//       <<"##################################################################\n"
//       <<"   Out of event with "<<m_Nprim<<"/"<<m_Ngen<<" from "<<m_Nladders
//       <<"\n"//<<(*blobs)<<"\n"
//       <<"##################################################################\n";

//     if (m_analyse) {
//       m_histograms[string("N_ladder_prim")]->Insert(m_Nprim);
//       m_histograms[string("N_ladder_true")]->Insert(m_Ngen);
//       m_histograms[string("N_ladder_sec")]->Insert(m_Ngen-m_Nprim);
//       m_histograms[string("N_ladder1_B")]->Insert(m_B,m_Nprim);
//       m_histograms[string("N_ladder_all_B")]->Insert(m_B,m_Ngen);
//     }
//     return 0;
//   }
//   msg_Tracking()<<"   @@@ undefined, Nprim = "<<m_Nprim<<", return -1.\n";
//   return -1;
// }



// bool Inelastic_Event_Generator::
// CreateBlob(Blob_List * blobs,const double & xsec) {
//   Vec4D pos(p_ladder->Position()*rpa->hBar()*rpa->c());
//   Blob * blob(blobs->FindFirst(btp::Soft_Collision));
//   //msg_Out()<<METHOD<<"("<<blobs->size()<<", blob = "<<blob
//   //	   <<", status = "<<blob->Status()<<"): \n";
//   //if (blobs->size()<=2) msg_Out()<<(*blobs)<<"\n";
//   bool add(true);
//   if (blob && blob->Status()==blob_status::needs_minBias) {
//     if (blob->NInP()>0)  {
//       msg_Error()<<"Error in "<<METHOD<<": blob has particles."<<endl
// 		 <<(*blob)<<endl;
//       blob->DeleteInParticles();
//     }
//     if (blob->NOutP()>0) {
//       msg_Error()<<"Error in "<<METHOD<<": blob has particles."<<endl
// 		 <<(*blob)<<endl;
//       blob->DeleteOutParticles();
//     }    
//     blob->UnsetStatus(blob_status::needs_minBias);
//     add = false;
//   }
//   else {
//     blob = new Blob();
//     blob->SetId();
//     blobs->push_back(blob);
//   }

//   blob->SetType(btp::Hard_Collision);
//   if (m_isUE) blob->SetTypeSpec("UnderlyingEvent");
//          else blob->SetTypeSpec("MinBias");    
//   blob->SetStatus(blob_status::needs_showers);
//   blob->SetPosition(pos);
//   Blob_Data_Base *winfo((*blob)["Weight"]);
//   if (!winfo) blob->AddData("Weight",new ATOOLS::Blob_Data<double>(1.));
//   Blob_Data_Base *wninfo((*blob)["Weight_Norm"]);
//   if (!wninfo) blob->AddData("Weight_Norm",new ATOOLS::Blob_Data<double>(1.));
//   Blob_Data_Base *tinfo((*blob)["Trials"]);
//   if (!tinfo) blob->AddData("Trials",new ATOOLS::Blob_Data<double>(1.));

//   Particle * part;
//   for (LadderMap::iterator liter=p_ladder->GetEmissionsBegin();
//        liter!=p_ladder->GetEmissionsEnd();liter++) {
//     part = liter->second.GetParticle();
//     blob->AddToOutParticles(part);
//   }
//   double shat((p_ladder->GetIn1()->GetParticle()->Momentum()
// 	      +p_ladder->GetIn2()->GetParticle()->Momentum()).Abs2());
//   m_rescatterhandler.FillInitialStateIntoBlob(blob,p_ladder);
//   blob->SetCMS();  

//   if (blob->CheckMomentumConservation().Abs2()/shat>1.e-6 ||
//       blob->CheckMomentumConservation()[0]/sqrt(shat)>1.e-3 ||
//       blob->CheckMomentumConservation()[3]/sqrt(shat)>1.e-3) {
//     msg_Error()<<"Problem in "<<METHOD<<":\n"
// 	       <<"   Scattering blob ("<<blob->Id()<<") seems fishy: "
// 	       <<blob->CheckMomentumConservation()<<".\n"
// 	       <<(*blob)<<"\n"<<(*p_ladder)<<"\n";
//   }
//   if (!blob->CheckColour()) {
//     msg_Error()<<"Problem in "<<METHOD<<":\n"
// 	       <<"   Scattering blob ("<<blob->Id()<<") seems fishy: "
// 	       <<"Bad colour configuration.\n"
// 	       <<(*blob)<<"\n"<<(*p_ladder)<<"\n";
//     return false;
//   }
//   //msg_Out()<<METHOD<<":\n"<<(*blob)<<"\n"
//   //	   <<"--> Hand over to shower now.\n"
//   //	   <<"===============================================\n";
//   return true;
// }


// double Inelastic_Event_Generator::Smin() const {
//   double smin(m_luminosity.Smin()*m_Nladders);
//   if (!p_ladder) return smin;
//   //   smin *= m_kt2fac;
//   if (p_ladder->IsHardDiffractive() && p_ladder->Size()==2) smin *= m_difffac;
//   return smin; 
// }

// bool Inelastic_Event_Generator::IsLastRescatter() const {
//   if (!p_ladder) return false;
//   return p_ladder->IsRescatter();
// }
/////////////////////////////////////////////////////////////////////////


//void Inelastic_Event_Generator::
//TestNumberOfLadders(Omega_ik * eikonal,const double & B){
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
//}


