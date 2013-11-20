#include "SHRiMPS/Event_Generation/Ladder_Generator.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "SHRiMPS/Event_Generation/Quark_Replace.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Ladder_Generator::
Ladder_Generator(Parton_Luminosity * lumi,const int & test) :
  m_IS(Initial_State(lumi)), m_FS(Final_State(test)),
  m_originalY(MBpars("originalY")), m_cutoffY(MBpars("deltaY")), 
  p_ladder(0), m_output(false),
  m_Nprim(0),m_Nsec(0),m_Ndd_p(0),m_Ndd_s(0),m_Nsd_p(0),m_Nsd_s(0),
  m_Ncep_p(0),m_Ncep_s(0)
{    
  Vec4D cms(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
  m_pplus   = cms.PPlus();
  m_pminus  = cms.PMinus();

  if (m_output) {
    m_N1h = m_N1 = m_N8 = m_resc0 = m_resc1 = 0;
    m_wtover1 = m_wtover2 = 0.;  
    m_wt = 0;
    m_histograms[string("N_inel")]    = new Histogram(0,0.0,20.0,20);
    m_histograms[string("N_diff")]    = new Histogram(0,0.0,20.0,20);
    m_histograms[string("N_hdiff")]   = new Histogram(0,0.0,20.0,20);
    m_histograms[string("Nemit1")]    = new Histogram(0,0.0,10.0,10);
    m_histograms[string("Nemit2")]    = new Histogram(0,0.0,10.0,10);
    m_histograms[string("Delta1")]    = new Histogram(0,0.0,5.0,50);
    m_histograms[string("Delta2")]    = new Histogram(0,0.0,5.0,50);
    m_histograms[string("D_y1")]      = new Histogram(0,-10.0,10.0,100);
    m_histograms[string("D_y2")]      = new Histogram(0,-10.0,10.0,100);
    m_histograms[string("KT1")]       = new Histogram(0, 0.0,50.0,200);
    m_histograms[string("KT1mid")]    = new Histogram(0, 0.0,50.0,200);
    m_histograms[string("KT2")]       = new Histogram(0, 0.0,50.0,200);
    m_histograms[string("KT2mid")]    = new Histogram(0, 0.0,50.0,200);
    m_histograms[string("y1")]        = new Histogram(0,-10.0,10.0,100);
    m_histograms[string("LadderWt")]  = new Histogram(0,0.0,2.0,200);
    m_histograms[string("QT")]        = new Histogram(0, 0.0,10.0,200);
    m_histograms[string("DDmass_p")]  = new Histogram(0, 0.0,1000.0,100);
    m_histograms[string("SDmass_p")]  = new Histogram(0, 0.0,10.0,100);
    m_histograms[string("CEPmass_p")] = new Histogram(0, 0.0,10.0,100);
    m_histograms[string("DDmass_s")]  = new Histogram(0, 0.0,1000.0,100);
    m_histograms[string("SDmass_s")]  = new Histogram(0, 0.0,10.0,100);
    m_histograms[string("CEPmass_s")] = new Histogram(0, 0.0,10.0,100);
    
    m_histogram2ds[string("pt2vspt2othat")] = 
      new Histogram_2D(0, 0.0, 100, 100, 0.0, 2.0, 20);
  }
}

Ladder_Generator::~Ladder_Generator() {
  if (m_output) {
    msg_Info()
      <<METHOD<<":\n"
      <<"   Ratio of diffractive ladders: "
      <<double(m_N1)/double(m_N1+m_N8)<<", of which hard diffractive = "
      <<(100.*m_N1h)/m_N1<<"%;\n"
      <<"   mean number of extra emissions in primary ladders: "
      <<(m_histograms[string("Nemit1")]->Average()-2)<<", "
      <<"Delta = "<<m_histograms[string("Delta1")]->Average()<<", "
      <<" average kt = "<<m_histograms[string("KT1")]->Average()
      <<"(mid-y:"<<m_histograms[string("KT1mid")]->Average()<<");\n";
    if (MBpars.RescMode()!=resc_mode::off &&
	m_histograms[string("Nemit2")]->Integral()>1.e-6) 
      msg_Info()<<"   mean number of extra emissions in secondary ladders: "
		<<(m_histograms[string("Nemit2")]->Average()-2)<<", "
		<<"Delta = "<<m_histograms[string("Delta2")]->Average()<<", "
		<<" average kt = "<<m_histograms[string("KT2")]->Average()
		<<"(mid-y:"<<m_histograms[string("KT2mid")]->Average()<<");\n";
    if (m_resc1>0) 
      msg_Info()<<"   had to enforce "<<(m_resc0/(m_resc1+m_resc0))<<" "
		<<"secondary ladders to be singlets.\n";
    if (m_wt>0)
      msg_Info()<<"    overshooting weights: "
		<<(100.*(m_wtover1+m_wtover2)/m_wt)<<" % "
		<<" in total, "<<(1.*m_wtover2/(m_wtover1+m_wtover2)*100)
		<<" % in secondary ladders, in total "<<m_wt<<" calls.\n";
  
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator 
	   hit=m_histograms.begin();hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = string("Ladder_Analysis/")+hit->first+string(".dat");
      histo->Finalize();
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  
    Histogram_2D * histo2d;
    for (map<string,Histogram_2D *>::iterator 
	   hit=m_histogram2ds.begin();hit!=m_histogram2ds.end();hit++) {
      histo2d = hit->second;
      name  = string("Ladder_Analysis/")+hit->first+string(".dat");
      histo2d->Finalize();
      histo2d->Output(name);
      delete histo2d;
    }
    m_histograms.clear();
  }
  if (p_ladder) delete p_ladder; 
}

void Ladder_Generator::
InitCollision(Omega_ik * eikonal,const double & B) {
  p_eikonal = eikonal;
  m_B       = B;
  m_IS.InitNewCollision(p_eikonal,m_B);
}


Ladder * Ladder_Generator::
operator()(Particle * part1,Particle * part2,const bool & rescatter,
	   const bool & first,const bool & weighted)
{
  double weight(0.),isweight(0.);  
  int trials(0);
  
  do {
    if (p_ladder) { delete p_ladder; p_ladder = NULL; }
    if (trials==1000) {
      msg_Tracking()<<METHOD<<" needs too many trials.  Return 0.\n";
      return p_ladder; 
    }
    isweight = weight  = InitialiseLadder(part1,part2,rescatter);
    if (p_ladder->IsRescatter()) {
      p_ladder->SetKtCut2(m_prodkt2/double(m_Nprops));
    }
    else 
      p_ladder->SetKtCut2(0.);
    weight  *= 
      //Max(0.,m_FS(p_ladder,m_cutoffY,first,trials==0))*
      Max(0.,m_FS(p_ladder,0.,first,trials==0))*
      Weight(isweight);
  } while ((trials++)<1000 && weight<ran->Get());
  if (p_ladder->Size()<2) {
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   Ladder with size = "<<p_ladder->Size()<<" survives.\n"
	       <<"   Return NULL and hope for the best.\n";
    return NULL;
  }

  //Quark_Replace doquarks;
  //doquarks.ReplaceWithQuarks(p_ladder);
  
  //check four momentum conservation
  if (!p_ladder->CheckFourMomentum()) 
    msg_Out()<<METHOD<<" Four Momentum violation in ladder"<<std::endl
             <<(*p_ladder)<<std::endl;
  
  if (!p_ladder->IsRescatter()) {
    m_prodkt2 = 1.;
    m_Nprops  = 0;
  }
  p_ladder->ProdProps(m_prodkt2,m_Nprops);


  if (m_output) Analyse(!p_ladder->IsRescatter());
  msg_Tracking()<<(*p_ladder);
  msg_Tracking()
    <<"    ---> accepted ladder with total weight="<<weight<<".\n"
    <<"   ------------------------------------------------------\n";
  return p_ladder;
}

double Ladder_Generator::
InitialiseLadder(Particle * part1,Particle * part2,const bool & rescatter) {
  Ladder_Particle * lpart1(new Ladder_Particle(part1));
  Ladder_Particle * lpart2(new Ladder_Particle(part2));
  m_IS.DefineIS(lpart1,lpart2,rescatter); 
  m_FS.Init(p_eikonal,m_IS.B1(),m_IS.B2(),m_IS.B());
  p_ladder = new Ladder(m_IS.Pos());
  p_ladder->SetInParticles(lpart1,lpart2);
  p_ladder->SetRescatter(rescatter);
  double kt2min(0.); 
  p_ladder->SetMinKT2(kt2min);
  if (FixFirstOutgoings()) {
    //msg_Out()<<"      ...............................................\n"
    //	     <<"      "<<METHOD<<" start ladder with:\n"<<(*p_ladder)
    //	     <<"      ...............................................\n";
    return m_IS.Weight();
  }
  return -1.;
}

bool Ladder_Generator::FixFirstOutgoings() {
  Vec4D inmom1(p_ladder->GetIn1()->m_mom);
  Vec4D inmom2(p_ladder->GetIn2()->m_mom);
  Vec4D outmom1, outmom2,qt;
  Flavour outflav1,outflav2;
  bool keep(true);
  if (!Fix2To2Outgoing(inmom1,inmom2,outmom1,outmom2,keep)) return false;
  outflav1 = p_ladder->GetIn1()->m_flav;
  outflav2 = p_ladder->GetIn2()->m_flav;
  qt       = inmom1-outmom1;

  Ladder_Particle part1(outflav1,outmom1);
  Ladder_Particle part2(outflav2,outmom2);

  p_ladder->AddParticle(outmom1.Y(),part1);
  p_ladder->AddParticle(outmom2.Y(),part2);
  T_Prop prop(colour_type::octet,qt,m_FS.Q02((outmom1.Y()+outmom2.Y())/2.));
  p_ladder->GetProps()->push_back(prop);
  p_ladder->SetMaxKT2(Max(outmom1.PPerp2(),outmom2.PPerp2()));
  return true;
}

bool Ladder_Generator::
Fix2To2Outgoing(const ATOOLS::Vec4D & inmom1,const ATOOLS::Vec4D & inmom2,
		ATOOLS::Vec4D & outmom1,ATOOLS::Vec4D & outmom2,bool keep) {
  if (keep) {
    outmom1 = inmom1; outmom2 = inmom2;
  }
  else {
    ATOOLS::Vec4D cms(inmom1+inmom2);  
    double shat  = cms.Abs2(), QT2min(0.0), QT2max = ATOOLS::Min(4.,shat/4.);
    double meany = (inmom1.Y()+inmom2.Y())/2.;
    double q1T2  = p_eikonal->FF1()->SelectQT2(QT2max,QT2min);
    double q2T2  = p_eikonal->FF2()->SelectQT2(QT2max,QT2min);
    double QT2   = sqrt(q1T2*q2T2); 
    double QT    = sqrt(QT2);
    if (m_output) m_histograms[string("QT")]->Insert(QT);
    double dy    = acosh(sqrt(shat/(4.*QT2)));
    double phi1  = 2.*M_PI*ATOOLS::ran->Get();
    double cphi1 = cos(phi1), cphi2 = -cphi1; 
    double sphi1 = sin(phi1), sphi2 = -sphi1; 
    outmom1      = QT*ATOOLS::Vec4D(cosh(dy),cphi1,sphi1,sinh(dy));
    outmom2      = QT*ATOOLS::Vec4D(cosh(dy),cphi2,sphi2,sinh(-dy));
    if (cms.PPerp()!=0.) {
      ATOOLS::Poincare rotate(cms,ATOOLS::Vec4D(1.,0.,0.,1.));
      rotate.RotateBack(outmom1);
      rotate.RotateBack(outmom2);
    }
    ATOOLS::Poincare boost(cms);
    boost.BoostBack(outmom1);
    boost.BoostBack(outmom2);
    if (((inmom1.Y()<inmom2.Y()) && (outmom1.Y()>outmom2.Y())) || 
        ((inmom1.Y()>inmom2.Y()) && (outmom1.Y()<outmom2.Y()))) { 
      Vec4D help = outmom1;
      outmom1 = outmom2;
      outmom2 = help;
    }
  }
  return true;
}

double Ladder_Generator::Weight(const double & isweight) {
  if (!p_ladder->ExtractHardest()) {
    msg_Error()<<"Error in "<<METHOD<<": "<<endl
	       <<"   Could not extract hardest 2->2 scatter in ladder:\n"
	       <<(*p_ladder)<<endl;
    return 0.;
  }
  double weight(1.);
  if (p_ladder->Size()>2) {
    //double tmin(p_ladder->Mu2());
    double tmin(Smin()/4.);
    double that(dabs(p_ladder->That()));
    Flavour in1,in2,out1,out2;
    if (!p_ladder->ReconstructMEFlavours(in1,in2,out1,out2)) return 0.;
    weight *=pow(tmin/(tmin+that),1.);
    /*
      double aslog(log(m_FS.AlphaS(Max(that,smin/4.))/m_FS.AlphaS(smin/4.)));
      double colfac1(out1.IsGluon()?3.:4./3.);
      double gamma1(out1.IsGluon()?-11./6.+10./9.:-3./2.);
      double colfac2(out2.IsGluon()?3.:4./3.);
      double gamma2(out2.IsGluon()?-11./6.+10./9.:-3./2.);
      double term1(colfac1/(2.*M_PI)*(sqr(aslog)+gamma1*aslog));
      double term2(colfac2/(2.*M_PI)*(sqr(aslog)+gamma2*aslog));
      weight *= exp(-term1-term2);

      if (p_ladder->IsHardDiffractive()) {
      //double Yhat(p_ladder->Yhat());
      double slad((p_ladder->GetIn1()->m_mom+p_ladder->GetIn2()->m_mom).Abs2());
      //double shat(p_ladder->Shat()), uhat(dabs(p_ladder->Uhat()));
      //double pt2(2.*that*uhat*shat/(shat*shat+uhat*uhat+that*that));
      double pt2(4.*that*uhat*shat/(shat*shat+uhat*uhat+that*that));
      weight *= sqr(m_FS.AlphaS(pt2)/m_FS.AlphaSMax());
      }
      //    else if (p_ladder->Size()==3) 
      //	  weight *= m_FS.AlphaS(pt2)/m_FS.AlphaSMax();
    */
  }
  if (m_output) m_histograms[string("LadderWt")]->Insert(weight);
  return weight;
}


// Ladder * Ladder_Generator::MakeLadder(Blob * blob) {
//   p_ladder = new Ladder(Vec4D(0.,0.,0.,0.));
// 
//   Particle * part;
//   Ladder_Particle * lpart, * lpart1(NULL), *lpart2(NULL);
//   Vec4D pextra(Vec4D(0.,0.,0.,0.));
//   list<Particle *> lpextra;
//   for (int i=0;i<blob->NOutP();i++) {
//     part = blob->OutParticle(i);    
//     if (part->Status()==part_status::active && part->Info()=='F') {
//       msg_Tracking()<<part->Number()<<" ("<<part->Momentum().Y()<<") ";
//       lpart = new Ladder_Particle(part);
//       p_ladder->AddParticle(part->Momentum().Y(),(*lpart));
//     }
//     if (part->Status()==part_status::decayed && part->Info()=='H') {
//       if (part->Flav().Strong()) {
// 	msg_Error()<<"Error in "<<METHOD<<":\n"
// 		   <<"   Did not expect strongly interacting particle with "
// 		   <<"decay blob here!\n"<<(*part)<<"\n";
// 	if (part->DecayBlob()) msg_Error()<<(*part->DecayBlob())<<"\n";
// 	exit(1);
//       }
//       lpextra.push_back(part);
//       pextra += part->Momentum();
//     }
//   }
//   lpart = new Ladder_Particle(Flavour(kf_cluster),pextra);
//   p_ladder->AddParticle(pextra.Y(),(*lpart));
//   if (lpextra.size()>0) {
//     for (list<Particle *>::iterator piter=lpextra.begin();
// 	 piter!=lpextra.end();piter++) msg_Tracking()<<(*piter)->Number()<<" ";
//   }
//   for (int i=0;i<blob->NInP();i++) {
//     part = blob->InParticle(i);    
//     if (part->Info()=='I' && part->ProductionBlob()==NULL) {
//       lpart = new Ladder_Particle(part->Flav(),part->Momentum());
//       if (rpa->gen.PBeam(0)[3]*lpart->m_mom[3]>0.) lpart1 = lpart;
//                                                      else lpart2 = lpart;
//     }
//   }
//   p_ladder->SetInParticles(lpart1,lpart2);
//   int trials(0);
//   double weight(1.);
//   do {
//     /*
//       if (m_IS.ProvideInitialState(p_ladder)) {
//       m_pos    = m_IS.Pos();
//       m_b1     = m_IS.B1();
//       m_b2     = m_IS.B2();
//       weight   = m_IS.Weight();
//       
//       m_FS.Init(p_eikonal,m_b1,m_b2);  
//       }
//     */
//   } while (trials<1000 && weight<ran->Get());
//   return p_ladder;
// }



void Ladder_Generator::Analyse(const bool & isprimary) {
  if (!m_output || !p_ladder) return;
  double y1(p_ladder->GetEmissionsBegin()->second.m_mom.Y());
  double y2(p_ladder->GetEmissionsRBegin()->second.m_mom.Y());
  double Deltay(dabs(y2-y1));
  int Nemit(p_ladder->Size());
  if (!isprimary) {
    m_Nsec++;
    m_histograms[string("D_y2")]->Insert(Deltay);
    m_histograms[string("Nemit2")]->Insert(Nemit);
    m_histograms[string("Delta2")]->Insert(double(Nemit-1)/Deltay);
  }
  if (isprimary) {
    m_Nprim++;
    m_histograms[string("D_y1")]->Insert(Deltay);
    m_histograms[string("Nemit1")]->Insert(Nemit);
    m_histograms[string("Delta1")]->Insert(double(Nemit-1)/Deltay);
    if (p_ladder->IsHardDiffractive()) {
      m_N1h++;
      m_N1++;
      m_histograms[string("N_hdiff")]->Insert(Nemit);
    }
    else if (p_ladder->IsDiffractive()) {
      m_N1++;
      m_histograms[string("N_hdiff")]->Insert(Nemit);
    }
    else {
      m_N8++;
      m_histograms[string("N_hdiff")]->Insert(Nemit);
    }
  }
  Histogram * histokt(m_histograms[isprimary?string("KT1"):string("KT2")]);
  Histogram * histoktmid(m_histograms[isprimary?string("KT1mid"):
				      string("KT2mid")]);
  Histogram * histoy(m_histograms[string("y1")]);
  double kt,y;
  Vec4D momtot(0.,0.,0.,0.);
  for (LadderMap::iterator piter=p_ladder->GetEmissionsBegin();
       piter!=p_ladder->GetEmissionsEnd();piter++) {
    kt = piter->second.m_mom.PPerp();
    y  = piter->second.m_mom.Y();
    histokt->Insert(kt);
    histoy->Insert(y);
    if (dabs(y)<4.) histoktmid->Insert(kt);
    momtot+=piter->second.m_mom;
  }
  if (p_ladder->GetEmissions()->size()==2 
      && p_ladder->GetPropsBegin()->m_col==colour_type::singlet) {
    if (isprimary) {
      m_Ndd_p++;
      m_histograms[string("DDmass_p")]->Insert(momtot.Abs());  
    }
    else {
      m_Ndd_s++;
      m_histograms[string("DDmass_s")]->Insert(momtot.Abs());  
    }
  }
  else {
    if (p_ladder->GetPropsBegin()->m_col==colour_type::singlet &&
        p_ladder->GetPropsRBegin()->m_col==colour_type::singlet) {
      momtot-=p_ladder->GetEmissionsBegin()->second.m_mom;
      momtot-=p_ladder->GetEmissionsRBegin()->second.m_mom;
      if (isprimary) {
        m_Ncep_p++;
        m_histograms[string("CEPmass_p")]->Insert(momtot.Abs());  
      }
      else {
        m_Ncep_s++;
        m_histograms[string("CEPmass_s")]->Insert(momtot.Abs());  
      }
    }
    else if (p_ladder->GetPropsBegin()->m_col==colour_type::singlet) {
      momtot-=p_ladder->GetEmissionsBegin()->second.m_mom;
      if (isprimary) {
        m_Nsd_p++;
        m_histograms[string("SDmass_p")]->Insert(momtot.Abs());  
      }
      else {
        m_Nsd_s++;
        m_histograms[string("SDmass_s")]->Insert(momtot.Abs());  
      }
    }
    else if (p_ladder->GetPropsRBegin()->m_col==colour_type::singlet) {
      momtot-=p_ladder->GetEmissionsRBegin()->second.m_mom;
      if (isprimary) {
        m_Nsd_p++;
        m_histograms[string("SDmass_p")]->Insert(momtot.Abs());  
      }
      else {
        m_Nsd_s++;
        m_histograms[string("SDmass_s")]->Insert(momtot.Abs());  
      }
    }      
  }
}

