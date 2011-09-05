/*
This is an example setup for an analysis of B-mixing properties in SHERPA.
It has to include Main_FullDecay.C, and the executable has to be
started with the following specialities:

  - ./HadronFullDecayMixingAnalysis DECAYER=300553 EVENTS=100000 ANALYSIS=1
  - set the decaydata directory in Datafiles/Fragmentation.dat to ./Decaydata
  - in there, HadronDecays.dat should have three lines, one for the Upsilon,
    one for the first B decay, and one for the second B decay, e.g.
    300553  ->   Upsilon4S/   Decays.dat;
    511     ->   Btag/        Decays.dat;
    511     ->   Bsignal/     Decays.dat;
  - the Upsilon should decay to B0 B0bar
  - the first B decay which happens should be chosen as a tagging decay
    mode (setup the Btag/Decays.dat file appropriately)
  - the second B decay should have the choice between the signal decay
    mode (e.g. J/psi KS) and another (unimportant) decay mode, such that
    we can see a rate asymmetry in the signal decay mode (setup the
    Bsignal/Decays.dat file appropriately).
*/

#include "Main_FullDecay.C"
#include "ATOOLS/Org/Shell_Tools.H"
#include "SHERPA/Initialization/Initialization_Handler.H"

static Hadron_Decay_Channel* signal_hdc;
#ifdef USING__ROOT
static TFile* rootfile;
static TH1D* B_tag_candidates;
static TH1D* Bbar_tag_candidates;
static TH1D* unmixed_B_events;
static TH1D* mixed_B_events;
static TH1D* unmixed_Bbar_events;
static TH1D* mixed_Bbar_events;
#endif

void InitialiseAnalysis()
{
#ifdef USING__ROOT
  Hadron_Decay_Map* decaymap = p_sherpa->GetInitHandler()
                               ->GetHadronDecayHandler("Hadrons")
                               ->GetHadrons()->DecayMap();
  Hadron_Decay_Table* decaytable = (*decaymap)[Flavour(kf_B)][1];
  signal_hdc = (Hadron_Decay_Channel*) decaytable->GetDecayChannel(0);

  std::string adir = "";
  ATOOLS::MakeDir(adir,493);
  rootfile = new TFile(string(adir+"/CPasymmetry_"+
                              signal_hdc->FileName()+".root").c_str(), "RECREATE");
  B_tag_candidates =    makeTH1D("B_tag_candidates", "",
                                 48, -12.0, 12.0,
                                 "#Deltat(ps)","Events [ps^{-1}]");
  Bbar_tag_candidates = makeTH1D("Bbar_tag_candidates", "",
                                 48, -12.0, 12.0,
                                 "#Deltat(ps)","Events [ps^{-1}]");
  unmixed_Bbar_events =      makeTH1D("unmixed_Bbar_events", "",
                                 48, -12.0, 12.0,
                                 "#Deltat(ps)", "Events [ps^{-1}]");
  mixed_Bbar_events =        makeTH1D("mixed_Bbar_events", "",
                                 48, -12.0, 12.0,
                                 "#Deltat(ps)", "Events [ps^{-1}]");
  unmixed_B_events =      makeTH1D("unmixed_B_events", "",
                                 48, -12.0, 12.0,
                                 "#Deltat(ps)", "Events [ps^{-1}]");
  mixed_B_events =        makeTH1D("mixed_B_events", "",
                                 48, -12.0, 12.0,
                                 "#Deltat(ps)", "Events [ps^{-1}]");
#endif
}


int Btag(Blob* blob)
{
  // tag > 0: contains l+ => B
  // tag < 0: contains l- => Bbar
  int tag=0;
  Particle_Vector out = blob->GetOutParticles();
  for(size_t i=0; i<out.size(); i++) {
    if(out[i]->Flav().IsLepton()) tag+=int(out[i]->Flav().Charge());
  }
  return int(tag);
}


void AnalyseEvent(Blob_List* blobs)
{
#ifdef USING__ROOT
  if(blobs->Find(btp::Hadron_Decay).size()<2)
    THROW(fatal_error, "not two B decay blobs.");
  Blob_List hadrondecayblobs = blobs->Find(btp::Hadron_Decay);
  bool foundfirst(false);
  Blob* tagBlob(NULL);
  Blob* decayBlob(NULL);
  for(size_t i=0; i<hadrondecayblobs.size(); ++i) {
    if(hadrondecayblobs[i]->InParticle(0)->Flav().Kfcode()==kf_B) {
      if(foundfirst) {
        decayBlob = hadrondecayblobs[i];
        break;
      }
      else {
        tagBlob = hadrondecayblobs[i];
        foundfirst = true;
      }
    }
  }
  double time = decayBlob->InParticle(0)->Time()-tagBlob->InParticle(0)->Time();
  if((*decayBlob)["hdc"]->Get<Hadron_Decay_Channel*>()==signal_hdc) {
    if(Btag(tagBlob)<0) Bbar_tag_candidates->Fill(time*1e12);
    else                B_tag_candidates->Fill(time*1e12);
  }

  if (tagBlob->InParticle(0)->Flav().IsAnti()) {
    if(blobs->FindFirst(btp::Hadron_Mixing))
      mixed_B_events->Fill(time*1e12);
    else
      unmixed_B_events->Fill(time*1e12);
  }
  else {
    if(blobs->FindFirst(btp::Hadron_Mixing))
      mixed_Bbar_events->Fill(time*1e12);
    else
      unmixed_Bbar_events->Fill(time*1e12);
  }
#endif
}


double Pmix(double* x, double* par)
{
  // par[0] < 0  ==> anti-meson mixing probability
  double t = x[0];
  double dm = Flavour(kf_B).DeltaM()/rpa->hBar()/1e12;
  double dG = Flavour(kf_B).DeltaGamma()/rpa->hBar()/1e12;
  double qoverp2=Flavour(kf_B).QOverP2();
  if(par[0]<0.0) qoverp2 = 1.0/qoverp2;
  double pmix = qoverp2 * (cosh(dG*t/2.0)-cos(dm*t));
  double pnotmix = (cosh(dG*t/2.0)+cos(dm*t));
  return pmix/(pmix+pnotmix);
}

double af(double* x, double* par)
{
  double t = x[0];
  double dm = Flavour(kf_B).DeltaM()/rpa->hBar()/1e12;
  double dG = Flavour(kf_B).DeltaGamma()/rpa->hBar()/1e12;
  double l2 = par[0];
  double Im_lambda = par[1];
  double Re_lambda = par[2];
  return (2.0*Im_lambda/(1.0+l2)*sin(dm*t) - (1.0-l2)/(1.0+l2)*cos(dm*t))/
    (cosh(dG*t/2.0) - 2.0*Re_lambda/(1.0+l2)*sinh(dG*t/2.0));
}


void FinishAnalysis()
{
#ifdef USING__ROOT
  TCanvas* c1 = new TCanvas("candidates","B and anti-B tag candidates",800,400);
  c1->SetLeftMargin(0.0);
  c1->SetRightMargin(0.0);
  c1->SetTopMargin(0.0);
  c1->SetBottomMargin(0.0);
  c1->Divide(1,2,0.0,0.0);
  c1->cd(1);
  c1->GetPad(1)->SetLeftMargin(0.06);
  c1->GetPad(1)->SetRightMargin(0.01);
  c1->GetPad(1)->SetTopMargin(0.05);
  c1->GetPad(1)->SetBottomMargin(0.0);
  B_tag_candidates->SetLineColor(kBlue);
  B_tag_candidates->SetFillStyle(3004);
  B_tag_candidates->SetFillColor(kBlue);
  B_tag_candidates->GetXaxis()->SetTitle("");
  B_tag_candidates->GetYaxis()->SetTitleOffset(0.35);
  B_tag_candidates->GetYaxis()->SetTitleSize(0.08);
  B_tag_candidates->Sumw2();
  B_tag_candidates->Draw("hist");
  B_tag_candidates->Print();
  B_tag_candidates->Write();
  Bbar_tag_candidates->SetLineColor(kRed);
  Bbar_tag_candidates->SetFillStyle(3005);
  Bbar_tag_candidates->SetFillColor(kRed);
  Bbar_tag_candidates->GetXaxis()->SetTitle("");
  Bbar_tag_candidates->Sumw2();
  Bbar_tag_candidates->Draw("hist same");
  Bbar_tag_candidates->Print();
  Bbar_tag_candidates->Write();
  TLegend* candidates_legend = new TLegend(0.6,0.6,0.89,0.89);
  candidates_legend->SetFillColor(kWhite);
  candidates_legend->AddEntry(B_tag_candidates, "B^{0} tag candidates");
  candidates_legend->AddEntry(Bbar_tag_candidates, "#bar{B}^{0} tag candidates");
  candidates_legend->Draw();

  c1->cd(2);
  c1->GetPad(2)->SetLeftMargin(0.06);
  c1->GetPad(2)->SetRightMargin(0.01);
  c1->GetPad(2)->SetTopMargin(0.0);
  c1->GetPad(2)->SetBottomMargin(0.1);
  c1->GetPad(2)->SetFillStyle(0);
  c1->GetPad(2)->SetGridx();
  c1->GetPad(2)->SetGridy();
  TH1D* candidate_asymmetry = new TH1D(((*B_tag_candidates)-(*Bbar_tag_candidates))/
                                       ((*B_tag_candidates)+(*Bbar_tag_candidates)));
  candidate_asymmetry->Sumw2();
  candidate_asymmetry->SetNameTitle("candidate_asymmetry",";#Deltat(ps)");
  candidate_asymmetry->SetLineColor(kBlack);
  candidate_asymmetry->SetFillStyle(0);
  candidate_asymmetry->GetYaxis()->SetTitle("Candidate Asymmetry");
  candidate_asymmetry->GetYaxis()->SetTitleOffset(0.35);
  candidate_asymmetry->GetYaxis()->SetTitleSize(0.08);
  candidate_asymmetry->GetYaxis()->SetLabelSize(0.06);
  candidate_asymmetry->GetXaxis()->SetTitleOffset(0.3);
  candidate_asymmetry->GetXaxis()->SetTitleSize(0.08);
  candidate_asymmetry->GetXaxis()->SetLabelSize(0.06);
  candidate_asymmetry->Draw();
  candidate_asymmetry->Write();
  string delta_m_ps = ToString(Flavour(kf_B).DeltaM()/rpa->hBar()/1e12);
  string delta_G_ps = ToString(Flavour(kf_B).DeltaGamma()/rpa->hBar()/1e12);
  string C = ToString(signal_hdc->CPAsymmetryC());
  string S = ToString(signal_hdc->CPAsymmetryS());
  TF1* candidate_asymmetry_theory = new TF1("candidate_asymmetry_theory",
    af,
    candidate_asymmetry->GetXaxis()->GetXmin(),
    candidate_asymmetry->GetXaxis()->GetXmax(),3);
  Complex lambda = signal_hdc->CPAsymmetryLambda();
  candidate_asymmetry_theory->SetParameters(sqr(abs(lambda)),
                                            lambda.imag(),
                                            lambda.real());
  candidate_asymmetry_theory->SetLineColor(kBlue);
  candidate_asymmetry_theory->SetLineWidth(1);
  candidate_asymmetry_theory->Draw("same");
  TLegend* candidates_asymmetry_legend = new TLegend(0.75,0.7,0.98,0.95);
  candidates_asymmetry_legend->SetFillColor(kWhite);
  candidates_asymmetry_legend->AddEntry(candidate_asymmetry, "Sherpa","L");
  candidates_asymmetry_legend->AddEntry(candidate_asymmetry_theory,
                                        ("Theory C="+C+" S="+S).c_str(),"L");
  candidates_asymmetry_legend->Draw();

  c1->Write();


  
  TCanvas* c2 = new TCanvas("Bmixing","Mixed and Unmixed B Events",800,400);
  c2->SetLeftMargin(0.0);
  c2->SetRightMargin(0.0);
  c2->SetTopMargin(0.0);
  c2->SetBottomMargin(0.0);
  c2->Divide(1,2,0.0,0.0);
  c2->cd(1);
  c2->GetPad(1)->SetLeftMargin(0.06);
  c2->GetPad(1)->SetRightMargin(0.01);
  c2->GetPad(1)->SetTopMargin(0.05);
  c2->GetPad(1)->SetBottomMargin(0.0);
  unmixed_B_events->SetLineColor(kBlue);
  unmixed_B_events->SetFillColor(kBlue);
  unmixed_B_events->SetFillStyle(3004);
  unmixed_B_events->GetXaxis()->SetTitle("");
  unmixed_B_events->GetYaxis()->SetTitleOffset(0.35);
  unmixed_B_events->GetYaxis()->SetTitleSize(0.08);
  unmixed_B_events->Sumw2();
  unmixed_B_events->Draw("hist");
  unmixed_B_events->Print();
  unmixed_B_events->Write();
  mixed_B_events->SetLineColor(kRed);
  mixed_B_events->SetFillColor(kRed);
  mixed_B_events->SetFillStyle(3005);
  mixed_B_events->Sumw2();
  mixed_B_events->Draw("hist same");
  mixed_B_events->Print();
  mixed_B_events->Write();
  TLegend* bmixing_legend = new TLegend(0.6,0.6,0.89,0.89);
  bmixing_legend->SetFillColor(kWhite);
  bmixing_legend->AddEntry(unmixed_B_events, "Unmixed Events");
  bmixing_legend->AddEntry(mixed_B_events, "Mixed Events");
  bmixing_legend->Draw();

  c2->cd(2);
  c2->GetPad(2)->SetLeftMargin(0.06);
  c2->GetPad(2)->SetRightMargin(0.01);
  c2->GetPad(2)->SetTopMargin(0.0);
  c2->GetPad(2)->SetBottomMargin(0.1);
  gPad->SetFillStyle(0);
  gPad->SetGridx();
  gPad->SetGridy();
  TH1D* bmixing_fraction = new TH1D(((*mixed_B_events))/
                                    ((*unmixed_B_events)+(*mixed_B_events)));
  bmixing_fraction->SetNameTitle("bmixing_fraction",";#Deltat(ps)");
  bmixing_fraction->SetLineColor(kBlack);
  bmixing_fraction->SetFillStyle(0);
  bmixing_fraction->GetYaxis()->SetTitle("#frac{mixed}{all} events");
  bmixing_fraction->GetYaxis()->SetTitleOffset(0.35);
  bmixing_fraction->GetYaxis()->SetTitleSize(0.08);
  bmixing_fraction->GetYaxis()->SetLabelSize(0.06);
  bmixing_fraction->GetXaxis()->SetTitleOffset(0.3);
  bmixing_fraction->GetXaxis()->SetTitleSize(0.08);
  bmixing_fraction->GetXaxis()->SetLabelSize(0.06);
  bmixing_fraction->Draw();
  bmixing_fraction->Write();
  TF1* bmixing_fraction_theory = new TF1("bmixing_fraction_theory",
    Pmix,
    bmixing_fraction->GetXaxis()->GetXmin(),
    bmixing_fraction->GetXaxis()->GetXmax(), 1);
  bmixing_fraction_theory->SetParameters(1.0, 0.0);
  bmixing_fraction_theory->SetLineColor(kBlue);
  bmixing_fraction_theory->SetLineWidth(1);
  bmixing_fraction_theory->Draw("same");
  TLegend* bmixing_fraction_legend = new TLegend(0.7,0.7,0.89,0.89);
  bmixing_fraction_legend->SetFillColor(kWhite);
  bmixing_fraction_legend->AddEntry(bmixing_fraction, "Sherpa","L");
  bmixing_fraction_legend->AddEntry(bmixing_fraction_theory, "Theory","L");
  bmixing_fraction_legend->Draw();

  c2->Write();


  TCanvas* c3 = new TCanvas("Bbarmixing","Mixed and Unmixed Bbar Events",800,400);
  c3->SetLeftMargin(0.0);
  c3->SetRightMargin(0.0);
  c3->SetTopMargin(0.0);
  c3->SetBottomMargin(0.0);
  c3->Divide(1,2,0.0,0.0);
  c3->cd(1);
  c3->GetPad(1)->SetLeftMargin(0.06);
  c3->GetPad(1)->SetRightMargin(0.01);
  c3->GetPad(1)->SetTopMargin(0.05);
  c3->GetPad(1)->SetBottomMargin(0.0);
  unmixed_Bbar_events->SetLineColor(kBlue);
  unmixed_Bbar_events->SetFillColor(kBlue);
  unmixed_Bbar_events->SetFillStyle(3004);
  unmixed_Bbar_events->GetXaxis()->SetTitle("");
  unmixed_Bbar_events->GetYaxis()->SetTitleOffset(0.35);
  unmixed_Bbar_events->GetYaxis()->SetTitleSize(0.08);
  unmixed_Bbar_events->Sumw2();
  unmixed_Bbar_events->Draw("hist");
  unmixed_Bbar_events->Print();
  unmixed_Bbar_events->Write();
  mixed_Bbar_events->SetLineColor(kRed);
  mixed_Bbar_events->SetFillColor(kRed);
  mixed_Bbar_events->SetFillStyle(3005);
  mixed_Bbar_events->Sumw2();
  mixed_Bbar_events->Draw("hist same");
  mixed_Bbar_events->Print();
  mixed_Bbar_events->Write();
  TLegend* bbarmixing_legend = new TLegend(0.6,0.6,0.89,0.89);
  bbarmixing_legend->SetFillColor(kWhite);
  bbarmixing_legend->AddEntry(unmixed_Bbar_events, "Unmixed Events");
  bbarmixing_legend->AddEntry(mixed_Bbar_events, "Mixed Events");
  bbarmixing_legend->Draw();

  c3->cd(2);
  c3->GetPad(2)->SetLeftMargin(0.06);
  c3->GetPad(2)->SetRightMargin(0.01);
  c3->GetPad(2)->SetTopMargin(0.0);
  c3->GetPad(2)->SetBottomMargin(0.1);
  gPad->SetFillStyle(0);
  gPad->SetGridx();
  gPad->SetGridy();
  TH1D* bbarmixing_fraction = new TH1D((*mixed_Bbar_events)/
                                    ((*unmixed_Bbar_events)+(*mixed_Bbar_events)));
  bbarmixing_fraction->SetNameTitle("bbarmixing_fraction",";#Deltat(ps)");
  bbarmixing_fraction->SetLineColor(kBlack);
  bbarmixing_fraction->SetFillStyle(0);
  bbarmixing_fraction->GetYaxis()->SetTitle("#frac{mixed}{all} events");
  bbarmixing_fraction->GetYaxis()->SetTitleOffset(0.35);
  bbarmixing_fraction->GetYaxis()->SetTitleSize(0.08);
  bbarmixing_fraction->GetYaxis()->SetLabelSize(0.06);
  bbarmixing_fraction->GetXaxis()->SetTitleOffset(0.3);
  bbarmixing_fraction->GetXaxis()->SetTitleSize(0.08);
  bbarmixing_fraction->GetXaxis()->SetLabelSize(0.06);
  bbarmixing_fraction->Draw();
  bbarmixing_fraction->Write();
  TF1* bbarmixing_fraction_theory = new TF1("bbarmixing_fraction_theory",
    Pmix,
    bbarmixing_fraction->GetXaxis()->GetXmin(),
    bbarmixing_fraction->GetXaxis()->GetXmax(), 1);
  bbarmixing_fraction_theory->SetParameters(-1.0, 0.0);
  bbarmixing_fraction_theory->SetLineColor(kBlue);
  bbarmixing_fraction_theory->SetLineWidth(1);
  bbarmixing_fraction_theory->Draw("same");
  TLegend* bbarmixing_fraction_legend = new TLegend(0.7,0.7,0.89,0.89);
  bbarmixing_fraction_legend->SetFillColor(kWhite);
  bbarmixing_fraction_legend->AddEntry(bbarmixing_fraction, "Sherpa","L");
  bbarmixing_fraction_legend->AddEntry(bbarmixing_fraction_theory, "Theory","L");
  bbarmixing_fraction_legend->Draw();

  c3->Write();
#endif
}
