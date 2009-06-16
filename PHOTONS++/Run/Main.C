#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Particle.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Standard_Model.H"
#include "ATOOLS/Phys/Blob.H"
#include "PHASIC++/Channels/Rambo.H"
#include "PHOTONS++/Main/Photons.H"
#include "PHOTONS++/Main/Dipole_Type.H"

// #define NEW__RANDOM
// #undef USING__ROOT

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__ROOT
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TApplication.h"
#include "TPaveLabel.h"
#include "TF1.h"
#endif

#include "ATOOLS/Org/Data_Reader.H"

using namespace std;
using namespace PHOTONS;
using namespace ATOOLS;
using namespace MODEL;

int main(int argc, char *argv[])
{
  if(argc<2) {
    cout<<"Usage: ./Photons <OUTPUT=OUTPUTLEVEL>"<<endl;
    abort();
  }
  ATOOLS::exh->Init();
  set_terminate(ATOOLS::Terminate);
  set_unexpected(ATOOLS::Terminate);
  signal(SIGSEGV,ATOOLS::SignalHandler);
  signal(SIGINT,ATOOLS::SignalHandler);
  signal(SIGBUS,ATOOLS::SignalHandler);
  signal(SIGFPE,ATOOLS::SignalHandler);
  signal(SIGABRT,ATOOLS::SignalHandler);
  signal(SIGTERM,ATOOLS::SignalHandler);
  signal(SIGXCPU,ATOOLS::SignalHandler);
  try {
    //ATOOLS::ParticleInit("./");
    has to be replaced by model and fragmentation handler init, because that's
    where particle initialization is done now.
    for (int i=1; i<argc;++i) {
      string par = string(argv[i]);
      string key,value;
      int equal  = par.find("=");
      if (equal!=-1) {
      value = par.substr(equal+1);
      key   = par = par.substr(0,equal);
      }
      if (equal!=-1) {
	Read_Write_Base::AddCommandLine(key+" = "+value+"; ");
      }
    }
    rpa.Init("./","Run.dat",argc,argv);
    MODEL::s_model = new Standard_Model("./","Model.dat",true);
    rpa.gen.SetEcms(1.0);
    msg->SetModifiable(true);

    short unsigned int nin  = 1;
    short unsigned int nout = 3;                                                // <--
    Flavour myflav[10] ;                                                         // <--
    myflav[0] = ATOOLS::Flavour(kf_tau);
    myflav[1] = ATOOLS::Flavour(kf_e);
    myflav[2] = ATOOLS::Flavour(kf_nue).Bar();
    myflav[3] = ATOOLS::Flavour(kf_nutau);
//     myflav[4] = ATOOLS::Flavour(kf_e).Bar();
//     myflav[5] = ATOOLS::Flavour(kf_pi_plus);
//     myflav[6] = ATOOLS::Flavour(kf_pi_plus).Bar();
//     myflav[7] = ATOOLS::Flavour(kf_pi);
//     myflav[8] = ATOOLS::Flavour(kf_nutau).Bar();
//     myflav[9] = ATOOLS::Flavour(kf_nutau).Bar();

    Dipole_Type::code dtype = Dipole_Type::fi;                                  // <--

    Vec4D mom[10];
    double p = 1; 
    double M    = myflav[0].HadMass();
    double sum  = 0;
    for (int i=nin; i<nout+nin-1; i++)  sum = sum + myflav[i].HadMass();
    double kmax = (M/2) * ( M/sum - sum/M ); 
    int count   = 0;
    mom[0] = Vec4D( sqrt(M*M + p*p), 0., 0., p );

    Poincare boost(mom[0]);
    boost.Boost(mom[0]);

    PHASIC::Rambo myrambo(nin,nout,&myflav[0],true);
    myrambo.GeneratePoint(&mom[0]);

    for (unsigned short int i=0; i<nout+1; i++) {
      boost.BoostBack(mom[i]);
    }

    Particle * parts[10];
    // let particle decay from its CMS
    Blob* blob = new Blob();
    blob->SetId(1);
    blob->SetStatus(blob_status::needs_extraQED);
    for (short unsigned int i=  0; i<     nin; i++) {
      parts[i] = new Particle(i+1,myflav[i],mom[i],'I');
      parts[i]->SetFinalMass(-1,-1);
      blob->AddToInParticles(parts[i]);
    }
    for (short unsigned int i=nin; i<nin+nout; i++) {
      parts[i] = new Particle(i+1,myflav[i],mom[i],'F');
      parts[i]->SetFinalMass(-1,-1);
      blob->AddToOutParticles(parts[i]);
    }
    msg_Out()<<*blob<<endl;

//     Data_Read * reader = new Data_Read(YFS.dat);
    Photons * photons = new Photons(/*reader,*/true);
    
#ifdef USING__ROOT
    TApplication myapp(string("Photons").c_str(),&argc,&(*argv));
#endif

    unsigned int nevents, bins;
    cout<<"Enter the number of events: "<<endl;
    cin>>nevents;
#ifdef USING__ROOT
    cout<<"Enter the number of bins: "<<endl;
    cin>>bins;
#endif

#ifdef USING__ROOT
//     gROOT->Reset();
//     gROOT->SetStyle("Plain");
    std::string filename = "test_Photon_Analysis";                      // <--
    cout<<"Filling histogram for process "<<filename<<endl;
    TFile * rootfile = NULL;
    if (rootfile) rootfile->Close();
    TH1D * hist1, * hist2, * hist3, * hist4, * hist5, * hist6, * hist7;
    rootfile = new TFile(("Analysis_Multi_Final/"+filename+".root").c_str(),"RECREATE");
    int kcount = 0;
    if (dtype == Dipole_Type::ff) {
      hist1    = new TH1D("E_p","E in decay-CMS",bins,0,M);
      hist2    = new TH1D("E_Q","E in original dipole-CMS",bins,0,M);
      hist3    = new TH1D("E_P","E in dipole-CMS",bins,0,kmax);
      hist4    = new TH1D("theta_orig_dip","theta in original dipole CMS",bins,0,M_PI);
      hist5    = new TH1D("theta_dip","theta in dipole CMS",bins,0,M_PI);
      hist6    = new TH1D("theta_dec","theta in decay frame",bins,0,M_PI);
      hist7    = new TH1D("theta_lab","theta in lab frame",bins,0,M_PI);
    }
    else if (dtype == Dipole_Type::fi) {
      hist1    = new TH1D("E","E in dipole-CMS",bins,0,M);
      hist2    = new TH1D("theta_dip","theta in dipole-CMS",bins,0,M_PI);
      hist3    = new TH1D("theta_dec","theta in decay-CMS",bins,0,M_PI);
      hist4    = new TH1D("theta_lab","theta in lab frame",bins,0,M_PI);
      hist5  = new TH1D("EE","E in decay-CMS",bins,0,M);
    }
#endif

#ifdef NEW__RANDOM
    time_t zeit ;
    ran.SetSeed(time(&zeit)) ;
#endif

    for(unsigned int i=1; i<nevents+1; i++) {
      Blob blobcopy = Blob(blob);
      bool anythingadded = photons->AddRadiation(&blobcopy);
      if (photons->DoneSuccessfully() == false) count++;
        
      msg_Info()<<"i: "<<i<<endl;
      msg_Events()<<"-------------------------------------------"<<endl;
      msg_Events()<<(blobcopy)<<endl;
      msg_Events()<<"-------------------------------------------"<<endl;
      msg_Info()<<"-------------------------------------------"<<endl;
        
#ifdef USING__ROOT
      if (anythingadded == true) {

      if (dtype == Dipole_Type::ff) {
        Vec3D zaxis = Vec3D(0.,0.,1.);
        Vec4D kmom;
        Vec3D k;

        // dipole rest frame before radiation
        double K0 = 0;
        Vec4D odip = Vec4D(0.,0.,0.,0.);
        if (dtype == Dipole_Type::fi)   odip = odip + blob->InParticle(0)->Momentum();
        for (unsigned int j=0; j<nout; j++) {
          if (blob->OutParticle(j)->Flav().Charge() != 0)
            odip = odip + blob->OutParticle(j)->Momentum();
        }
        Vec4D op1 = blob->OutParticle(0)->Momentum();
        Poincare odipboost(odip);
        odipboost.Boost(op1);
        // p1 in z-direction
        Poincare odiprot(op1,Vec4D(0.,0.,0.,1.));
        for (unsigned int j=nout; j<blobcopy.GetOutParticles().size(); j++) {
          if (blobcopy.OutParticle(j)->Flav().IsPhoton() == true) {
            kmom = blobcopy.OutParticle(j)->Momentum();
            odipboost.Boost(kmom);
            odiprot.Rotate(kmom);
            k = Vec3D(kmom);
            double theta = acos((k*zaxis)/k.Abs());
            hist4->Fill(theta);
            K0 = K0 + kmom[0];
            kcount++;
          }
        }
        if (K0 != 0) hist2->Fill(K0);

        // dipole rest frame after radiation
        K0 = 0;
        Vec4D dip = Vec4D(0.,0.,0.,0.);
        if (dtype == Dipole_Type::fi)   dip = dip + blobcopy.InParticle(0)->Momentum();
        for (unsigned int j=0; j<nout; j++) {
          if (blobcopy.OutParticle(j)->Flav().Charge() != 0)
            dip = dip + blobcopy.OutParticle(j)->Momentum();
        }
        Vec4D p1 = blobcopy.OutParticle(0)->Momentum();
        Poincare dipboost(dip);
        dipboost.Boost(p1);
        // p1 in z-direction
        Poincare diprot(p1,Vec4D(0.,0.,0.,1.));
        for (unsigned int j=nout; j<blobcopy.GetOutParticles().size(); j++) {
          if (blobcopy.OutParticle(j)->Flav().IsPhoton() == true) {
            kmom = blobcopy.OutParticle(j)->Momentum();
            dipboost.Boost(kmom);
            diprot.Rotate(kmom);
            k = Vec3D(kmom);
            double theta = acos((k*zaxis)/k.Abs());
            hist5->Fill(theta);
            K0 = K0 + kmom[0];
          }
        }
        if (K0 != 0) hist3->Fill(K0);

        // decay rest frame
        K0 = 0;
        Vec4D dec = blob->InParticle(0)->Momentum();
        // decaying particle in z-direction
        Poincare decrot(dec,Vec4D(0.,0.,0.,1.));
        decrot.Rotate(dec);
        Poincare decboost(dec);
        for (unsigned int j=nout; j<blobcopy.GetOutParticles().size(); j++) {
          if (blobcopy.OutParticle(j)->Flav().IsPhoton() == true) {
            kmom = blobcopy.OutParticle(j)->Momentum();
            decrot.Rotate(kmom);
            decboost.Boost(kmom);
            k = Vec3D(kmom);
            double theta = acos((k*zaxis)/k.Abs());
            hist6->Fill(theta);
            K0 = K0 + kmom[0];
          }
        }
        if (K0 != 0) hist1->Fill(K0);

        // lab frame
        // z-axis as defined by the momentum of decaying particle
        for (unsigned int j=nout; j<blobcopy.GetOutParticles().size(); j++) {
          if (blobcopy.OutParticle(j)->Flav().IsPhoton() == true) {
            kmom = blobcopy.OutParticle(j)->Momentum();
            k = Vec3D(kmom);
            double theta = acos((k*zaxis)/k.Abs());
            hist7->Fill(theta);
          }
        }
      }



      else if (dtype == Dipole_Type::fi) {
        msg_Info()<<"Dipole_Type::fi"<<endl;
        // only for charged->charged + N neutral decays
        Particle_Vector dip;
        dip.push_back(blobcopy.InParticle(0));
        for (short unsigned int j=0; j<nout; j++) {
          if (blobcopy.OutParticle(j)->Flav().Charge() != 0)
            dip.push_back(blob->OutParticle(j));
        }
        Vec3D zaxis = Vec3D(0.,0.,1.);
        Vec4D kmom;
        // lab angle and total photon energy in lab frame
        Vec4D plab = blobcopy.InParticle(0)->Momentum();
        double K0 = 0;
        Poincare labrot(plab,Vec4D(0.,0.,0.,1.));
        for (short unsigned int j=nout; j<blobcopy.GetOutParticles().size(); j++) {
          if (blobcopy.OutParticle(j)->Flav().IsPhoton() == true) {
            kmom = blobcopy.OutParticle(j)->Momentum();
            labrot.Rotate(kmom);
            Vec3D k = Vec3D(kmom);
            double theta_lab = acos((k*zaxis)/k.Abs());
            hist4->Fill(theta_lab);
            K0 = kmom[0] +K0;
          }
        }
        // dec angle and total photon energy in decay particle'rest frame
        Vec4D pdec = blobcopy.InParticle(0)->Momentum();
        K0 = 0;
        Poincare decboost(pdec);
        decboost.Boost(pdec);
        Poincare decrot(pdec,Vec4D(0.,0.,0.,1.));
        for (short unsigned int j=nout; j<blobcopy.GetOutParticles().size(); j++) {
          if (blobcopy.OutParticle(j)->Flav().IsPhoton() == true) {
            kmom = blobcopy.OutParticle(j)->Momentum();
            decboost.Boost(kmom);
            decrot.Rotate(kmom);
            Vec3D k = Vec3D(kmom);
            double theta_dec = acos((k*zaxis)/k.Abs());
            hist3->Fill(theta_dec);
            K0 = kmom[0] +K0;
          }
        }
        if (K0 != 0) hist5->Fill(K0);
        // dip angle and total photon energy in dipole rest frame
        Vec4D pdip1 = dip.at(0)->Momentum();
        Vec4D pdip2 = dip.at(1)->Momentum();
        K0 = 0;
        Poincare dipboost(pdip1+pdip2);
        dipboost.Boost(pdip1);
        dipboost.Boost(pdip2);
        Poincare diprot(pdip1-pdip2,Vec4D(0.,0.,0.,1.));
        for (short unsigned int j=nout; j<blobcopy.GetOutParticles().size(); j++) {
          if (blobcopy.OutParticle(j)->Flav().IsPhoton() == true) {
            kmom = blobcopy.OutParticle(j)->Momentum();
            dipboost.Boost(kmom);
            diprot.Rotate(kmom);
            Vec3D k = Vec3D(kmom);
            double theta_dip = acos((k*zaxis)/k.Abs());
            hist2->Fill(theta_dip);
            K0 = kmom[0] +K0;
            kcount++;
          }
        }
        if (K0 != 0)  hist1->Fill(K0);
      }

      }
#endif

#ifdef USING__ROOT
      if(i%10000==0 && i%100000!=0) msg_Out()<<i<<" events passed."<<std::endl;
      if(i%100000==0) msg_Out()<<i<<" events passed....    "<<hist1->GetEntries()<<" photons generated..."<<std::endl;
#else
      if(i%10000==0) msg_Out()<<i<<" events passed."<<std::endl;
#endif
    }
    msg_Out()<<"Generated "<<nevents<<" events"<<endl;
    msg_Out()<<"Number of failed events: "<<count<<endl;
    delete photons;
    delete MODEL::s_model;
    delete blob;
#ifdef USING__ROOT
    if (dtype == Dipole_Type::unknown) {
      TCanvas * ch = new TCanvas("ch","generated distribution",0,0,1000,1000);
      TPad * pad0  = new TPad("pad0","title",0.01,0.96,0.99,0.99);
      pad0->Draw();
      TPad * pad1  = new TPad("pad1","photon energy distribution in decay frame",0.01,0.52,0.49,0.95);
      pad1->SetLogy(1);
      pad1->Draw();
      TPad * pad2  = new TPad("pad2","angular distribution in dipole CMS",0.01,0.08,0.49,0.51);
      pad2->Draw();
      TPad * pad3  = new TPad("pad3","angular distribution in decay frame",0.51,0.52,0.99,0.95);
      pad3->Draw();
      TPad * pad4  = new TPad("pad4","angular distribution in Lab",0.51,0.08,0.99,0.51);
      pad4->Draw();
      TPad * pad5  = new TPad("pad5","",0.01,0.01,0.99,0.07);
      pad5->Draw();
      gPad = pad0;
      char name[100];
      sprintf(name,filename.c_str());
      TPaveLabel * title = new TPaveLabel(0.1,0.1,0.9,0.9,name);
      title->SetBorderSize(0);
//       title->BorderMode(0);
      title->Draw();
      gPad = pad1;
      hist1->Draw();
      gPad = pad2;
      hist5->Draw();
      gPad = pad3;
      hist6->Draw();
      gPad = pad4;
      hist7->Draw();
      gPad = pad5;
      char string1[100];
      sprintf(string1,"number of events: %d",nevents);
      char string2[100];
      sprintf(string2,"number of photons generated: %d",kcount);
      char string3[100];
      sprintf(string3,"3-momentum in lab-frame of decaying particle: ( 0 , 0 , %f )",p);
      TPaveLabel * text1 = new TPaveLabel(0.1,0.51,0.4,0.99,string1);
      TPaveLabel * text2 = new TPaveLabel(0.6,0.51,0.9,0.99,string2);
      TPaveLabel * text3 = new TPaveLabel(0.1,0.01,0.9,0.49,string3);
      text1->Draw();
      text2->Draw();
      text3->Draw();

      TCanvas * ch2 = new TCanvas("ch2","energy distributions in different frames",0,0,1000,600);
      gPad = ch2;
      TPad * pad00  = new TPad("pad00","energy in different frames",0.01,0.91,0.99,0.99);
      pad00->Draw();
      TPad * pad6   = new TPad("pad6","energy in decay CMS",0.01,0.46,0.49,0.89);
      pad6->SetLogy(1);
      pad6->Draw();
      TPad * pad7   = new TPad("pad7","energy in dipole CMS",0.51,0.46,0.99,0.89);
      pad7->SetLogy(1);
      pad7->Draw();
      TPad * pad8   = new TPad("pad8","energy in original dipole CMS",0.26,0.01,0.74,0.44);
      pad8->SetLogy(1);
      pad8->Draw();
      gPad = pad6;
      hist1->Draw();
      gPad = pad00;
      title->Draw();
      gPad = pad6;
      hist1->Draw();
      gPad = pad7;
      hist3->Draw();
      gPad = pad8;
      hist2->Draw();
    }

    else if (dtype == Dipole_Type::unknown) {
      TCanvas * ch = new TCanvas("ch","generated distribution",0,0,1000,1000);
      TPad * pad0  = new TPad("pad0","title",0.01,0.96,0.99,0.99);
      pad0->Draw();
      TPad * pad1  = new TPad("pad1","photon energy distribution",0.01,0.52,0.49,0.95);
      pad1->SetLogy(1);
      pad1->Draw();
      TPad * pad2  = new TPad("pad2","angular distribution in primary decay products' CMS",0.01,0.08,0.49,0.51);
      pad2->Draw();
      TPad * pad3  = new TPad("pad3","angular distribution in decayed particle's CMS",0.51,0.52,0.99,0.95);
      pad3->Draw();
      TPad * pad4  = new TPad("pad4","angular distribution in Lab",0.51,0.08,0.99,0.51);
      pad4->Draw();
      TPad * pad5  = new TPad("pad5","",0.01,0.01,0.99,0.07);
      pad5->Draw();
      gPad = pad0;
      char name[100];
      sprintf(name,filename.c_str());
      TPaveLabel * title = new TPaveLabel(0.1,0.1,0.9,0.9,name);
      title->SetBorderSize(0);
//       title->BorderMode(0);
      title->Draw();
      gPad = pad1;
      hist5->Draw();
      gPad = pad2;
      hist2->Draw();
      gPad = pad3;
      hist3->Draw();
      gPad = pad4;
      hist4->Draw();
      gPad = pad5;
      char string1[100];
      sprintf(string1,"number of events: %d",nevents);
      char string2[100];
      sprintf(string2,"number of photons generated: %d",kcount);
      char string3[100];
      sprintf(string3,"3-momentum in lab-frame of decaying particle: ( 0 , 0 , %f )",p);
      TPaveLabel * text1 = new TPaveLabel(0.1,0.51,0.4,0.99,string1);
      TPaveLabel * text2 = new TPaveLabel(0.6,0.51,0.9,0.99,string2);
      TPaveLabel * text3 = new TPaveLabel(0.1,0.01,0.9,0.49,string3);
      text1->Draw();
      text2->Draw();
      text3->Draw();

      TCanvas * ch2 = new TCanvas("ch2","photon energy distribution in different frames",0,0,1000,500);
      gPad = ch2;
      TPad * pad6 = new TPad("pad6","energy distribution in dipole-CMS",0.01,0.01,0.49,0.89);
      pad6->SetLogy(1);
      pad6->Draw();
      TPad * pad7 = new TPad("pad7","energy distribution in decaying particle's CMS",0.51,0.01,0.99,0.89);
      pad7->SetLogy(1);
      pad7->Draw();
      TPad * pad8 = new TPad("pad8","photon energy distribution in different frames",0.01,0.91,0.99,0.99);
      pad8->Draw();
      gPad = pad6;
      hist1->Draw();
      gPad = pad7;
      hist5->Draw();
      gPad = pad8;
      TPaveLabel * title1 = new TPaveLabel(0.1,0.1,0.9,0.9,name);
      title1->SetBorderSize(0);
      title1->Draw();
    }

//     myapp.Run(kTRUE);
    if (dtype == Dipole_Type::ff) {
      hist1->Write();
      hist2->Write();
      hist3->Write();
      hist4->Write();
      hist5->Write();
      hist6->Write();
      hist7->Write();
    }
    else if (dtype == Dipole_Type::fi) {
      hist1->Write();
      hist2->Write();
      hist3->Write();
      hist4->Write();
      hist5->Write();
    }
    rootfile->Close();
#endif
  }
  catch (Exception Sherpaexception) {
    Sherpaexception.UpdateLogFile();
    msg_Error()<<Sherpaexception<<endl;
    terminate();
  }
  catch (std::exception stdexception) {
    cout<<"Sherpa: throws exception "
	     <<stdexception.what()<<" ..."<<endl;
    terminate();
  }
}
