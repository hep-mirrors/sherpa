#include "Root_Histogram.H"
#include "Message.H"
#include "TFile.h"
#include "TCanvas.h"

using namespace ATOOLS;

Root_Histogram::Root_Histogram() :
  Histogram_Base(), 
  m_finished(false), m_type(100), p_roothisto(NULL), p_rootcomp(NULL)
{
  // std::cout<<"In Root_Histogram."<<std::endl;
  abort();
}

Root_Histogram::Root_Histogram(const std::string filename,
			       const std::string dataname) :
  m_finished(false), m_type(100), p_roothisto(NULL), p_rootcomp(NULL)
{
  std::cout<<"In Root_Histogram: "<<filename<<"/"<<dataname<<std::endl;
  bool found(false);
  ifstream test((filename).c_str());
  if (test.good()) found=true;
  test.close();
  if (found) {
    // std::cout<<"Look for the file."<<std::endl;
    TFile file((filename).c_str());
    p_rootcomp  = (TH1D*)(file.Get(dataname.c_str()));
    p_rootcomp->SetDirectory(0);
    file.Close();
    p_roothisto = (TH1D*)(p_rootcomp->Clone("SHERPA"));
    p_roothisto->Reset();
    return;
  }
  else {
    msg.Error()<<"ERROR in Root_Histogram("<<filename<<") : "<<std::endl
	       <<"   File not found, abort."<<std::endl;
    abort();
  }
}

Root_Histogram::Root_Histogram(const Root_Histogram * rhisto) :
  m_finished(false), m_type(100), p_roothisto(NULL), p_rootcomp(NULL)
{
  p_rootcomp  = new TH1D(*rhisto->p_rootcomp);
  p_roothisto = new TH1D(*rhisto->p_roothisto);
}

Root_Histogram::~Root_Histogram()
{
  if (p_roothisto) { delete p_roothisto; p_roothisto = NULL; }
  if (p_rootcomp)  { delete p_rootcomp;  p_rootcomp  = NULL; }
  // std::cout<<"Deleted the root-histogram."<<std::endl;
}


Root_Histogram & Root_Histogram::operator+=(const Root_Histogram & histo)
{
  p_roothisto->Add(histo.p_roothisto,1.);
  return (*this);
}

void Root_Histogram::Output() {
  if (!msg.LevelIsDebugging()) return;
  p_roothisto->Print("all");
}


void Root_Histogram::Output(const std::string name) 
{
  // std::cout<<"Output root file : "<<name<<std::endl;
  TFile file(name.c_str(),"NEW","test");
  p_roothisto->Write();
  file.Close();
  if (m_finished) {
    // std::cout<<"roothisto->Nbinsx: "<<p_roothisto->GetNbinsX()<<std::endl
    //	     <<"comphisto->Nbinsx: "<<p_rootcomp->GetNbinsX()<<std::endl;
    TCanvas * canvas = new TCanvas();
    p_roothisto->SetLineColor(kRed);
    p_roothisto->Draw();
    p_rootcomp->Draw("same");
    canvas->Print("test.eps","eps");
    delete canvas;
  }
}

void Root_Histogram::Insert(double x)
{
  // std::cout<<METHOD<<": Insert "<<x<<std::endl;
  p_roothisto->Fill(x); 
}

void Root_Histogram::Insert(double x, double weight, int ntimes) 
{
  // std::cout<<METHOD<<": Insert "<<x<<"("<<weight<<")."<<std::endl;
  p_roothisto->Fill(x,weight); 
}

void Root_Histogram::Finalize()
{ 
  if (!m_finished) {
    m_finished=true;
    p_roothisto->Sumw2();
    double sum(1./p_roothisto->Integral(""));
    p_roothisto->Scale(sum);    
    // std::cout<<METHOD<<":"<<sum<<std::endl;
  }
}
