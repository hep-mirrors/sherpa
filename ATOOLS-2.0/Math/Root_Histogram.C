#include "Root_Histogram.H"
#include "Message.H"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace ATOOLS;

Root_Histogram::Root_Histogram() :
  Histogram_Base(), 
  m_finished(false), m_written(false), 
  m_type(100), p_roothisto(NULL), m_mustdeletecomps(false),
  m_legend("Off"), m_lxmin(-1.), m_lymin(-1.), m_lxmax(-1.), m_lymax(-1.)
{
  // std::cout<<"In Root_Histogram."<<std::endl;
  abort();
}

Root_Histogram::Root_Histogram(const std::string filename,const std::string outfile,
			       std::vector<std::string> & data) :
  m_finished(false), m_written(false), 
  m_outfile(outfile), m_type(100), p_roothisto(NULL), m_mustdeletecomps(true),
  m_legend("Off"), m_lxmin(-1.), m_lymin(-1.), m_lxmax(-1.), m_lymax(-1.)
{
  for (std::vector<std::string>::iterator dit=data.begin(); 
       dit!=data.end();dit++) {
    m_data.push_back((*dit));
  }
  bool found(false);
  ifstream test((filename).c_str());
  if (test.good()) found=true;
  test.close();
  if (found) {
    // std::cout<<"Look for the file."<<std::endl;
    TFile file((filename).c_str());
    TH1D * th1d(NULL);
    for (std::vector<std::string>::iterator dit=m_data.begin(); 
	 dit!=m_data.end();dit++) {
      th1d = (TH1D*)(file.Get((*dit).c_str()));
      th1d->SetDirectory(0);
      th1d->SetTitle((*dit).c_str());
      m_rootcomps.push_back(th1d);
    }
    file.Close();
    p_roothisto = (TH1D*)((*m_rootcomps.begin())->Clone("SHERPA"));
    p_roothisto->Reset();
    return;
  }
  else {
    msg.Error()<<"ERROR in Root_Histogram("<<filename<<") : "<<std::endl
	       <<"   File not found, abort."<<std::endl;
    abort();
  }
}

Root_Histogram::Root_Histogram(Root_Histogram * rhisto) :
  m_finished(false), m_written(false), 
  m_outfile(rhisto->m_outfile), m_type(100), p_roothisto(NULL), m_mustdeletecomps(false),
  m_legend(rhisto->m_legend), 
  m_lxmin(rhisto->m_lxmin), m_lymin(rhisto->m_lymin),
  m_lxmax(rhisto->m_lxmax), m_lymax(rhisto->m_lymax)

{
  for (std::vector<std::string>::iterator dit=rhisto->m_data.begin(); 
       dit!=rhisto->m_data.end();dit++) m_data.push_back((*dit));
  for (std::vector<TH1D*>::iterator th1dit=rhisto->m_rootcomps.begin(); 
       th1dit!=rhisto->m_rootcomps.end();th1dit++) m_rootcomps.push_back((*th1dit));
  p_roothisto = (TH1D*)((*m_rootcomps.begin())->Clone("SHERPA"));
  p_roothisto->Reset();
}

Root_Histogram::~Root_Histogram()
{
  if (p_roothisto) { delete p_roothisto; p_roothisto = NULL; }
  if (!m_rootcomps.empty()) {
    if (m_mustdeletecomps) {
      for (std::vector<TH1D*>::iterator th1dit=m_rootcomps.begin(); 
	   th1dit!=m_rootcomps.end();th1dit++) {
	if ((*th1dit)) { delete (*th1dit); (*th1dit)=NULL; }
      }
    }
    m_rootcomps.clear();
  }
}



Root_Histogram & Root_Histogram::operator+=(const Root_Histogram & histo)
{
  p_roothisto->Add(histo.p_roothisto,1.);
  return (*this);
}

void Root_Histogram::Output() {
  if (!msg.LevelIsDebugging()) return;
}


void Root_Histogram::Output(const std::string name)
{
  if (!m_written) {
    TFile file1((name+".root").c_str(),"NEW","test");
    p_roothisto->Write();
    file1.Close();

    TFile file((name+".total.root").c_str(),"NEW","test");
    if (m_finished) {
      TCanvas canvas((name).c_str());
      p_roothisto->SetTitle(m_title.c_str());
      p_roothisto->SetLineColor(kRed);
      TLegend legend(m_lxmin,m_lymin,m_lxmax,m_lymax);
      p_roothisto->Draw("E1");
      legend.AddEntry(p_roothisto,"Sherpa");
      for (std::vector<TH1D*>::iterator th1dit=m_rootcomps.begin(); 
	   th1dit!=m_rootcomps.end();th1dit++) {
	(*th1dit)->Draw("same E1");
        legend.AddEntry((*th1dit),(*th1dit)->GetTitle());
	(*th1dit)->Write();
      }
      p_roothisto->Write();
      if (m_legend=="On") {legend.Draw();legend.Write();}
      canvas.Write();
      canvas.Print((name+".eps").c_str(),"eps");
      canvas.Close();
    }
    file.Close();
    m_written=true;
  }
}

void Root_Histogram::Insert(double x)
{
  p_roothisto->Fill(x); 
}

void Root_Histogram::Insert(double x, double weight, int ntimes) 
{
  p_roothisto->Fill(x,weight); 
}

void Root_Histogram::Finalize()
{ 
  if (!m_finished) {
    m_finished=true;
    p_roothisto->Sumw2();
    double sum(1./p_roothisto->Integral(""));
    p_roothisto->Scale(sum);    
  }
}
