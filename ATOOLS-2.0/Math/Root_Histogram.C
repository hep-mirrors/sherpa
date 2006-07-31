#include "Root_Histogram.H"
#include "Message.H"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace ATOOLS;

Root_Histogram::Root_Histogram() :
  Histogram_Base(), 
  m_finished(false), m_written(false), m_legend(false), m_logy(false), 
  m_type(100), p_roothisto(NULL), m_mustdeletecomps(false),
  m_lxmin(-1.), m_lymin(-1.), m_lxmax(-1.), m_lymax(-1.)
{
  abort();
}

Root_Histogram::Root_Histogram(const std::string filename,const std::string outfile,
			       std::vector<std::string> * data) :
  m_finished(false), m_written(false), m_legend(false), m_logy(false), 
  m_outfile(outfile), m_type(100), p_roothisto(NULL), m_mustdeletecomps(true),
  m_lxmin(-1.), m_lymin(-1.), m_lxmax(-1.), m_lymax(-1.)
{
  for (std::vector<std::string>::iterator dit=data->begin(); 
       dit!=data->end();dit++) {
    m_data.push_back((*dit));
  }
  bool found(false);
  ifstream test((filename).c_str());
  if (test.good()) found=true;
  test.close();
  if (found) {
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
    p_roothisto->Reset("ICE");
    return;
  }
  else {
    msg.Error()<<"ERROR in Root_Histogram("<<filename<<") : "<<std::endl
	       <<"   File not found, abort."<<std::endl;
    abort();
  }
}

Root_Histogram::Root_Histogram(Root_Histogram * rhisto) :
  m_finished(false), m_written(false), m_legend(rhisto->m_legend), m_logy(rhisto->m_logy), 
  m_outfile(rhisto->m_outfile), m_type(100), p_roothisto(NULL), m_mustdeletecomps(false),
  m_lxmin(rhisto->m_lxmin), m_lymin(rhisto->m_lymin),
  m_lxmax(rhisto->m_lxmax), m_lymax(rhisto->m_lymax)

{
  for (std::vector<std::string>::iterator dit=rhisto->m_data.begin(); 
       dit!=rhisto->m_data.end();dit++) m_data.push_back((*dit));
  for (std::vector<TH1D*>::iterator th1dit=rhisto->m_rootcomps.begin(); 
       th1dit!=rhisto->m_rootcomps.end();th1dit++) m_rootcomps.push_back((*th1dit));
  p_roothisto = (TH1D*)((*m_rootcomps.begin())->Clone("SHERPA"));
  p_roothisto->Reset("ICE");
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


void Root_Histogram::Output(const std::string filename)
{
  if (!m_written) {
    TFile file1((filename+".root").c_str(),"NEW","test");
    p_roothisto->Write();
    file1.Close();

    TFile file((filename+".total.root").c_str(),"NEW","test");
    if (m_finished) {
      std::string canvasstring = filename.substr(filename.rfind("/")+1);
      TCanvas canvas( canvasstring.c_str() );
      p_roothisto->SetTitle(m_title.c_str());
      p_roothisto->SetLineColor(kRed);
      p_roothisto->SetStats(kFALSE);
      TLegend legend(m_lxmin,m_lymin,m_lxmax,m_lymax);
      if (m_logy) gPad->SetLogy();
      p_roothisto->Draw("E1");
      legend.AddEntry(p_roothisto,"Sherpa");
      for (std::vector<TH1D*>::iterator th1dit=m_rootcomps.begin();
           th1dit!=m_rootcomps.end();th1dit++) {
        double compintegral = (*th1dit)->Integral("width");
        if( (*th1dit)->GetEntries() > 0 && compintegral>0.0) {
          (*th1dit)->Scale(1./compintegral);
          (*th1dit)->Draw("same E1");
          legend.AddEntry((*th1dit),(*th1dit)->GetTitle());
          (*th1dit)->Write();
        }
      }
      p_roothisto->Write();
      if (m_legend) {legend.Draw(); legend.Write();}
      canvas.Print((filename+".eps").c_str(),"eps");
      canvas.Write();
      canvas.Close();
    }
    file.Close();
    m_written=true;
  }
}

void Root_Histogram::Insert(double x)
{
  double binwidth(p_roothisto->GetBinWidth(p_roothisto->FindBin(x)));
  if(binwidth>0.0 ) p_roothisto->Fill(x,1./binwidth);
  else              p_roothisto->Fill(x);
}

void Root_Histogram::Insert(double x, double weight, int ntimes) 
{
  double binwidth(p_roothisto->GetBinWidth(p_roothisto->FindBin(x)));
  p_roothisto->Fill(x,weight/binwidth); 
}

void Root_Histogram::Insert(std::string bin, double weight, int ntimes) 
{
  double binwidth = 1.;
  p_roothisto->Fill(bin.c_str(),weight/binwidth); 
}

void Root_Histogram::Finalize()
{
  if (!m_finished) {
    m_finished=true;
    p_roothisto->Sumw2();
    double integral = p_roothisto->Integral("width");
    if(integral>0.0) p_roothisto->Scale(1./integral);
  }
}
