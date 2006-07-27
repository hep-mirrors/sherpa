#ifndef Root_Histogram_H
#define Root_Histogram_H

#include "Histogram_Base.H"
#include "TH1D.h"
#include "Message.H"

namespace ATOOLS {

  class Root_Histogram : public Histogram_Base {
    bool                     m_finished, m_written, m_legend, m_logy;
    std::string              m_outfile;
    std::vector<std::string> m_data;

    int                      m_type, m_depth, m_fills;
    TH1D                   * p_roothisto;
    std::vector<TH1D *>      m_rootcomps;
    bool                     m_mustdeletecomps;
    std::string              m_title;
    double                   m_lxmin,m_lymin,m_lxmax,m_lymax;
  public :
    Root_Histogram();
    Root_Histogram(const std::string,const std::string,
		   std::vector<std::string> *);
    Root_Histogram(Root_Histogram *);
    Root_Histogram(const Histogram *)                  { std::cout<<METHOD<<std::endl; abort(); }
    ~Root_Histogram();

    void   Reset()                                     { p_roothisto->Reset(); }
    void   Insert(double x);                            
    void   Insert(double x, double weight, int ntimes); 
    void   Insert(std::string bin, double weight, int ntimes);
    double Mean()     const                            { return p_roothisto->GetMean(); }
    double Integral() const                            { return p_roothisto->Integral(); }
    double Ymax()     const                            { return p_roothisto->GetMaximum(); }
    double Ymin()     const                            { return p_roothisto->GetMinimum(); }

    void   InsertRange(double start, double end, double value) { std::cout<<METHOD<<std::endl; abort(); }
    double LogCoeff() const { abort(); }

    void   Scale(double factor)                        {  p_roothisto->Scale(factor); }
    void   Output();
    void   Output(const std::string);
    void   Finalize();
    void   Restore()                                   { std::cout<<METHOD<<std::endl; abort(); } 
    void   SetBin(int i, double x)                     { p_roothisto->SetBinContent(i,x); }
    void   Extrapolate(double x, double * ys,int mode) { std::cout<<METHOD<<std::endl; abort(); }

    // basic access methods
    int    Type()           const { return m_type; }
    int    Depth()          const { return m_depth; }
    int    Nbins()          const { return p_roothisto->GetNbinsX(); }
    double Xmin()           const { return p_roothisto->GetXaxis()->GetXmin(); }
    double Xmax()           const { return p_roothisto->GetXaxis()->GetXmax(); }
    double Value(int i)     const { return p_roothisto->GetBinContent(i); }
    double Fills()          const { return m_fills; }
    double BinSize(int i=1) const { return p_roothisto->GetBinWidth(i); }

    inline void SetTitle(const std::string title) { m_title=title; }
    inline std::string Title() const              { return m_title; }
    inline void SetLogY()                         { m_logy=true; }
    inline void SetLegend(const bool legend,
			  const double xmin,const double ymin,
			  const double xmax,const double ymax) 
    {
      m_legend=legend;m_lxmin=xmin;m_lymin=ymin;m_lxmax=xmax;m_lymax=ymax;
    }

    Root_Histogram & operator+=(const Root_Histogram & histo);
    Root_Histogram & operator+=(const Histogram & histo) { std::cout<<METHOD<<std::endl; abort(); }
  };
}

#endif 
