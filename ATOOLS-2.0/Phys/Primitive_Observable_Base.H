#ifndef Primitive_Observable_Base_H
#define Primitive_Observable_Base_H

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "Selector.H"
#include "Vector.H"
#include "Flavour.H"
#include "Histogram.H"
#include "Message.H"

#include "Parton_List.H"
#include "Blob_List.H"

namespace APHYTOOLS {

  class Primitive_Observable_Base {  
  protected:
    int         type,nbins;
    double      xmin,xmax;

    std::string name;
    AMATOOLS::Histogram      * histo;
    APHYTOOLS::Selector_Base * sel;
    int          nout;
    Flavour    * flavs;
    AMATOOLS::Vec4D * moms;
  public:
    Primitive_Observable_Base() {
      xmin  = xmax  = 0.;
      type  = nbins = 0;
      histo = 0;
      sel   = 0;
      name  = std::string("noname");
      nout  = 0;
      flavs = 0;
      moms  = 0;
    };

    Primitive_Observable_Base(Primitive_Observable_Base * old) {
      type  = old->type;
      sel   = old->sel;
      xmin  = old->xmin;
      xmax  = old->xmax;
      nbins = old->nbins;
      name  = old->name;
      histo = new AMATOOLS::Histogram(old->histo);
    }

    Primitive_Observable_Base(int _type,double _xmin,double _xmax,int _nbins,
			      APHYTOOLS::Selector_Base * _sel) :
      type(_type), nbins(_nbins), xmin(_xmin), xmax(_xmax), sel(_sel) 
    { 
      AORGTOOLS::msg.Error()<<"Primitive_Observable_Base::Primitive_Observable_Base."<<std::endl;
    };

    virtual ~Primitive_Observable_Base() {
      if (histo!=0) { delete histo; histo = 0; }
    }

    virtual void Evaluate(double ) {
      AORGTOOLS::msg.Error()<<"Error in Primitive_Observable_Base::Evaluate : 0 "<<name<<std::endl;
    }

    virtual void SetFlavInfo(int _nout,AMATOOLS::Vec4D * _moms,APHYTOOLS::Flavour * _flavs) {
      nout=_nout;
      moms=_moms;
      flavs=_flavs;
    }

    virtual void Evaluate(int,AMATOOLS::Vec4D *,APHYTOOLS::Flavour *,double) {
      AORGTOOLS::msg.Error()<<"Error in Primitive_Observable_Base::Evaluate : A "<<name<<std::endl;
    }

    virtual void Evaluate(const APHYTOOLS::Parton_List &,double) {
      AORGTOOLS::msg.Error()<<"Error in Primitive_Observable_Base::Evaluate : B "<<name<<std::endl;
    }

    virtual void Evaluate(const APHYTOOLS::Blob_List &,double) {
      AORGTOOLS::msg.Error()<<"Error in Primitive_Observable_Base::Evaluate : C "<<name<<std::endl;
    }

    virtual void EndEvaluation() {
      histo->Finalize();
      AORGTOOLS::msg.Events()<<"  "<<name<<" : "<<std::endl;
      histo->Output();
      AORGTOOLS::msg.Events()<<std::endl<<std::endl;
    }

    virtual void EndEvaluation(double scale) {
      histo->Finalize();
      histo->Scale(scale);
      AORGTOOLS::msg.Events()<<"  "<<name<<" : "<<std::endl;
      histo->Output();
      AORGTOOLS::msg.Events()<<std::endl<<std::endl;
    }

    virtual Primitive_Observable_Base * GetCopy() {
      AORGTOOLS::msg.Error()<<"Error in Primitive_Observable_Base::GetCopy."<<std::endl;
      return 0;
    }

    virtual void Output(std::string pname) {
      int  mode_dir = 448;
      mkdir((pname).c_str(),mode_dir); 
      histo->Output((pname+std::string("/")+name).c_str());
    }

    int                        Type()  { return type; }
    int                        Nbins() { return nbins; }
    double                     Xmin()  { return xmin; }
    double                     Xmax()  { return xmax; }
    std::string                Name()  { return name; }
    AMATOOLS::Histogram      * Histo() { return histo; }
    APHYTOOLS::Selector_Base * Sel()   { return sel; }
  };
}

#endif
