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

#include "Particle_List.H"
#include "Blob_List.H"

namespace ATOOLS {

  class Primitive_Observable_Base {  
  protected:
    int         type,nbins;
    double      xmin,xmax;

    std::string name;
    Histogram      * histo;
    Selector_Base * sel;
    int          nout;
    Flavour    * flavs;
    Vec4D * moms;
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
      histo = new Histogram(old->histo);
    }

    Primitive_Observable_Base(int _type,double _xmin,double _xmax,int _nbins,
			      Selector_Base * _sel) :
      type(_type), nbins(_nbins), xmin(_xmin), xmax(_xmax), sel(_sel) 
    { 
      msg.Error()<<"Primitive_Observable_Base::Primitive_Observable_Base."<<std::endl;
    };

    virtual ~Primitive_Observable_Base() {
      if (histo!=0) { delete histo; histo = 0; }
    }

    virtual void Evaluate(double ) {
      msg.Error()<<"Error in Primitive_Observable_Base::Evaluate : 0 "<<name<<std::endl;
    }

    virtual void SetFlavInfo(int _nout,Vec4D * _moms,Flavour * _flavs) {
      nout=_nout;
      moms=_moms;
      flavs=_flavs;
    }

    virtual void Evaluate(int,Vec4D *,Flavour *,double) {
      msg.Error()<<"Error in Primitive_Observable_Base::Evaluate : A "<<name<<std::endl;
    }

    virtual void Evaluate(const Particle_List &,double) {
      msg.Error()<<"Error in Primitive_Observable_Base::Evaluate : B "<<name<<std::endl;
    }

    virtual void Evaluate(const Blob_List &,double) {
      msg.Error()<<"Error in Primitive_Observable_Base::Evaluate : C "<<name<<std::endl;
    }

    virtual void EndEvaluation() {
      histo->Finalize();
      msg.Events()<<"  "<<name<<" : "<<std::endl;
      histo->Output();
      msg.Events()<<std::endl<<std::endl;
    }

    virtual void EndEvaluation(double scale) {
      histo->Finalize();
      histo->Scale(scale);
      msg.Events()<<"  "<<name<<" : "<<std::endl;
      histo->Output();
      msg.Events()<<std::endl<<std::endl;
    }

    virtual Primitive_Observable_Base * GetCopy() {
      msg.Error()<<"Error in Primitive_Observable_Base::GetCopy."<<std::endl;
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
    Histogram      * Histo() { return histo; }
    Selector_Base * Sel()   { return sel; }
  };
}

#endif
