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
    int             type,nbins;
    double          xmin,xmax;
    Histogram     * histo;
    Selector_Base * sel;
    int             nout;
    Flavour       * flavs;
    Vec4D         * moms;
    std::string     name;
    std::string     m_blobtype;
    bool            m_blobdisc;
  public:
    Primitive_Observable_Base();
    Primitive_Observable_Base(int,double,double,int,Selector_Base *);
    Primitive_Observable_Base(Primitive_Observable_Base *);
    virtual ~Primitive_Observable_Base();

    virtual Primitive_Observable_Base * GetCopy()          { return NULL; }

    virtual void Evaluate(double)                          { }
    virtual void Evaluate(int,Vec4D *,Flavour *,double=1.) { }
    virtual void Evaluate(const Particle_List &,double=1.) { }
    virtual void Evaluate(const Blob_List &,double=1.);
    virtual void EndEvaluation(double=1.);


    virtual void Output(std::string);


    virtual void            SetFlavInfo(int,Vec4D *,Flavour *);
    virtual void            SetBlobType(std::string=std::string(""));

    virtual int             Type();
    virtual int             Nbins();
    virtual double          Xmin();
    virtual double          Xmax();
    virtual std::string     Name();
    virtual Histogram     * Histo();
    virtual Selector_Base * Sel();
  };
}

#endif
