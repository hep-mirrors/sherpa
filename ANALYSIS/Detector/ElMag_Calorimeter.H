#ifndef ElMag_Calorimeter_h
#define ElMag_Calorimeter_h

#include "Detector_Element.H"
#include "Histogram.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class ElMag_Calorimeter : public Detector_Element {
  private:
    double Plaquette(Cell * cell,double & eta,double & phi,
		     const int strip1,const int strip2,const int cell1,const int cell2);
    void   FillPlaquette(Cell * cell,std::vector<Cell *> & cells,
			 const int etastrip,const int cellno,const int dim);
  public:
    ElMag_Calorimeter(Primitive_Analysis *);
    virtual ~ElMag_Calorimeter();

    void Evaluate(const ATOOLS::Blob_List &,double,int) {}
    Analysis_Object * GetCopy() const;

    void Reset();
    bool Fill(const double E,const double eta,const double phi,ATOOLS::Particle * part);

    bool BuildCluster(Cell * cell,std::vector<Cell *> & cluster,
		      const int dim,double & E,double & eta,double & phi);
  };
}

#endif
