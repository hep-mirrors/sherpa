#include "Histogram.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>

using namespace AMATOOLS;
using namespace AORGTOOLS;

Histogram::Histogram(int _type,double _lower,double _upper,int _nbins) :
  type(_type), lower(_lower), upper(_upper), nbins(_nbins), bins(0), fills(0)
{
  logarithmic = int(type/10);
  depth       = type-logarithmic*10+1;

  logbase = 1;
  switch(logarithmic) {
    case 1: 
      logbase = log(10.);
      upper   = log(upper)/logbase; lower = log(lower)/logbase;
      break;
    case 2: 
      upper   = log(upper); lower = log(lower);
      break;
    default: break;
  }
  binsize     = (upper-lower)/(double(nbins));

  if (binsize<=0.) {
    msg.Error()<<"Error in Histogram : Tried to initialize a histogram with  binsize <= 0 !"<<std::endl;
    active = 0;
    return;
  }

  nbins += 2;
  active = 1;
  bins   = new double * [nbins];

  for (int i=0;i<nbins;i++) {
    bins[i]   = new double[depth];
    for (int j=0;j<depth;j++) bins[i][j] = 0.;
  }

  msg.Debugging()<<"Histogram initialized : "
			    <<nbins-2<<" bins in "<<lower<<" ... "<<upper<<std::endl;
};

Histogram::Histogram(Histogram * histo) {
  lower   = histo->lower;
  upper   = histo->upper;
  logbase = histo->logbase;
  logarithmic = histo->logarithmic;
  nbins   = histo->nbins;
  depth   = histo->depth;
  type    = histo->type;
  fills   = histo->fills;

  binsize = histo->binsize;
  active  = 1;

  bins    = new double*[nbins];
  for (int i=0;i<nbins;i++) {
    bins[i]  = new double[depth]; 
    for (int j=0;j<depth;j++) bins[i][j] = histo->bins[i][j];
  }
  msg.Debugging()<<"Histogram initialized(copy) : "
			    <<nbins-2<<" bins in "<<lower<<" ... "<<upper<<std::endl;
}


Histogram::Histogram(std::string _pID) {
  char pID[100];
  sprintf(pID,"%s",_pID.c_str());

  std::ifstream ifile(pID);

  int    _type, _nbins;
  double _lower, _upper;
  ifile>>_type>>_nbins>>_lower>>_upper;

  type = _type; nbins = _nbins; lower = _lower; upper = _upper;

  logarithmic = int(type/10);
  depth       = type-logarithmic*10+1;

  logbase = 1;
  switch(logarithmic) {
    case 1: 
      logbase = log(10.);
      break;
    default: break;
  }
  binsize     = (upper-lower)/(double(nbins-2));

  if (binsize<=0.) {
    msg.Error()<<"Error in Histogram : "
	       <<"Tried to initialize a histogram with binsize <= 0 !"
	       <<std::endl;
    active = 0;
    return;
  }

  //  nbins += 2;
  active = 1;
  bins   = new double * [nbins];

  for (int i=0;i<nbins;i++) bins[i]   = new double[depth];

  for (int j=0;j<depth;j++) ifile>>bins[0][j];
  for (int j=0;j<depth;j++) ifile>>bins[nbins-1][j];

  double value;
  for (int i=0;i<nbins-1;i++) {
    ifile>>value;
    for (int j=0;j<depth;j++) ifile>>bins[i+1][j];
  }
  ifile.close();
}


Histogram::~Histogram() {
    
  if (bins!=0) { 
    for (int i=0;i<nbins;i++) 
      delete [] bins[i];
    delete [] bins; bins = 0; 
  }
    
};


void Histogram::Finalize() {
  double total = 0;
  for (int i=0;i<nbins;i++) { 
    //    total += bins[i][0]; 
    total += bins[i][0]/=fills; 
  }
  msg.Debugging()<<"Finalize histogram : "<<total<<"/"<<fills<<"="<<total/fills<<std::endl;
}

void Histogram::Reset() {
  double total = 0;
  for (int i=0;i<nbins;i++) { 
    bins[i][0]=0;
  }
  fills=0;
}

void Histogram::Scale(double scale) {
  double total = 0;
  for (int i=0;i<nbins;i++) { 
    bins[i][0]*= scale;
    total += bins[i][0]; 
  }
  fills = int(double(fills)/scale);
  msg.Debugging()<<"Scale histogram : "<<scale<<std::endl;
}

void Histogram::Output() {
  msg.Tracking()<<"----------------------------------------"<<std::endl;
  msg.Tracking()<<"    "<<bins[0][0]<<std::endl;
  msg.Tracking()<<"----------------------------------------"<<std::endl;
  double result = 0.;
  for (int i=0;i<nbins-2;i++) {
    msg.Tracking()<<lower+i*binsize<<"  ";
    for (int j=0;j<depth;j++) msg.Tracking()<<bins[i+1][j]<<"  ";
    result += bins[i+1][0];
    msg.Tracking()<<std::endl;
  }
  msg.Tracking()<<lower+(nbins-2)*binsize<<" == "<< upper<<std::endl;
  msg.Tracking()<<"----------------------------------------"<<std::endl;
  msg.Tracking()<<"    "<<bins[nbins-1][0]<<std::endl;
  msg.Tracking()<<"----------------------------------------"<<std::endl;
  msg.Tracking()<<"Inside the range : "<<result<<std::endl;
}


void Histogram::Output(std::string name) {
  std::ofstream ofile;
  ofile.open(name.c_str());

  ofile<<type<<" "<<nbins<<" "<<lower<<" "<<upper<<" ";
  for (int j=0;j<depth;j++) ofile<<bins[0][j]<<"  ";
  for (int j=0;j<depth;j++) ofile<<bins[nbins-1][j]<<"  ";
  ofile<<std::endl;
  for (int i=0;i<nbins-1;i++) {
    ofile<<lower+i*binsize<<"  ";
    for (int j=0;j<depth;j++) ofile<<bins[i+1][j]<<"  ";
    ofile<<std::endl;
  }
  ofile.close();

  msg.Out()<<"written file "<<name<<std::endl;
}







void Histogram::Insert(double coordinate) {
  if (!active) {
    msg.Error()<<"Error in Histogram : Tried to access a "
			  <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }

  fills++;

  if (logarithmic>0) coordinate = log(coordinate)/logbase;
  if (coordinate<lower) { bins[0][0]       += double(1); return; }
  if (coordinate>upper) { bins[nbins-1][0] += double(1); return; }
  for (int i=1;i<nbins-1;i++) {
    if ( (coordinate >= lower + (i-1)*binsize) &&
	 (coordinate <  lower + i*binsize) ) {
      bins[i][0] += double(1); 
      return; 
    }
  }
}

void Histogram::Insert(double coordinate,double value) {
  if (!active) {
    msg.Error()<<"Error in Histogram : Tried to access a "
			  <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }
  if (logarithmic>0) coordinate = log(coordinate)/logbase;

  fills++;

  if (coordinate<lower) { 
    bins[0][0] += value;
    if (depth>1) {
      if (value>bins[0][1]) bins[0][1] = value;
    }
    return; 
  }

  if (coordinate>upper) { 
    bins[nbins-1][0] += value; 
    if (depth>1) {
      if (value>bins[nbins-1][1]) bins[nbins-1][1] = value;
    }
    return; 
  }

  double low,up;
  low = lower; up = lower+binsize;
  for (int i=1;i<nbins-1;i++) {
    if ( (coordinate >= low) && (coordinate < up) ) {
      bins[i][0] += value;
      if (depth>1) {
	if (value>bins[i][1]) bins[i][1] = value;
      }
      return; 
    }
    low = up;
    up += binsize;
  }
}

void Histogram::InsertRange(double start, double end, double value) {
  if (!active) {
    msg.Error()<<"Error in Histogram : Tried to access a "
			  <<"histogram with binsize <= 0 !"<<std::endl;
    return;
  }
  if (logarithmic>0) {
    if (start>0)
      start = log(start)/logbase;
    else
      start = -30;
    if (end>0)
      end = log(end)/logbase;
    else 
      end = -30;
  }
  fills++;

  // underrun
  if (start<lower) { 
    bins[0][0] += value;
    if (depth>1) {
      if (value>bins[0][1]) bins[0][1] = value;
    }
  }

  // overflow
  if (start>upper) { 
    bins[nbins-1][0] += value; 
    if (depth>1) {
      if (value>bins[nbins-1][1]) bins[nbins-1][1] = value;
    }
  }

  double low,up;
  int hit=0;
  low = lower; up = lower+binsize;
  for (int i=1;i<nbins-1;i++) {
    if ((start < up) && (end >= low) ) {
      double fac=1;
      if ((start<=low)&&(up<=end )) {
	bins[i][0] += value;
	//	std::cout<<" case a "<<std::endl;
      } 
      else if ((low<start)&&(up<=end)) {
	fac = (start-low)/binsize;
	bins[i][0] += value *fac;
	//	std::cout<<" case b "<<std::endl;
      }
      else if ((start<=low)&&(end < up)) {
	fac = (up-end)/binsize;
	bins[i][0] += value *fac;
	//	std::cout<<" case c "<<std::endl;
      }
      else if ((low<start)&&(end <up)) {
	fac = (end-start)/binsize;
	bins[i][0] += value *fac;
// 	std::cout<<" case d "<<std::endl;
// 	std::cout<<" fill ("<<start<<","<<end<<") in ["<<low<<","<<up<<"]  fac="<<fac<<std::endl;
      }
    

      hit=1;
      if (depth>1) {
	if (value>bins[i][1]) bins[i][1] = value;
      }
    }
    low = up;
    up += binsize;
  }

  //  if (end==0. && hit==0) std::cout<<" not stored "<<std::endl;
}



double * Histogram::Bin(int bin) { return bins[bin]; }

double * Histogram::Bin(double coordinate) { 
  if (!active) {
    msg.Error()<<"Error in Histogram : Tried to access a histogram wit binsize <= 0 ! Return 0.."<<std::endl;
    return NULL;
  }
  else {
    if (logarithmic>0) coordinate = log(coordinate)/logbase;

    if (coordinate<lower) return bins[0];
    if (coordinate>upper) return bins[nbins-1];
    for (int i=1;i<nbins+1;i++) {
      if ( (coordinate >= lower + (i-1)*binsize) &&
	   (coordinate <  lower + i*binsize) ) 
	return bins[i];
    }
  }
  return NULL;
}


void Histogram::Extrapolate(double coordinate,double * res,int mode) {
  if (!active) {
    msg.Error()<<"Error in Histogram : Tried to access a histogram with binsize <= 0 ! Return 0.."<<std::endl;
  }
  else {
    if (logarithmic>0) coordinate = log(coordinate)/logbase;
    std::cout<<" in Extrapolate("<<coordinate<<","<<mode<<")"<<std::endl;
    std::cout<<" "<<lower<<" < "<<coordinate<<" <   "<<upper<<std::endl;

    for (int i=1;i<nbins;i++) {
      if ( (coordinate >= lower + (i-1)*binsize) &&
	   (coordinate <  lower + i*binsize) ) {
	// coordinate falls in bin no. i
	std::cout<<" i="<<i<<std::endl;
	// start of bin is lower + (i-1)*binsize 
	// end   of bin is lower +  i*binsize    

	res[0] = bins[i-1][0] +
	  (bins[i][0]-bins[i-1][0])/binsize *
	  (coordinate - lower - (i-1) * binsize);
	double over = 0.; double under = 0.;
	switch (mode) {
	case 1  : 
	  under = (coordinate - lower - (i-1) * binsize)/binsize * bins[i][0];
	  over = (lower + (i) * binsize - coordinate)/binsize * bins[i][0];
	  for (int j=0;j<i;j++) {
	    under += bins[j][0];
	  }
	  for (int j=i;j<nbins-1;j++) {
	    over += bins[j+1][0];
	  }
	  res[0]  = over;

	  if (depth) {
	    res[1]=0.;
	    for (int j=i;j<nbins;j++) {
	      res[1] = AMATOOLS::Max(res[1],bins[j][1]);
	    }
	  }


	  break;
	case -1 :
	  for (int j=i-1;j>0;j--) {
	    over  += bins[j][0];
	    under += bins[j-1][0];
	  }
	  over    += bins[0][0];
	  res[0]  += (over+under)/2;
	  break;
	default : break;
	}
	/*
	if (depth>0) {
	  for (int j=1;j<=depth;j++) 
	    res[j] = AMATOOLS::Max(bins[i][j],bins[i-1][j]);
	}
	*/

      }
    }
  }
}
