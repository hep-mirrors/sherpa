#include "Cut_Data.H"
#include "Run_Parameter.H"

using namespace APHYTOOLS;

Cut_Data::Cut_Data() {
  energymin = 0;
  energymax = 0;
  cosmin    = 0;
  cosmax    = 0;
  scut      = 0;
  ncut      = 0;
}

Cut_Data::~Cut_Data() {
  for (short int i=0;i<ncut;i++) {
    delete[] cosmin[i];
    delete[] cosmax[i];
    delete[] scut[i];
  }
  delete[] cosmin;
  delete[] cosmax;
  delete[] scut;
  delete[] energymin;
  delete[] energymax;
}

void Cut_Data::Init(int _ncut,Flavour * _fl) {
  if (energymin != 0) return;
  ncut      = _ncut;
  fl        = _fl;
  double E  = AORGTOOLS::rpa.gen.Ecms();
  energymin = new double[ncut];
  energymax = new double[ncut];
  cosmin    = new double*[ncut];
  cosmax    = new double*[ncut];
  scut      = new double*[ncut];
  scut_save = new double*[ncut];

  for (int i=0;i<ncut;i++) {
    cosmin[i]      = new double[ncut];
    cosmax[i]      = new double[ncut];
    scut[i]        = new double[ncut];
    scut_save[i]   = new double[ncut];
    energymin[i]   = AMATOOLS::Max(0.,fl[i].Mass());
    energymax[i]   = E;
  }

  for (int i=0;i<ncut;i++) {
    for (int j=i;j<ncut;j++) {
      cosmin[i][j] = cosmin[j][i] = -1.;
      cosmax[i][j] = cosmax[j][i] =  1.;
      double sc =
	+AMATOOLS::sqr(fl[i].Mass())+AMATOOLS::sqr(fl[j].Mass())
	+2.*energymin[i]*energymin[j]
	-2.*sqrt(AMATOOLS::dabs(AMATOOLS::sqr(energymin[i])-AMATOOLS::sqr(fl[i].Mass())))
	*sqrt(AMATOOLS::dabs(AMATOOLS::sqr(energymin[j])-AMATOOLS::sqr(fl[j].Mass())))
	*cosmax[i][j];
      scut[i][j] = scut[j][i] = scut_save[i][j] =
	AMATOOLS::Max(sc,1.e-12*AMATOOLS::sqr(AORGTOOLS::rpa.gen.Ecms()));
    }
  }  
}

void Cut_Data::Init(Cut_Data * _cuts,Flavour * _fl) {
  ncut      = _cuts->ncut;
  fl        = _fl;
  double E  = AORGTOOLS::rpa.gen.Ecms();
  if (!energymin) {
    std::cout<<"Initialize new cuts."<<std::endl;
    energymin           = new double[ncut];
    energymax           = new double[ncut];
    cosmin              = new double*[ncut];
    cosmax              = new double*[ncut];
    scut                = new double*[ncut];
    scut_save           = new double*[ncut];
    for (int i=0;i<ncut;i++) {
      energymin[i]      = _cuts->energymin[i];
      energymax[i]      = _cuts->energymax[i];
      cosmin[i]         = new double[ncut];
      cosmax[i]         = new double[ncut];
      scut[i]           = new double[ncut];
      scut_save[i]      = new double[ncut];
      for (int j=i;j<ncut;j++) {
	cosmin[i][j]    = cosmin[j][i]    = _cuts->cosmin[i][j];
	cosmax[i][j]    = cosmax[j][i]    = _cuts->cosmax[i][j];
	scut[i][j]      = scut[j][i]      = _cuts->scut[i][j];
	scut_save[i][j] = scut_save[j][i] = scut[i][j];
      }
    }
  }
  else {
    std::cout<<"Compare cuts."<<std::endl;
    for (int i=0;i<ncut;i++) {
      energymin[i]      = AMATOOLS::Min(energymin[i],_cuts->energymin[i]);
      energymax[i]      = AMATOOLS::Max(energymax[i],_cuts->energymax[i]);
      for (int j=i;j<ncut;j++) {
	cosmin[i][j]    = cosmin[j][i]    = AMATOOLS::Min(cosmin[i][j],_cuts->cosmin[i][j]);
	cosmax[i][j]    = cosmax[j][i]    = AMATOOLS::Max(cosmax[i][j],_cuts->cosmax[i][j]);
	scut[i][j]      = scut[j][i]      = AMATOOLS::Min(scut[i][j],_cuts->scut[i][j]);
	scut_save[i][j] = scut_save[j][i] = scut[i][j];
      }
    }
  }
}




void Cut_Data::Update(double sprime,double y) {
  double E = sqrt(sprime);
  for (int i=0;i<ncut;i++) {
    energymin[i]   = AMATOOLS::Max(0.,fl[i].Mass());
    energymax[i]   = E;
    for (int j=i+1;j<ncut;j++) {
      cosmin[i][j] = cosmin[j][i] = -1.;
      cosmax[i][j] = cosmax[j][i] =  1.;
    }
  }
  for (int i=0;i<ncut;i++) {
    for (int j=i+1;j<ncut;j++) {
      scut[i][j] = scut[j][i] = AMATOOLS::Max(scut_save[i][j],1.e-12*sprime);
    }
  }  
}







