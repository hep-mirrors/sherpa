#include "Lorentz_Function.H"

using namespace AMEGIC;

void Lorentz_Function::AddPermutation(int sign,int a,int b=-1,int c=-1,int d=-1)
{
  int* newperm = new int[NofIndex()];
  newperm[0] = partarg[a];
  if (NofIndex()>1) newperm[1] = partarg[b];
  if (NofIndex()>2) newperm[2] = partarg[c];
  if (NofIndex()>3) newperm[3] = partarg[d];

  permlist.push_back(newperm);
  signlist.push_back(sign);
}

void Lorentz_Function::InitPermutation()
{
  if (!permlist.empty()) {
    for (short int i=0;i<permlist.size();i++) delete[] permlist[i]; 
    permlist.clear();
    signlist.clear();
  }

  switch (type) {
  case lf::Gab   : 
    AddPermutation(1,0,1);
    AddPermutation(1,1,0);  
    break;
  case lf::Gauge3:
    AddPermutation( 1,0,1,2);
    AddPermutation(-1,0,2,1);  
    AddPermutation(-1,1,0,2);
    AddPermutation(-1,2,1,0);  
    AddPermutation( 1,1,2,0);
    AddPermutation( 1,2,0,1);  
    break;
  case lf::Gauge4: 
    AddPermutation( 1,0,1,2,3);
    AddPermutation( 1,1,0,2,3);
    AddPermutation( 1,0,1,3,2);
    AddPermutation( 1,1,0,3,2);
    AddPermutation( 1,2,3,1,0);    
    AddPermutation( 1,3,2,1,0);
    AddPermutation( 1,2,3,0,1);
    AddPermutation( 1,3,2,0,1);
    break;
  case lf::Gluon4: 
    AddPermutation( 1,0,1,2,3);
    AddPermutation(-1,2,1,0,3);
    AddPermutation(-1,0,3,2,1);
    AddPermutation( 1,2,3,0,1);
    AddPermutation( 1,1,0,3,2);
    AddPermutation(-1,3,0,1,2);
    AddPermutation(-1,1,2,3,0);
    AddPermutation( 1,3,2,1,0);
    break;
  case lf::SSV   : 
    AddPermutation( 1,0,1,2);
    AddPermutation(-1,1,0,2);
    break;
  case lf::Photon3_NC:
    AddPermutation( 1,0,1,2);
    AddPermutation(-1,0,2,1);  
    AddPermutation(-1,1,0,2);
    AddPermutation(-1,2,1,0);  
    AddPermutation( 1,1,2,0);
    AddPermutation( 1,2,0,1);  
    break;
  case lf::Photon4_NC:
    //still to be done
    break;
  }
  permcount = 0;
}

int Lorentz_Function::ResetPermutation() 
{
  permcount=0;
  for (short int i=0;i<NofIndex();i++) partarg[i]  = permlist[permcount][i];
}

int Lorentz_Function::NextPermutation()
{
  if (NofIndex()<2) return 0;
  permcount++;
  if (permcount==permlist.size()) return 0;
  
  for (short int i=0;i<NofIndex();i++) partarg[i]  = permlist[permcount][i];
					 
}

int Lorentz_Function::GetSign() 
{
  if (signlist.empty()) return 1;
  return signlist[permcount];
}









