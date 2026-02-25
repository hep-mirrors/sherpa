#ifndef LDME_H
#define LDME_H

#include "ATOOLS/Org/Run_Parameter.H"
#include <unordered_map>

void msgLDME(kf_code kfc, double LDME);
void LoadLDME();
void WarnZeroLDME(kf_code kfc);
double GetLDME(kf_code kfc);

extern std::unordered_map<kf_code, double> ldme_map;



#endif