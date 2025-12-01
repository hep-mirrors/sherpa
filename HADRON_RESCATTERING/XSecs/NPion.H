#ifndef HADRON_RESCATTERING_NPION_H
#define HADRON_RESCATTERING_NPION_H

#include "HADRON_RESCATTERING/XSecs/Helper.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include "HADRON_RESCATTERING/XSecs/HPR1R2.H"
#include <math.h>

class NPion 
{
    private:
    double mProton, mNeutron, mPionPlus, mPionMinus, mPionZero; 
    // mProton = 0.938272;
    // mPionPlus = 0.13957;
    const double GevToMB = 0.389379;   // GeV^-2 -> mb
    HADRON_RESCATTERING::HPR1R2  m_hpr1r2;


    public:
    NPion();
    ~NPion();
    //total scattering cross sections
    //pPi
    double pPiMinus(const double & s);
    double pPiPlus(const double & s);

    //npi
    double nPiMinus(const double & s);
    double nPiPlus(const double & s);
};
#endif