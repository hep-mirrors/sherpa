/**
 * @file
 * @author  Hannes Jung <hannes.jung@desy.de>
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * https://www.gnu.org/copyleft/gpl.html
 *
 *
 */


#ifndef TMDlib_H
#define TMDlib_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <math.h>
#include <vector>
#include <iostream>
#include <map>

#include "factories.h"

//the following are UBUNTU/LINUX, and MacOS ONLY terminal color codes.
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */





namespace TMDlib {

double ipow(double,int);
extern int iipow(int,int);


class TMD {
public:
    TMD() ;  // This is the constructor
int icount_d = 0;

    void TMDinit(const std::string name, int irep, int imode  ) ;
    void TMDinit(const std::string name, int irep ) ;
    void TMDinit(const std::string name) ;
    void TMDinit(int & iset) ;

    /// Strip leading and trailing spaces (not in-place) from LHAPDF
    inline std::string trim(const std::string& s) {
        const size_t firstnonspacepos = s.find_first_not_of(" ");
        const size_t lastnonspacepos = s.find_last_not_of(" ");
        if (firstnonspacepos == std::string::npos) return "";
        return s.substr(firstnonspacepos, lastnonspacepos-firstnonspacepos+1);
    }

    ~TMD()
    {
        icount_d = icount_d + 1;
        try
        {
            if (f_grid0)
            {
                delete[] f_grid0;
                std::cout << "deallocate f_grid0 " << icount_d << std::endl;
            }
            if (f_grid1)
            {
                delete[] f_grid1;
                std::cout << "deallocate f_grid1" << std::endl;
            }
            if (f_grid2)
            {
                delete[] f_grid2;
                std::cout << "deallocate f_grid2" << std::endl;
            }
            if (f_grid3)
            {
                delete[] f_grid3;
                std::cout << "deallocate f_grid3" << std::endl;
            }
            if (f_grid4)
            {
                delete[] f_grid4;
                std::cout << "deallocate f_grid4" << std::endl;
            }
            if (f_grid5)
            {
                delete[] f_grid5;
                std::cout << "deallocate f_grid5" << std::endl;
            }
            if (f_grid6)
            {
                delete[] f_grid6;
                std::cout << "deallocate f_grid6" << std::endl;
            }
            if (f_grid7)
            {
                delete[] f_grid7;
                std::cout << "deallocate f_grid7" << std::endl;
            }
            if (f_grid8)
            {
                delete[] f_grid8;
                std::cout << "deallocate f_grid8" << std::endl;
            }
            if (f_grid9)
            {
                delete[] f_grid9;
                std::cout << "deallocate f_grid9" << std::endl;
            }
            if (f_grid10)
            {
                delete[] f_grid10;
                std::cout << "deallocate f_grid10" << std::endl;
            }
            if (f_grid11)
            {
                delete[] f_grid11;
                std::cout << "deallocate f_grid11" << std::endl;
            }
            // Added missing deletions:
            if (f_grid1m)
            {
                delete[] f_grid1m;
                std::cout << "deallocate f_grid1m" << std::endl;
            }
            if (f_grid2m)
            {
                delete[] f_grid2m;
                std::cout << "deallocate f_grid2m" << std::endl;
            }
            if (f_grid3m)
            {
                delete[] f_grid3m;
                std::cout << "deallocate f_grid3m" << std::endl;
            }
            if (f_grid4m)
            {
                delete[] f_grid4m;
                std::cout << "deallocate f_grid4m" << std::endl;
            }
            if (f_grid5m)
            {
                delete[] f_grid5m;
                std::cout << "deallocate f_grid5m" << std::endl;
            }
            if (f_grid6m)
            {
                delete[] f_grid6m;
                std::cout << "deallocate f_grid6m" << std::endl;
            }
            if (a)
            {
                delete[] a;
                std::cout << "deallocate a" << std::endl;
            }
            if (xa)
            {
                delete[] xa;
                std::cout << "deallocate xa" << std::endl;
            }
            if (px)
            {
                delete[] px;
                std::cout << "deallocate px" << std::endl;
            }
            if (q2x)
            {
                delete[] q2x;
                std::cout << "deallocate q2x" << std::endl;
            }
            if (xx)
            {
                delete[] xx;
            }
        }
        catch(const std::exception& e)
        {
            std::cout << "crash " << e.what();
        }
    }
    


    // TMDpdf 1st to return TMD value at x,xbar,kt,mu as a vector
    std::vector<double> TMDpdf( double x, double xbar, double kt, double mu) ;
    // TMDpdf 3rd to return TMD value at x,xbar,kt,mu as a vector
    void TMDpdf( double x, double xbar, double kt, double mu, std::vector<double>& xpq ) ;
    // TMDpdf 5th to return TMD value at x,xbar,kt,mu for flavor
    void TMDpdf( double x, double xbar, double kt, double mu, double& uval, double& dval, double& sea, double& charm, double& bottom, double& gluon,  double& photon);
    // TMDpdf 7th to return TMD value at x,xbar,kt,mu for flavor
    void TMDpdf( double x, double xbar, double kt, double mu, double& up, double& ubar, double& dn, double& dbar, double& strange,  double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& gluon,  double& photon );
    // TMDpdf 9th to return TMD value at x,xbar,kt,mu for flavor
    void TMDpdf( double x, double xbar, double kt, double mu, double& up, double& ubar, double& dn, double& dbar, double& strange,  double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& gluon,  double& photon, double& z0, double&  wplus,  double& wminus,  double& higgs);
    // TMDpdf 11th to return TMD value at x,xbar,kt,mu for flavor
    void TMDpdf( double x, double xbar, double kt, double mu, double& up, double& ubar, double& dn, double& dbar, double& strange,  double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& top, double& tbar, double& gluon,  double& photon, double& z0, double&  wplus,  double& wminus,  double& higgs);
    // main TMD routine (short)
    void TMDdensity( double x, double xbar, double kt, double mu, double& up, double& ubar, double& dn, double& dbar, double& strange,  double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& gluon, double& photon, double& z0, double&  wplus,  double& wminus,  double& higgs);

    // main TMD routine (full)
    void TMDdensity( double x, double xbar, double kt, double mu, double& up, double& ubar, double& dn, double& dbar, double& strange,  double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& top, double& tbar, double& gluon, double& photon, double& z0, double&  wplus,  double& wminus,  double& higgs);

    double TMDalphas(double mu);

    double TMDgetLam4( );
    int TMDgetNf( );
    int TMDgetNumMembers( );
    int TMDgetOrderAlphaS( );
    int TMDgetOrderPDF( );
    int TMDverbosity();
    double TMDgetXmin();
    double TMDgetXmax();
    double TMDgetQ2min();
    double TMDgetQ2max();
    double TMDgetQmin();
    double TMDgetQmax();
    double TMDgetKtmin();
    double TMDgetKtmax();
    std::string TMDgetDesc();
    std::string TMDgetScheme();
    std::string TMDgetIndex();
    int TMDnumberPDF(std::string name);
    std::string TMDstringPDF(int index);
    void TMDinfo(const std::string name);
    std::string TMDgetExtrapolation_Q2();
    std::string TMDgetExtrapolation_x();
    std::string TMDgetExtrapolation_kt();

    bool DoesFileExist (const std::string name);

    double get_key_val_as_double(const std::string name);
    int get_key_val_as_int(const std::string name);


    //double Cdhfint(int narg, double* arg[500], int* nent[500], double* ent[500], double* table[1000000]);
    //double Cdhfint(int narg, double arg[500], int nent[500], double ent[500], double table[1000000]);
    double Cdhfint(int narg, double arg[], int nent[], double ent[], double table[]);
    void polint(double xa[], double ya[], int n, double x, double& y, double& dy) ;
    void polin3(double x1a[], double x2a[], double x3a[], double ya[][1][4], int m, int n, int o, double x1, double x2, double x3, double& y, double& dy);
// PB stuff
    void allFlavuPDF(double x, double kt, double p, double& up, double& ubar, double& dn, double& dbar, double& strange, double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& top, double& tbar, double& phot, double& gluon ) ;
    void allFlavuPDF_old(double x, double kt, double p, double& up, double& ubar, double& dn, double& dbar, double& strange, double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& top, double& tbar, double& phot, double& gluon ) ;
    void allFlavuPDFn(double x, double kt, double p, double& up, double& ubar, double& dn, double& dbar, double& strange, double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& top, double& tbar, double& phot, double& gluon ) ;
    void allFlavuPDFew(double x, double kt, double p, double& up, double& ubar, double& dn, double& dbar, double& strange, double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& top, double& tbar, double& phot, double& gluon, double& z0, double&  wplus, double& wminus, double& higgs) ;
    void ccfm_gluon(double x, double kt, double p, double& up, double& ubar, double& dn, double& dbar, double& strange, double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& top, double& tbar, double& phot, double& gluon ) ;

// GBW
    void gbwuPDFlight(double x, double kt, double& gluon ) ;
    void gbwuPDFcharm(double x, double kt, double& gluon ) ;
    // double blueml(double x, double kt, double p ) ;


// SBRS
    void sbrsPDF(int irep, int imode, double x, double kt, double p, double& uval, double& dval, double& sea, double& charm);

// Kutak Sapeta
    void ksPDF(double x, double kt, double p, double& upl, double& dn, double& sea, double& charm, double& bottom, double& gluon);
    void ksuPDFgrid(int kf, double x, double kt, double p, double& up, double& dn, double& sea, double& charm, double& bottom, double& gluon ) {};
    void ksDLCPDF(double x, double kt, double p, double& up, double& ubar, double& down, double& dbar, double& strange, double& sbar, double& charm,double& cbar, double& bottom, double& bbar, double& gluon);
    void ksBHKSPDF(double x, double kt, double p, double& up, double& ubar, double& down, double& dbar, double& strange, double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& gluon);

// M. Echevaria TMDgluon
    void tmd_ME(double x, double kt, double p, double& uval, double& dval, double& sea, double& charm, double& bottom, double& gluon ) ;

// Pavia TMDPDFs and TMDFFs
    void Pavia(double x, double kt, double p, double& up, double& ubar, double& dn, double& dbar, double& strange, double& sbar, double& charm, double& cbar, double& bottom, double& bbar, double& top, double& tbar, double& phot, double& gluon ) ;


    void setVerbosity(int noiselevel);

private:
    int iset, irep, imode ;
    int noiselevel ;
    // const std::string pdfpath;
    bool first ;
    int ncallTMDdensity;
    TMDlib::TMDGrid* TMDs;

    int NewFormat;
    std::string TMD_Name, TMD_Dir, TMD_Mem;

    double qscal, qg0, q2;
    int n1min, n1max, n2min, n2max, n3min, n3max, ncall;
    int iqqbar, iglu, ikincut, ipgg, ns_sel;
    double* a = new double[153];
    double* xa= new double[3];
    double* px = new double[52];
    double* xx = new double[52];
    double* q2x = new double[52];
    double* f_grid0 = new double[132651];
    double* f_grid1 = new double[132651];
    double* f_grid2 = new double[132651];
    double* f_grid3 = new double[132651];
    double* f_grid4 = new double[132651];
    double* f_grid5 = new double[132651];
    double* f_grid6 = new double[132651];
    double* f_grid7 = new double[132651];
    double* f_grid1m = new double[132651];
    double* f_grid2m = new double[132651];
    double* f_grid3m = new double[132651];
    double* f_grid4m = new double[132651];
    double* f_grid5m = new double[132651];
    double* f_grid6m = new double[132651];
    double* f_grid8 = new double[132651];
    double* f_grid9 = new double[132651];
    double* f_grid10 = new double[132651];
    double* f_grid11 = new double[132651];

    std::map<std::string,std::string> TMDdict;
    std::map<int, std::string> TMDindex;
   

};


extern "C" {
    // Fortran routines
    // remember FORTRAN passes pointer (instead of values) to subroutines
#define initalphas initalphas_
    void initalphas(int* iord, double* fr2, double* mur, double* asmur, double* mc, double* mb, double* mt);
#define unpol unpol_
    void unpol(int* in, double* scale, double* scale2, double* bc, double* lam1, double* lam2, double* x_ME, double* kt_ME,double* Qf, double* result );
#define blueml blueml_
    double blueml(double* x, double* kt, double* p ) ;
}



}

#endif
