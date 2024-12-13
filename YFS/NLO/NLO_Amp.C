#include "YFS/NLO/NLO_Amp.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace YFS;
using namespace std;

NLO_Amp::NLO_Amp(){

}


NLO_Amp::~NLO_Amp(){

}


double NLO_Amp::Calculate(const Vec4D_Vector &p){
	// copied for Janusz Gluza code
	const double pi = 3.141592653589793;
    
    double facini, facfin, facint, e, p1p2, p4p5, p1p4, p2p5, p1p5, p2p4;
    double p1p3, p4p3, p5p3, p2p3, p1p2pos, p4p5pos, p1p4pos, p2p5pos;
    double p1p5pos, p2p4pos, p1p3pos, p4p3pos, p5p3pos, p2p3pos;
    double matSQ1, matSQ2, s45;
    double Aini, Aint, Afin, num1, Ainipos, Aintpos, Afinpos;
    double flux, JcL3, JcL3Pos, flambda, me;
    double ecm, s, m[5], alpha0, betapi;
    Complex FPIs, FPIs45;

    // Initialize values (these should be set according to your context)
    FPIs = 1.0;
    FPIs45 = 1.0;
    e = 1.0; // Example initialization, update with actual values
    me = p[0].Mass();
    double mp = p[2].Mass();
    s = (p[0]+p[1]).Abs2();
    betapi = sqrt(1.-4*mp*mp/s);

    facini = pow(e, 6) * abs(FPIs45) * abs(FPIs45);
    facfin = pow(e, 6) * abs(FPIs) * abs(FPIs);
    facint = pow(e, 6) * 2.0 * (real(FPIs * conj(FPIs45)) + real(conj(FPIs) * FPIs45)) / 2.0;

    p1p2 = p[0]*p[1];
    p4p5 = p[3]*p[4];
    p1p4 = p[0]*p[3];
    p2p5 = p[1]*p[4];
    p1p5 = p[0]*p[4];
    p2p4 = p[1]*p[3];
    p1p3 = p[0]*p[2];
    p4p3 = p[3]*p[2];
    p5p3 = p[4]*p[2];
    p2p3 = p[1]*p[2];
    s45 = (p[3]+p[2]).Abs2();

    if (p1p2 == 0.0 || p4p5 == 0.0 || p1p4 == 0.0 || p2p5 == 0.0 ||
        p1p5 == 0.0 || p2p4 == 0.0 || p1p3 == 0.0 || p4p3 == 0.0 ||
        p5p3 == 0.0 || p2p3 == 0.0) {
        msg_Error()<< "Division by zero in "<< METHOD << endl;
        return -1; // Indicate error
    }

    Aini = - 2 * pow(me, 2) * pow(s45, -2) * p1p2 * p4p5 * pow(p1p3, -2) 
           + 4 * pow(me, 2) * pow(s45, -2) * p1p2 * p4p5 * pow(p1p3, -1) * pow(p2p3, -1)
           - 2 * pow(me, 2) * pow(s45, -2) * p1p2 * p4p5 * pow(p2p3, -2) 
           + 2 * pow(me, 2) * pow(s45, -2) * p1p4 * p2p5 * pow(p1p3, -2) 
           + 2 * pow(me, 2) * pow(s45, -2) * p1p4 * p2p5 * pow(p2p3, -2) 
           - 2 * pow(me, 2) * pow(s45, -2) * p1p4 * p2p4 * pow(p1p3, -2) 
           - 2 * pow(me, 2) * pow(s45, -2) * p1p4 * p2p4 * pow(p2p3, -2) 
           - 2 * pow(me, 2) * pow(s45, -2) * p2p5 * p1p5 * pow(p1p3, -2) 
           - 2 * pow(me, 2) * pow(s45, -2) * p2p5 * p1p5 * pow(p2p3, -2) 
           + 2 * pow(me, 2) * pow(s45, -2) * p1p5 * p2p4 * pow(p1p3, -2) 
           + 2 * pow(me, 2) * pow(s45, -2) * p1p5 * p2p4 * pow(p2p3, -2) 
           - 0.5 * pow(me, 2) * pow(s45, -1) * pow(betapi, 2) * p1p2 * pow(p1p3, -2) 
           + pow(me, 2) * pow(s45, -1) * pow(betapi, 2) * p1p2 * pow(p1p3, -1) * pow(p2p3, -1) 
           - 0.5 * pow(me, 2) * pow(s45, -1) * pow(betapi, 2) * p1p2 * pow(p2p3, -2) 
           + 0.5 * pow(me, 2) * pow(s45, -1) * p1p2 * pow(p1p3, -2) 
           - pow(me, 2) * pow(s45, -1) * p1p2 * pow(p1p3, -1) * pow(p2p3, -1) 
           + 0.5 * pow(me, 2) * pow(s45, -1) * p1p2 * pow(p2p3, -2) 
           - 2 * pow(me, 4) * pow(s45, -2) * p4p5 * pow(p1p3, -2) 
           - 2 * pow(me, 4) * pow(s45, -2) * p4p5 * pow(p2p3, -2) 
           - 0.5 * pow(me, 4) * pow(s45, -1) * pow(betapi, 2) * pow(p1p3, -2) 
           - 0.5 * pow(me, 4) * pow(s45, -1) * pow(betapi, 2) * pow(p2p3, -2) 
           + 0.5 * pow(me, 4) * pow(s45, -1) * pow(p1p3, -2) 
           + 0.5 * pow(me, 4) * pow(s45, -1) * pow(p2p3, -2) 
           - 4 * pow(s45, -2) * p1p2 * p1p4 * p2p5 * pow(p1p3, -1) * pow(p2p3, -1) 
           + 4 * pow(s45, -2) * p1p2 * p1p4 * p2p4 * pow(p1p3, -1) * pow(p2p3, -1) 
           + 4 * pow(s45, -2) * p1p2 * p2p5 * p1p5 * pow(p1p3, -1) * pow(p2p3, -1) 
           - 4 * pow(s45, -2) * p1p2 * p1p5 * p2p4 * pow(p1p3, -1) * pow(p2p3, -1) 
           + 4 * pow(s45, -2) * pow(p1p2, 2) * p4p5 * pow(p1p3, -1) * pow(p2p3, -1) 
           + pow(s45, -1) * pow(betapi, 2) * pow(p1p2, 2) * pow(p1p3, -1) * pow(p2p3, -1) 
           - pow(s45, -1) * pow(p1p2, 2) * pow(p1p3, -1) * pow(p2p3, -1);

    Aint =  - 2 * pow(me, 2) * pow(s, -1) * pow(s45, -1) * p4p5 * p1p4 * pow(p1p3, -1) * pow(p4p3, -1)
            - 2 * pow(me, 2) * pow(s, -1) * pow(s45, -1) * p4p5 * p2p5 * pow(p5p3, -1) * pow(p2p3, -1)
            + 2 * pow(me, 2) * pow(s, -1) * pow(s45, -1) * p4p5 * p1p5 * pow(p1p3, -1) * pow(p5p3, -1)
            + 2 * pow(me, 2) * pow(s, -1) * pow(s45, -1) * p4p5 * p2p4 * pow(p4p3, -1) * pow(p2p3, -1)
            - 0.5 * pow(me, 2) * pow(s, -1) * pow(betapi, 2) * p1p4 * pow(p1p3, -1) * pow(p4p3, -1)
            - 0.5 * pow(me, 2) * pow(s, -1) * pow(betapi, 2) * p2p5 * pow(p5p3, -1) * pow(p2p3, -1)
            + 0.5 * pow(me, 2) * pow(s, -1) * pow(betapi, 2) * p1p5 * pow(p1p3, -1) * pow(p5p3, -1)
            + 0.5 * pow(me, 2) * pow(s, -1) * pow(betapi, 2) * p2p4 * pow(p4p3, -1) * pow(p2p3, -1)
            + 0.5 * pow(me, 2) * pow(s, -1) * p1p4 * pow(p1p3, -1) * pow(p4p3, -1)
            + 0.5 * pow(me, 2) * pow(s, -1) * p2p5 * pow(p5p3, -1) * pow(p2p3, -1)
            - 0.5 * pow(me, 2) * pow(s, -1) * p1p5 * pow(p1p3, -1) * pow(p5p3, -1)
            - 0.5 * pow(me, 2) * pow(s, -1) * p2p4 * pow(p4p3, -1) * pow(p2p3, -1)
            - 2 * pow(s, -1) * pow(s45, -1) * p1p2 * p4p5 * p1p4 * pow(p1p3, -1) * pow(p4p3, -1)
            - 2 * pow(s, -1) * pow(s45, -1) * p1p2 * p4p5 * p2p5 * pow(p5p3, -1) * pow(p2p3, -1)
            + 2 * pow(s, -1) * pow(s45, -1) * p1p2 * p4p5 * p1p5 * pow(p1p3, -1) * pow(p5p3, -1)
            + 2 * pow(s, -1) * pow(s45, -1) * p1p2 * p4p5 * p2p4 * pow(p4p3, -1) * pow(p2p3, -1);

    Aint = Aint 
            - 2 * pow(s, -1) * pow(s45, -1) * p1p4 * p2p5 * p1p5 * pow(p1p3, -1) * pow(p4p3, -1)
            - 2 * pow(s, -1) * pow(s45, -1) * p1p4 * p2p5 * p1p5 * pow(p1p3, -1) * pow(p5p3, -1)
            - 2 * pow(s, -1) * pow(s45, -1) * p1p4 * p2p5 * p2p4 * pow(p4p3, -1) * pow(p2p3, -1)
            - 2 * pow(s, -1) * pow(s45, -1) * p1p4 * p2p5 * p2p4 * pow(p5p3, -1) * pow(p2p3, -1)
            + 2 * pow(s, -1) * pow(s45, -1) * p1p4 * pow(p2p5, 2) * pow(p5p3, -1) * pow(p2p3, -1)
            + 2 * pow(s, -1) * pow(s45, -1) * p1p4 * p1p5 * p2p4 * pow(p1p3, -1) * pow(p4p3, -1)
            + 2 * pow(s, -1) * pow(s45, -1) * p1p4 * p1p5 * p2p4 * pow(p1p3, -1) * pow(p5p3, -1)
            + 2 * pow(s, -1) * pow(s45, -1) * p1p4 * pow(p2p4, 2) * pow(p4p3, -1) * pow(p2p3, -1)
            + 2 * pow(s, -1) * pow(s45, -1) * pow(p1p4, 2) * p2p5 * pow(p1p3, -1) * pow(p4p3, -1)
            - 2 * pow(s, -1) * pow(s45, -1) * pow(p1p4, 2) * p2p4 * pow(p1p3, -1) * pow(p4p3, -1)
            + 2 * pow(s, -1) * pow(s45, -1) * p2p5 * p1p5 * p2p4 * pow(p4p3, -1) * pow(p2p3, -1)
            + 2 * pow(s, -1) * pow(s45, -1) * p2p5 * p1p5 * p2p4 * pow(p5p3, -1) * pow(p2p3, -1)
            + 2 * pow(s, -1) * pow(s45, -1) * p2p5 * pow(p1p5, 2) * pow(p1p3, -1) * pow(p5p3, -1)
            - 2 * pow(s, -1) * pow(s45, -1) * pow(p2p5, 2) * p1p5 * pow(p5p3, -1) * pow(p2p3, -1)
            - 2 * pow(s, -1) * pow(s45, -1) * p1p5 * pow(p2p4, 2) * pow(p4p3, -1) * pow(p2p3, -1);

    Aint = Aint 
            - 2 * pow(s, -1) * pow(s45, -1) * pow(p1p5, 2) * p2p4 * pow(p1p3, -1) * pow(p5p3, -1)
            - 0.5 * pow(s, -1) * pow(betapi, 2) * p1p2 * p1p4 * pow(p1p3, -1) * pow(p4p3, -1)
            - 0.5 * pow(s, -1) * pow(betapi, 2) * p1p2 * p2p5 * pow(p5p3, -1) * pow(p2p3, -1)
            + 0.5 * pow(s, -1) * pow(betapi, 2) * p1p2 * p1p5 * pow(p1p3, -1) * pow(p5p3, -1)
            + 0.5 * pow(s, -1) * pow(betapi, 2) * p1p2 * p2p4 * pow(p4p3, -1) * pow(p2p3, -1)
            + 0.5 * pow(s, -1) * p1p2 * p1p4 * pow(p1p3, -1) * pow(p4p3, -1)
            + 0.5 * pow(s, -1) * p1p2 * p2p5 * pow(p5p3, -1) * pow(p2p3, -1)
            - 0.5 * pow(s, -1) * p1p2 * p1p5 * pow(p1p3, -1) * pow(p5p3, -1)
            - 0.5 * pow(s, -1) * p1p2 * p2p4 * pow(p4p3, -1) * pow(p2p3, -1);

    Afin =  + 0.5 * pow(me, 2) * pow(s, -2) * s45 * pow(betapi, 2) * p4p5 * pow(p4p3, -2)
            + pow(me, 2) * pow(s, -2) * s45 * pow(betapi, 2) * p4p5 * pow(p4p3, -1) * pow(p5p3, -1)
            + 0.5 * pow(me, 2) * pow(s, -2) * s45 * pow(betapi, 2) * p4p5 * pow(p5p3, -2)
            - 0.5 * pow(me, 2) * pow(s, -2) * s45 * p4p5 * pow(p4p3, -2)
            - pow(me, 2) * pow(s, -2) * s45 * p4p5 * pow(p4p3, -1) * pow(p5p3, -1)
            - 0.5 * pow(me, 2) * pow(s, -2) * s45 * p4p5 * pow(p5p3, -2)
            - 0.25 * pow(me, 2) * pow(s, -2) * pow(s45, 2) * pow(betapi, 2) * pow(p4p3, -2)
            - 0.25 * pow(me, 2) * pow(s, -2) * pow(s45, 2) * pow(betapi, 2) * pow(p5p3, -2)
            + 0.125 * pow(me, 2) * pow(s, -2) * pow(s45, 2) * pow(betapi, 4) * pow(p4p3, -2)
            + 0.125 * pow(me, 2) * pow(s, -2) * pow(s45, 2) * pow(betapi, 4) * pow(p5p3, -2)
            + 0.125 * pow(me, 2) * pow(s, -2) * pow(s45, 2) * pow(p4p3, -2)
            + 0.125 * pow(me, 2) * pow(s, -2) * pow(s45, 2) * pow(p5p3, -2)
            + 4 * pow(me, 2) * pow(s, -2) * pow(p4p5, 2) * pow(p4p3, -1) * pow(p5p3, -1)
            + 0.5 * pow(s, -2) * s45 * pow(betapi, 2) * p1p2 * p4p5 * pow(p4p3, -2)
            + pow(s, -2) * s45 * pow(betapi, 2) * p1p2 * p4p5 * pow(p4p3, -1) * pow(p5p3, -1)
            + 0.5 * pow(s, -2) * s45 * pow(betapi, 2) * p1p2 * p4p5 * pow(p5p3, -2)
            - 0.5 * pow(s, -2) * s45 * pow(betapi, 2) * p1p4 * p2p5 * pow(p4p3, -2)
            - 0.5 * pow(s, -2) * s45 * pow(betapi, 2) * p1p4 * p2p5 * pow(p5p3, -2);

    Afin = Afin 
            + 0.5 * pow(s, -2) * s45 * pow(betapi, 2) * p1p4 * p2p4 * pow(p4p3, -2)
            + 0.5 * pow(s, -2) * s45 * pow(betapi, 2) * p1p4 * p2p4 * pow(p5p3, -2)
            + 0.5 * pow(s, -2) * s45 * pow(betapi, 2) * p2p5 * p1p5 * pow(p4p3, -2)
            + 0.5 * pow(s, -2) * s45 * pow(betapi, 2) * p2p5 * p1p5 * pow(p5p3, -2)
            - 0.5 * pow(s, -2) * s45 * pow(betapi, 2) * p1p5 * p2p4 * pow(p4p3, -2)
            - 0.5 * pow(s, -2) * s45 * pow(betapi, 2) * p1p5 * p2p4 * pow(p5p3, -2)
            - 0.5 * pow(s, -2) * s45 * p1p2 * p4p5 * pow(p4p3, -2)
            - pow(s, -2) * s45 * p1p2 * p4p5 * pow(p4p3, -1) * pow(p5p3, -1)
            - 0.5 * pow(s, -2) * s45 * p1p2 * p4p5 * pow(p5p3, -2)
            + 0.5 * pow(s, -2) * s45 * p1p4 * p2p5 * pow(p4p3, -2)
            + 0.5 * pow(s, -2) * s45 * p1p4 * p2p5 * pow(p5p3, -2)
            - 0.5 * pow(s, -2) * s45 * p1p4 * p2p4 * pow(p4p3, -2)
            - 0.5 * pow(s, -2) * s45 * p1p4 * p2p4 * pow(p5p3, -2)
            - 0.5 * pow(s, -2) * s45 * p2p5 * p1p5 * pow(p4p3, -2)
            - 0.5 * pow(s, -2) * s45 * p2p5 * p1p5 * pow(p5p3, -2)
            + 0.5 * pow(s, -2) * s45 * p1p5 * p2p4 * pow(p4p3, -2)
            + 0.5 * pow(s, -2) * s45 * p1p5 * p2p4 * pow(p5p3, -2)
            - 0.25 * pow(s, -2) * pow(s45, 2) * pow(betapi, 2) * p1p2 * pow(p4p3, -2)
            - 0.25 * pow(s, -2) * pow(s45, 2) * pow(betapi, 2) * p1p2 * pow(p5p3, -2)
            + 0.125 * pow(s, -2) * pow(s45, 2) * pow(betapi, 4) * p1p2 * pow(p4p3, -2)
            + 0.125 * pow(s, -2) * pow(s45, 2) * pow(betapi, 4) * p1p2 * pow(p5p3, -2);

    Afin = Afin 
            + 0.125 * pow(s, -2) * pow(s45, 2) * p1p2 * pow(p4p3, -2)
            + 0.125 * pow(s, -2) * pow(s45, 2) * p1p2 * pow(p5p3, -2)
            + 4 * pow(s, -2) * p1p2 * pow(p4p5, 2) * pow(p4p3, -1) * pow(p5p3, -1)
            - 4 * pow(s, -2) * p4p5 * p1p4 * p2p5 * pow(p4p3, -1) * pow(p5p3, -1)
            + 4 * pow(s, -2) * p4p5 * p1p4 * p2p4 * pow(p4p3, -1) * pow(p5p3, -1)
            + 4 * pow(s, -2) * p4p5 * p2p5 * p1p5 * pow(p4p3, -1) * pow(p5p3, -1)
            - 4 * pow(s, -2) * p4p5 * p1p5 * p2p4 * pow(p4p3, -1) * pow(p5p3, -1);

   matSQ1= Aini+Aint+Afin;    
   return m_alpha/2./pow(2 * M_PI, 3)*matSQ1/s/s;
}