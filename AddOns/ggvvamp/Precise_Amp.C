#include "Precise_Amp.H"
#include <iostream>
#include <cmath>
#include <ginac/ginac.h>
#include <cln/lfloat.h>

#if defined(DOUBLE_PRECISION) && defined(QUAD_PRECISION)

using namespace std;
using namespace GiNaC;

namespace Precision{

    // Precise_Amp methods
    Precise_Amp::Precise_Amp(int Nf, IRscheme ir) :
        maxrel(1e-5), maxabs(1e-20),
        d(Nf,ir), dx12(Nf,ir), dx34(Nf,ir), dx12x34(Nf,ir),
        q(Nf,ir), qx12(Nf,ir), qx34(Nf,ir), qx12x34(Nf,ir)
    #ifdef ARBITRARY_PRECISION
        , a(Nf,ir), ax12(Nf,ir), ax34(Nf,ir), ax12x34(Nf,ir)
    #endif
    {}


    // compute missing A coefficients from the "irreducible" ones and their crossed versions
    template<class T>
    void employ_crossing_relations(T A[nAcoeff+1][nloop][nreim],
			                 const T Ax12[nAcoeff+1][nloop][nreim],
			                 const T Ax34[nAcoeff+1][nloop][nreim],
			                 const T Ax12x34[nAcoeff+1][nloop][nreim])
    {
        for (int l = 0; l < nloop; ++l) {
            for (int m = 0; m < nreim; ++m) {
            A[3][l][m] = Ax12[2][l][m];
            A[6][l][m] = Ax12[5][l][m];
            A[7][l][m] = Ax12[4][l][m];
            A[10][l][m] = -Ax34[8][l][m];
            A[11][l][m] = -Ax34[9][l][m];
            A[12][l][m] = Ax12[9][l][m];
            A[13][l][m] = Ax12[8][l][m];
            A[14][l][m] = -Ax12x34[9][l][m];
            A[15][l][m] = -Ax12x34[8][l][m];
            A[19][l][m] = Ax12[18][l][m];
            A[20][l][m] = Ax12[17][l][m];
            }
        }
    }


    void Precise_Amp::print_statistics(std::ostream& os) const {
        os << "number of events per number of Digits used:" << endl;
        map<long, long>::const_iterator i;
        for (i = eventsWithDigits.begin(); i != eventsWithDigits.end(); ++i) {
            os << "  " << i->first << " -> " << i->second << endl;
        }
        os << "\n";
    }

    // check if two results are compatible within the given precision
    // and copy higher precision result to target array
    template<class T1, class T2>
    bool control_precision(const T1 Aloprec[nAcoeff+1][nloop][nreim],
		                   const T2 Ahiprec[nAcoeff+1][nloop][nreim],
		                   double Aresult[nAcoeff+1][nloop][nreim],
		                   double maxrel, double maxabs) {
        bool stable = true;
        for (int j = 1; j <= nAcoeff; ++j) {
            for (int l = 0; l < nloop; ++l) {
	            for (int r = 0 ; r < nreim; ++r) {
	                double lo = todouble(Aloprec[j][l][r]);
	                double hi = todouble(Ahiprec[j][l][r]);
	                if (fabs(lo-hi) > maxabs && fabs(lo-hi) > maxrel*(fabs(lo)+fabs(hi))) {stable = false;}
	            Aresult[j][l][r] = hi;
	            }
            }
        }
        return stable;
    }

    numeric tomanydigits(double x) {
        return numeric(cln::cl_float(x, cln::default_float_format));
    }



    // Form factor computing function with numerical accuracy check

    void Precise_Amp::compute(double s_d, double t_d, double ma2_d, double mb2_d) {
        using GiNaC::Digits;
        long maxdigits = 500;
        // disable the following two flags for more conservative (slower) evaluations
        bool recycle_gpls = true;
        bool recycle_coeff = true;

        Digits = maxdigits;
        numeric s = tomanydigits(s_d);
        numeric t = tomanydigits(t_d);
        numeric ma2 = tomanydigits(ma2_d);
        numeric mb2 = tomanydigits(mb2_d);
        numeric u = ma2 + mb2 - s - t;


        // run1: double precision GPL, quad precision algebra
        Digits = 17;
        quadtype (*nq)(const GiNaC::ex&) = &ggvvamp<quadtype>::n;
        q.compute(      nq(s), nq(t), nq(ma2), nq(mb2));
        qx12.compute(   nq(s), nq(u), nq(ma2), nq(mb2));
        qx34.compute(   nq(s), nq(u), nq(mb2), nq(ma2));
        qx12x34.compute(nq(s), nq(t), nq(mb2), nq(ma2));
        employ_crossing_relations(q.A, qx12.A, qx34.A, qx12x34.A);



        // run2: double precision algebra (recycling previous GPL evaluations)
        if (recycle_gpls) {
            size_t hi = ggvvamp<double>::rstartalgebraic;
            for (size_t i = 0 ; i <= hi; ++i) {
            d.r[i]       = todouble(q.r[i]);
            dx12.r[i]    = todouble(qx12.r[i]);
            dx34.r[i]    = todouble(qx34.r[i]);
            dx12x34.r[i] = todouble(qx12x34.r[i]);
            }
        }

        d.compute_gpls       = !recycle_gpls;
        dx12.compute_gpls    = !recycle_gpls;
        dx34.compute_gpls    = !recycle_gpls;
        dx12x34.compute_gpls = !recycle_gpls;
        double (*nd)(const GiNaC::ex&) = &ggvvamp<double>::n;
        d.compute(      nd(s), nd(t), nd(ma2), nd(mb2));
        dx12.compute(   nd(s), nd(u), nd(ma2), nd(mb2));
        dx34.compute(   nd(s), nd(u), nd(mb2), nd(ma2));
        dx12x34.compute(nd(s), nd(t), nd(mb2), nd(ma2));
        employ_crossing_relations(d.A, dx12.A, dx34.A, dx12x34.A);

        if (control_precision(d.A, q.A, A, maxrel, maxabs)) {
            ++eventsWithDigits[Digits];
            return;
        }


        // run3: quad precision GPL, quad precision algebra
        Digits = 32;
        //cerr << "ggvvamp: increasing precision to " << Digits << endl;
        q.compute(      nq(s), nq(t), nq(ma2), nq(mb2));
        qx12.compute(   nq(s), nq(u), nq(ma2), nq(mb2));
        qx34.compute(   nq(s), nq(u), nq(mb2), nq(ma2));
        qx12x34.compute(nq(s), nq(t), nq(mb2), nq(ma2));
        employ_crossing_relations(q.A, qx12.A, qx34.A, qx12x34.A);

        bool looks_ok = control_precision(A, q.A, A, maxrel, maxabs);
        if (recycle_coeff && looks_ok) {
            ++eventsWithDigits[Digits];
            return;
        }

    #ifdef ARBITRARY_PRECISION
        // run4 and more: arbitrary precision GPL and algebra
        for (Digits = 64; Digits <= maxdigits; Digits = Digits*2) {
            //cerr << "ggvvamp: increasing precision to " << Digits << endl;
            cln::cl_F (*na)(const GiNaC::ex&) = &ggvvamp<cln::cl_F>::n;
            a.compute(      na(s), na(t), na(ma2), na(mb2));
            ax12.compute(   na(s), na(u), na(ma2), na(mb2));
            ax34.compute(   na(s), na(u), na(mb2), na(ma2));
            ax12x34.compute(na(s), na(t), na(mb2), na(ma2));
            employ_crossing_relations(a.A, ax12.A, ax34.A, ax12x34.A);
            if (control_precision(A, a.A, A, maxrel, maxabs)) {
            ++eventsWithDigits[Digits];
            return;
            }
        }
    #endif
        // no convergence
        throw runtime_error("ggvvamp: could not converge to stable result !");
    /* } // else end */
        std::cout<<"no form factor defined..."<<std::endl;
    }
}
#endif

    // template instantiations
    #ifdef ARBITRARY_PRECISION
    template class ggvvamp<cln::cl_F>;
    #endif
    #ifdef QUAD_PRECISION
    template class ggvvamp<quadtype>;
    #endif
    #ifdef DOUBLE_PRECISION
    template class ggvvamp<double>;
    #endif
