#include "ggvvamp.h"
#include <map>

#if defined(DOUBLE_PRECISION) && defined(QUAD_PRECISION)

namespace Precision{
    class Precise_Amp {
    public:
        Precise_Amp(int Nf, IRscheme ir);
        void compute(double s, double t, double ma2, double mb2);
        void print_statistics(std::ostream&) const;
        double maxrel; // maximal relative target uncertainty
        double maxabs; // maximal absolute target uncertainty
        std::map<long, long> eventsWithDigits; // counters for some statistics
        double A[nAcoeff+1][nloop][nreim];
    private:
        ggvvamp<double>    d, dx12, dx34, dx12x34; // double precision evaluation
        ggvvamp<quadtype>  q, qx12, qx34, qx12x34; // quad precision evaluation
#ifdef ARBITRARY_PRECISION
        ggvvamp<cln::cl_F> a, ax12, ax34, ax12x34; // arbitrary precision evaluation
#endif
    };
}

# endif
