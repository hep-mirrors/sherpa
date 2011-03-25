# plot EPA spectra

set term postscript enhanced color solid
set output 'EPA_spectrum.ps'
set nogrid
set logscale x
set logscale y
set xlabel 'x_{/Symbol w}={/Symbol w}/{/Symbol w}_0'
set ylabel 'x_{/Symbol w}*n(x_{/Symbol w}) [GeV^{-1}]'
set key inside left bottom

set pointsize .5

# fits used to smooth hcs
f(x)=c*x**a
a=-0.059967
c=1.36549
#g(x)=d*x**b
# fit for x_omega not useful
# used EPA::selfTest with finer steps to determine cut-off:
# x_omega = 1.33086 with x_omega*n(x_omega)=9.924e-6

fit  [1.4e-3:5e-3] f(x) 'EPA_DebugOut_PionsHCSGauss1.log' via a,c
#fit  [1.27:1.3288] g(x) 'EPA_DebugOut_PionsHCSGauss1.log' via b,d

plot [1e-3:5] [1e-3:10] 'EPA_DebugOut_Pions1.log' \
                        with points title 'point form', \
     	      		'EPA_DebugOut_PionsHCSGauss1.log' \
			with points title 'hcs form', \
			'EPA_DebugOut_PionsHCSGauss2.log' \
			with points title 'gauss form', \
			'EPA_DebugOut_PionsSmoothHCS2.log' \
			with points title 'smooth hcs form', \
			f(x) \
			title 'fit for smoothing hcs at low x_{/Symbol w}'

#			g(x) \
#			title 'fit for smoothing hcs at high x_{/Symbol w}'
