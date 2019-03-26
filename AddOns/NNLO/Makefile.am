
GITTAG = AddOns/NNLO
include ../../git.make

pkglib_LTLIBRARIES    =	libNNLOqT.la

libNNLOqT_la_SOURCES = \
	Git_Info.C \
	polylog.c \
	Tools.C \
	QT_Selector.C \
	DIS_Selector.C \
	Scale_Setter.C \
	coeffqt.C \
	param.C \
	hard.C \
        DY_KFactor.C \
	higgsfullsm.C \
	H_KFactor.C \
	DISinclu.C \
	DIS_KFactor.C \
	DIS_Scale.C

NNLOqT_EXTRA_DIST = \
	Tools.H \
	QT_Selector.H \
	coeffqt.H \
	param.H \
	hard.H \
	higgsfullsm.H \
	DISinclu.H

EXTRA_DIST = $(NNLOqT_EXTRA_DIST)
