
SVNTAG = AddOns/NNLO
include ../../svn.make

pkglib_LTLIBRARIES    =	libNNLOqT.la

libNNLOqT_la_SOURCES = \
	SVN_Info.C \
	polylog.c \
	Tools.C \
	QT_Selector.C \
	Scale_Setter.C \
	coeffqt.C \
	param.C \
	hard.C \
        DY_KFactor.C \
	higgsfullsm.C \
	H_KFactor.C

NNLOqT_EXTRA_DIST = \
	Tools.H \
	QT_Selector.H \
	coeffqt.H \
	param.H \
	hard.H \
	higgsfullsm.H

EXTRA_DIST = $(NNLOqT_EXTRA_DIST)
