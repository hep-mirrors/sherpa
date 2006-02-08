MAKE       = make 
if CVS_INCLUDE
CVSINCLUDE = .svn
endif
if ISAJET_SUPPORT
ISAJETDIR     =	Isajet
else
ISAJETDIR     =
endif
if HDECAY_SUPPORT
HDECAYDIR     =	Hdecay
else
HDECAYDIR     =
endif
SUBDIRS    = Main $(HDECAYDIR) $(ISAJETDIR) 
EXTRA_DIST    = $(CVSINCLUDE)


