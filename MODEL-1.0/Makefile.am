MAKE       = make 

if ISAJET_SUPPORT
ISAJETSUBDIR     =	Isajet
else
ISAJETSUBDIR     =
endif

if HDECAY_SUPPORT
HDECAYSUBDIR     =	Hdecay
else
HDECAYSUBDIR     =
endif

SUBDIRS    = Main Interaction_Models $(HDECAYSUBDIR) $(ISAJETSUBDIR) 

EXTRA_DIST    = @SVNINCLUDE@
