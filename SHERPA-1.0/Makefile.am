MAKE = make 

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if MCATNLO_SUPPORT
MCATNLOMAKE   = HerwigTools
endif

SUBDIRS       =	\
		$(MCATNLOMAKE) \
		Tools \
		LundTools \
		Initialization \
		PerturbativePhysics \
		SoftPhysics \
		Single_Events \
		Main \
		Run
EXTRA_DIST    =	$(CVSINCLUDE)

