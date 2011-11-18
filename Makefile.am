ACLOCAL_AMFLAGS = -I m4

SUBDIRS       = \
		ATOOLS \
		BEAM \
		PDF \
		MODEL \
		METOOLS \
		PHASIC++ \
		MCATNLO \
		POWHEG \
		EXTRA_XS \
		AMEGIC++ \
		COMIX \
		CSSHOWER++ \
		AMISIC++ \
		AHADIC++ \
		HADRONS++ \
		PHOTONS++ \
		AddOns \
		SHERPA \
		Manual \
		Examples

EXTRA_DIST = @SVNINCLUDE@ GUIDELINES

fastinstall:
	awk 'function install(){ \
  if (old!="") system("cd "old"; make install; cd -"); \
  old=""; }BEGIN{ n=0; "pwd" | getline cur[n]; \
  cmd="$(MAKE); echo EXIT $$?"; \
while (cmd | getline) { \
  if ($$1=="EXIT") { ret=$$2; exit ret; }; \
  if ($$2=="Entering") { install(); \
    cur[n]=substr($$4,2,length($$4)-2); } \
  if ($$2=="Leaving") --n; print $$0; \
  if ((match($$0,"CXX")>0 ||  \
       match($$0,"CXXLD")>0 ||  \
       match($$0,"mode=compile")>0 ||  \
       match($$0,"mode=link")>0)) old=cur[n]; \
}}END{ close(cmd); if (ret==0) install(); }'

fi: fastinstall
