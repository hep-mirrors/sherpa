SUBDIRS       = \
		ATOOLS \
		METOOLS \
		BEAM \
		PDF \
		MODEL \
		PHASIC++ \
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
	$(MAKE) $(MAKEFLAGS) 2>&1 | $(AWK) ' \
BEGIN{ n=0; ruc[n]=0; } \
{ \
  print $$0; \
  if ($$2=="Entering") { \
    ruc[++n]=0; \
    cur[n]=substr($$4,2,length($$4)-2); \
  } \
  if ($$2=="Leaving") { \
    if (ruc[n]==1) system("cd "cur[n]"; make install; cd -"); \
    --n; \
  } \
  if ((match($$0,"mode=compile")>0 || \
       match ($$0,"mode=link")>0) && ruc[n]==0) ruc[n]=1; \
  if (ruc[n]==1 && match($$0," Error ")>0) ruc[n]=-1; \
}END{ if (ruc[n]==1) system("make install"); }'
