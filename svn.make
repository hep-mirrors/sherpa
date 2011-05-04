
MD5_EXCLUDE ?= 

SVN_Info.C: SVN_Info.C.in
	@if ! which svn > /dev/null || ! test -d $(srcdir)/.svn; then \
	  if test -f $(srcdir)/$@; then \
	    cp $(srcdir)/$@ $@.tmp; chmod u+rw $@.tmp; \
	  else \
	    echo "*********************************************"; \
	    echo "* Incomplete sources and no SVN information *"; \
	    echo "* This copy of Sherpa will not be supported *"; \
	    echo "* Please contact info@sherpa-mc.de for help *"; \
	    echo "*********************************************"; \
	    exit 1; \
	  fi; \
	else \
	  cur=$$(echo "/"$(SVNTAG) | sed -e's/[+]/[+]/g'); \
	  url=$$(svn info $(srcdir) | awk '{ if ($$1=="URL:") { \
	    split($$2,a,"svn/sherpa/"); \
	    sub("'$$cur'","",a[2]); print a[2]; } }'); \
	  rev=$$(svnversion $(srcdir)); \
	  echo '#include "ATOOLS/Org/SVN_Info.H"' > $@.tmp; \
	  echo 'static ATOOLS::SVN_Info initializer' >> $@.tmp; \
	  echo '("$(SVNTAG)","'$$url'","'$$rev'","X");' >> $@.tmp; \
	fi

	@if test -z $(NOMD5SUM); then \
	  mds=$$(cat $(addprefix $(srcdir)/, \
	    $(filter-out $@ $(CONFIG_HEADER) $(MD5_EXCLUDE), \
	    $(SOURCES) $(HEADERS))) | $(MD5COMMAND)); \
	  $(SEDCOMMAND) -e's/".?"\);/"'$$mds'");/g' $@.tmp; \
	fi;

	@if ! diff $@.tmp $@ > /dev/null 2>&1; then \
	  mv $@.tmp $@; \
	else \
	  rm $@.tmp; \
	fi;

.PHONY: SVN_Info.C.in

DISTCLEANFILES = SVN_Info.C
