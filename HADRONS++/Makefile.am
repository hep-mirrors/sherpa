SUBDIRS       =	Main ME_Library Current_Library PS_Library

DIST_SUBDIRS  = Main ME_Library Current_Library PS_Library Run

install-data-local:
	@echo "Installing Decaydata (this may take a while)."; \
	cd "$(srcdir)"; \
	for i in Decaydata/*; do \
	  if test -d "$$i"; then \
	    mkdir -p "$(DESTDIR)$(pkgdatadir)/$$i"; \
	    for j in $$i/*.dat; do \
	      cp "$$j" "$(DESTDIR)$(pkgdatadir)/$$i/"; \
	    done; \
	    if test -d "$$i/PhaseSpaceFunctions" ; then \
	      cp -rf "$$i/PhaseSpaceFunctions" "$(DESTDIR)$(pkgdatadir)/$$i/"; \
	      chmod -R +w "$(DESTDIR)$(pkgdatadir)/$$i/"; \
	      if test "x@SVNINCLUDE@" == "x"; then \
		rm -rf "$(DESTDIR)$(pkgdatadir)/$$i/PhaseSpaceFunctions/.svn"; \
	      fi; \
	    fi; \
	  fi; \
	done; \
	cp Decaydata/*.dat "$(DESTDIR)$(pkgdatadir)/Decaydata/"; \
	cd -;

uninstall-local:
	@echo "Removing Decaydata."; \
	rm -rf "$(DESTDIR)$(pkgdatadir)/Decaydata"; \
	if test "$$(ls -A $(DESTDIR)$(pkgdatadir))" = ""; then \
	rm -rf "$(DESTDIR)$(pkgdatadir)"; \
	fi;

dist-hook:
	cd $(srcdir); \
	for i in Decaydata/*; do \
	  if test -d $$i; then \
	    mkdir -p $(distdir)/$$i; \
	    for j in $$i/*.dat; do \
	      cp $$j $(distdir)/$$i/; \
	    done; \
	    if test -d $$i/PhaseSpaceFunctions; then \
	      cp -rf $$i/PhaseSpaceFunctions $(distdir)/$$i/; \
	      if test "x@SVNINCLUDE@" == "x"; then \
		rm -rf "$(distdir)/$$i/PhaseSpaceFunctions/.svn"; \
	      fi; \
	    fi; \
	  fi; \
	done; \
	cp Decaydata/*.dat $(distdir)/Decaydata/; \
	cd -;
