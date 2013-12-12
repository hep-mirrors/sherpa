
%.db: %
	@fl=$$(find $<); rm -f $@; \
	printf "Building '%s'(%i) " $@ $$(echo $$fl | wc -w); \
	sqlite3 $@ "create table path(file,content);"; \
	for i in $$fl; do \
	  test -d $$i && continue; printf "."; \
	  fn=$$(echo $$i | sed 's|'$<'||g;s|^/||g'); \
	  sed -e"s|'|''|g" -e "$$ s|$$|');|1" \
	    -e"1 s|^|insert into path values('$$fn','|1" \
	    $$i | sqlite3 $@; \
	done; echo " done"
