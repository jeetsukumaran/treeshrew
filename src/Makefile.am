AM_CPPFLAGS = -std=c++11 @CPPFLAGS@
AM_LDFLAGS  = @LDFLAGS@

SUBDIRS = .

BUILT_SOURCES = treeshrew_info.h
EXTRA_DIST = echo_treeshrew_info.sh

treeshrew_info.h: $(HEADERS) $(SOURCES)
	$(top_srcdir)/src/echo_treeshrew_info.sh > treeshrew_info.h
