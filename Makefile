CXXFLAGS = -O2 -std=c++11 -g -I/home/bober/include/flint -I$(CURDIR) -Wall -Wno-sign-compare -DUSE_ARB -fdiagnostics-color=auto -Wno-unused-function -Wno-unused-variable -pthread
CXXFLAGS_NOARB = -O2 -std=c++11 -g -I/home/bober/include/flint -I$(CURDIR) -Wall -Wno-sign-compare -fdiagnostics-color=auto
CFLAGS = -O2 -g -I/home/bober/include/flint -I$(CURDIR) -Wall -Wno-sign-compare -DUSE_ARB -fdiagnostics-color=auto
AT=@

PREFIX=$(HOME)

export CXXFLAGS CFLAGS

INCLUDES = $(wildcard *.h)
#MODFORM_BINARIES = bin/modp_newbasis bin/newform-dimension bin/newforms_acb bin/hecke-polynomials bin/print-mfdb bin/export-mfdb bin/read-single-mfdb-entry bin/read-single-polydb-entry bin/read-polydb-entries bin/print-traces bin/hecke-polynomials2
MODFORM_BINARIES = bin/modp_newbasis bin/newform-dimension bin/newforms_acb bin/hecke-polynomials bin/print-mfdb bin/export-mfdb bin/read-single-mfdb-entry bin/read-polydb-entries bin/print-traces bin/hecke-polynomials2 bin/read-polydb-entries2
MF_DIRS = classnumbers arb-extras cuspforms_acb mfformat cuspforms_modp flint-extras slint

LIBMF_OBJECTS = $(foreach dir, $(MF_DIRS), $(patsubst %.c,%.o,$(wildcard $(dir)/*.c)) $(patsubst %.cc,%.o,$(wildcard $(dir)/*.cc)) )

all: libmodform.a $(MODFORM_BINARIES)

libmodform.a: $(MF_DIRS) $(INCLUDES)
	$(AT)$(foreach dir, $(MF_DIRS), MOD_DIR=$(dir); export MOD_DIR; $(MAKE) -f ../Makefile.subdirs -C $(dir) lib || exit $$?;)
	ar rcs libmodform.a $(LIBMF_OBJECTS)

clean:
	$(AT)$(foreach dir, $(MF_DIRS), MOD_DIR=$(dir); export MOD_DIR; $(MAKE) -f ../Makefile.subdirs -C $(dir) clean || exit $$?;)
	-rm $(MODFORM_BINARIES)
	-rm libmodform.a
	-rm modform-programs/*.o

$(MODFORM_BINARIES): bin/% : modform-programs/%.o libmodform.a
	$(CXX) $(CXXFLAGS) -o $@ $< libmodform.a -lflint -larb -lsqlite3 -lgmp -lntl

link:
	ln -f -s $(realpath $(wildcard *.h)) $(PREFIX)/include/
	ln -f -s $(realpath $(MODFORM_BINARIES)) $(PREFIX)/bin/
	ln -f -s $(realpath libmodform.a) $(PREFIX)/lib
