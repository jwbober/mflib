CXXFLAGS = -O2 -std=c++11 -g -I/home/bober/include/flint -I$(CURDIR)/.. -Wall -Wno-sign-compare -DUSE_ARB -fdiagnostics-color=auto -Wno-unused-function -Wno-unused-variable -pthread
CFLAGS = -O2 -g -I/home/bober/include/flint -I$(CURDIR)/.. -Wall -Wno-sign-compare -DUSE_ARB -fdiagnostics-color=auto
AT=@
LDFLAGS = ../libmodform.a -lflint -larb -lgmp -lntl -lsqlite3

MISC_BINARIES = modp_basis modp_trace weight1 modp_newbasis print-traces read-polydb-entries read-polydb-entries2 read-single-mfdb-entry weight1-basis-modp weight1-hecke-matrix

all: $(MISC_BINARIES)

$(MISC_BINARIES): % : %.cc ../libmodform.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

build-factor-table: build-factor-table.cc
	g++ -O2 build-factor-factor-table.cc -std=c++11 -o build-factor-table

clean:
	-rm $(MISC_BINARIES)
	-rm build-factor-table
