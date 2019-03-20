.PHONY: all clean distclean

UNAME := $(shell uname)

OBJS = cmd_analyzer.o cmd.o command.o cmd_dtb.o error.o file.o pixel_dtb.o profiler.o protocol.o r4stest.o rpc_calls.o rpc.o rpc_error.o rpc_io.o scanner.o settings.o usb.o linux/rs232.o

ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs) # with Gui

ifeq ($(UNAME), Darwin)
CXXFLAGS = -g -Os -Wall -Wextra $(ROOTCFLAGS) -I/usr/local/include -Wno-logical-op-parentheses -I/usr/X11/include
LDFLAGS = -lftd2xx -lreadline -L/usr/local/lib $(ROOTGLIBS) -L/usr/X11/lib -lX11
endif

ifeq ($(UNAME), Linux)

#CXXFLAGS = -g -Os -Wall -I/usr/local/include -I/usr/X11/include -pthread
CXXFLAGS = -g -pg -O2 -Wall -Wextra $(ROOTCFLAGS) -I/usr/local/include -I/usr/X11/include -I/usr/include/libftdi1 -pthread

#LDFLAGS = -lftd2xx -lreadline -L/usr/local/lib -L/usr/X11/lib -lX11 -pthread -lrt
LDFLAGS = -lftd2xx -lftdi1 -lreadline -L/usr/local/lib -pthread -lrt $(ROOTGLIBS) -L/usr/X11/lib -lX11

endif

RPCGEN = ./rpcgen/rpcgen

#################
# PATTERN RULES #
#################
obj/%.o : %.cpp
	@mkdir -p obj/linux
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/%.d : %.cpp obj
	@mkdir -p obj/linux
	$(shell $(CXX) -MM $(CXXFLAGS) $< | awk -F: '{if (NF > 1) print "obj/"$$0; else print $0}' > $@)


###########
# TARGETS #
###########
all: bin/r4stest
	@true

obj:
	@mkdir -p obj/linux

bin:
	@mkdir -p bin

rpc_calls.cpp:
	make -C rpcgen
	$(RPCGEN) pixel_dtb.h -hrpc_calls.cpp > rpcgen.log

bin/r4stest: $(addprefix obj/,$(OBJS)) bin rpc_calls.cpp
	$(CXX) -o $@ $(addprefix obj/,$(OBJS)) -fopenmp $(LDFLAGS)

r2r: r2r.cc
	g++ $(ROOTCFLAGS) r2r.cc \
	-Wall -O2 -o r2r $(ROOTLIBS)
	@echo 'done: r2r'

c2r: c2r.cc
	g++ $(ROOTCFLAGS) c2r.cc \
	-Wall -O2 -o c2r $(ROOTLIBS)
	@echo 'done: c2r'

rdroi: rdroi.cc
	g++ $(ROOTCFLAGS) rdroi.cc \
	-Wall -O2 -o rdroi $(ROOTLIBS)
	@echo 'done: rdroi'

edroi: edroi.cc
	g++ $(ROOTCFLAGS) edroi.cc \
	-Wall -O2 -o edroi $(ROOTGLIBS)
	@echo 'done: edroi'

drei: drei.cc
	g++ $(ROOTCFLAGS) drei.cc \
	-Wall -O2 -o drei $(ROOTGLIBS)
	@echo 'done: drei'

drei_reader: drei_reader.cc
	g++ $(ROOTCFLAGS) drei_reader.cc \
	-Wall -O2 -o drei_reader $(ROOTGLIBS)
	@echo 'done: drei_reader'

drei_align: drei_align.cc
	g++ $(ROOTCFLAGS) drei_align.cc \
	-Wall -O2 -o drei_align $(ROOTGLIBS)
	@echo 'done: drei_align'

zweiAB: zweiAB.cc
	g++ $(ROOTCFLAGS) zweiAB.cc \
	-Wall -O2 -o zweiAB $(ROOTGLIBS)
	@echo 'done: zweiAB'

zweiBC: zweiBC.cc
	g++ $(ROOTCFLAGS) zweiBC.cc \
	-Wall -O2 -o zweiBC $(ROOTGLIBS)
	@echo 'done: zweiBC'

drawRun: drawRun.cc
	g++ $(ROOTCFLAGS) drawRun.cc \
	-Wall -O2 -o drawRun $(ROOTLIBS)
	@echo 'done: drawRun'

ped: ped.cc
	g++ $(ROOTCFLAGS) ped.cc \
	-Wall -O2 -o ped $(ROOTLIBS)
	@echo 'done: ped'

clean:
	rm -rf obj
#	rm -rf rpc_calls.cpp

distclean: clean
	rm -rf bin

################
# DEPENDENCIES #
################
-include $(addprefix obj/,$(OBJS:.o=.d))
