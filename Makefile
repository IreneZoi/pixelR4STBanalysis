.PHONY: all clean distclean

UNAME := $(shell uname)

OBJS = cmd_analyzer.o cmd.o command.o cmd_dtb.o error.o file.o pixel_dtb.o profiler.o protocol.o r4stest.o rpc_calls.o rpc.o rpc_error.o rpc_io.o scanner.o settings.o usb.o linux/rs232.o

ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
# root C flags =  -pthread -m64 -I/home/pitzl/ROOT/armin/root-cern/include

ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs) # with Gui

ifeq ($(UNAME), Darwin)
CXXFLAGS = -g -Os -Wall -I/usr/local/include -Wno-logical-op-parentheses -I/usr/X11/include
LDFLAGS = -lftd2xx -lreadline -L/usr/local/lib -L/usr/X11/lib -lX11
endif

ifeq ($(UNAME), Linux)

#CXXFLAGS = -g -Os -Wall -I/usr/local/include -I/usr/X11/include -pthread
CXXFLAGS = -g -pg -O2 -Wall -Wextra $(ROOTCFLAGS) -I/usr/local/include -I/usr/X11/include -pthread

#LDFLAGS = -lftd2xx -lreadline -L/usr/local/lib -L/usr/X11/lib -lX11 -pthread -lrt
LDFLAGS = -lftd2xx -lreadline -L/usr/local/lib -pthread -lrt $(ROOTGLIBS) -L/usr/X11/lib -lX11

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
	$(CXX) -o $@ $(addprefix obj/,$(OBJS)) $(LDFLAGS)

clean:
	rm -rf obj
#	rm -rf rpc_calls.cpp

distclean: clean
	rm -rf bin

################
# DEPENDENCIES #
################
-include $(addprefix obj/,$(OBJS:.o=.d))
