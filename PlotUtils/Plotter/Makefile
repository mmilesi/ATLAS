#DEST	=   RealFakeEffFitter
#OBJS1	=   RealFakeEffFitter.o
DEST	=   RealFakeEffFitter_Class
OBJS1	=   RealFakeEffFitter_Class.o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

# set the compiler options
CXX           = g++
CXXFLAGS      = -g -Wall -Wextra -O2
                                   #-g compile with debug information
                                   #-O optimize (-O2 is higher grade of optimalization)
                                   #-Wall and -Wextra show warnings for everything

# set the linker options
LD            = g++
LDFLAGS       = -O2

#######################
CXXFLAGS   = $(ROOTCFLAGS)
LIBS       = $(ROOTLIBS) -L $(ROOTSYS)/lib -lMinuit
GLIBS      = $(ROOTGLIBS)

all: $(DEST)

$(DEST): $(OBJS1)
	@$(LD) $(CXXFLAGS) $(LIBS) -o $(DEST) $(OBJS1)

clean:
	@rm -rf $(DEST) $(OBJS1)
