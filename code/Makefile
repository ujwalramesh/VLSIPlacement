MAKE=make
# CODE_ROOT=/root/thesis/Ideas/PlacementImprovement/code
CODE_ROOT=/home/rameshul/thesis/VLSIPlacement/code
MAKEFLAGS += -j1
OBJDIR=$(CODE_ROOT)/objdir
INCLUDE_DIR=$(CODE_ROOT)/include
LIB_DIR=$(CODE_ROOT)/libs
#INCLUDE_WNLIB=$(CODE_ROOT)/include/wnlib/
STATIC_LIB_PATH=/home/rameshul/thesis/VLSIPlacement/wnlib/acc/text.a

#LIBSOOQP= -lblas /usr/local/lib/libma27.a -lgfortran -lprofiler
#LIBSOOQP= -lblas -lgfortran -lprofiler 
#changed by rameshul
LIBSOOQP= -lblas -lgfortran  
ADDINCFLAGS=/home/rameshul/Bonmin-1.7.4/build/include/coin

LIBQUADPROG= #-lQuadProgpp
LIBMOSEK= #-lmosek64

#LIBS= ./libs/text.a -L/usr/local/lib -lc -ldl -lm $(LIBMOSEK) $(LIBQUADPROG) $(LIBSOOQP)
#LIBS= -L/usr/local/lib -lc -ldl -lm $(LIBMOSEK) $(LIBQUADPROG) $(LIBSOOQP)
LIBS= -L/usr/local/lib -lc -ldl -lm $(LIBMOSEK) $(LIBQUADPROG)
MYLIBS= -L/home/rameshul/Bonmin-1.7.4/build/lib 
BOOSTLIBS=/home/rameshul/boost_1_57_0/
LD_LIBRARY_PATH= -L/home/rameshul/Bonmin-1.7.4/build/lib  
#MYBONPKG= -L/home/rameshul/Bonmin-1.7/lib -lPkg  
#MYBONPKG= `PKG_CONFIG_PATH=/home/rameshul/Bonmin-1.7/lib/pkgconfig:/usr/bin/pkg-config -libs pkg`

SRCDIRS:= Pin Cell Net Design HyperGraph common Stat PhysData \
	Flags Legalize TimingModels Bin \
	Plot ConjGrad Env Path lists FDP Grid
# Removed modules 
# CellSpread AnalyticalSolve		

OBJFILES:= $(SRCDIRS)
DEBUG = 

ifeq ($(MAKECMDGOALS), debug)
	OBJFILES:= $(SRCDIRS:=-g)
	DEBUG = debug
else
	CFLAGS= -lpthread -O3 -Wuninitialized -w -pipe -DNDEBUG -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas 
        CXXLINKFLAGS=  -wl,--rpath -Wl,/home/rameshul/Bonmin-1.7.4/build/lib
endif
OBJFILES:= $(OBJFILES:=.o)
OBJFILES:= $(addprefix $(OBJDIR)/, $(OBJFILES))



#### Below variables required to interface bonmin correctly######
LBZ2= -L/usr/lib/x86_64-linux-gnu -lbz2
LAPACK= -L/usr/lib -llapack
LIBZ= -L/usr/lib/x86_64-linux-gnu -lz
LIBBLAS= -L/usr/lib -lblas
LIBDL= -L/usr/lib32 -ldl

BONMINLIB=/home/rameshul/Bonmin-1.7.4/build/lib
BONSTLIB= $(wildcard $(BONMINLIB)/*.so)
null= 
BONLIBSUFFIX= $(null) -lm
BONLIBPREFIX= -L$(BONMINLIB)$(null) -
BONLIBNAMES= $(BONSTLIB:/home/rameshul/Bonmin-1.7.4/build/lib/lib%=l%)
BONLIBS= $(patsubst %.so,%,$(BONLIBNAMES))
BONFILES:= $(addprefix $(BONLIBPREFIX),$(BONLIBS))
#BONFILES:= $(addsuffix $(BONLIBSUFFIX),$(BONFILES))



export CFLAGS
export LD_LIBRARY_PATH 

LFLAGS= 

ifeq ($(MAKECMDGOALS), debug)
	OFILE=the_exec-g
else
	OFILE=the_exec
endif


CC=g++

TARGETS = all 

all: COMPILE_ALL GEN_EXEC

all-g: COMPILE_ALL GEN_EXEC-g

COMPILE_ALL: $(SRCDIRS)
	echo $(OBJDIR)
	echo $(BONFILES)
	@+for dir in $(SRCDIRS); do \
		make $(DEBUG) -C $$dir; \
	done	
     
GEN_EXEC: main.cxx $(OBJFILES)   
#	$(CC) -c main.cxx -o $(OBJDIR)/main.o -I $(INCLUDE_DIR)
#	$(CC) -c variables.cxx -o $(OBJDIR)/variables.o -I $(INCLUDE_DIR)
#	ld -r $(CXXLINKFLAGS) -I $(ADDINCFLAGS) -o $(OBJDIR)/main.o 
	$(CC) $(CFLAGS) $(LFLAGS) $(OBJFILES) -I $(ADDINCFLAGS) main.cxx variables.cxx -o $(OFILE) -I $(INCLUDE_DIR) $(BONFILES) $(LIBS) $(STATIC_LIB_PATH) -lm $(LBZ2) $(LAPACK) $(LIBZ) $(LIBBLAS) $(LIBDL) -I $(BOOSTLIBS) 
#	$(CC) $(CFLAGS) $(LFLAGS) $(OBJFILES)  main.cxx variables.cxx -o $(OFILE) -I $(INCLUDE_DIR) $(BONFILES) $(LIBS) $(STATIC_LIB_PATH) -lm  

GEN_EXEC-g: main.cxx $(OBJFILES) 
	$(CC) $(CFLAGS) -g $(LFLAGS) $(OBJFILES) -I $(ADDINCFLAGS) main.cxx variables.cxx -o $(OFILE) -I $(INCLUDE_DIR) $(LIBS) $(BONFILES) $(STATIC_LIB_PATH) -lm $(LBZ2) $(LAPACK) $(LIBZ) $(LIBBLAS) $(LIBDL) -I $(BOOSTLIBS)	


#$(OBJDIR)/main.o:    
#	ld -r $(MYLIBS) -I $(ADDINCFLAGS) -o $(OBJDIR)/main.o

tag: 	
	ctags . -R $(CODE_ROOT)
	sed 's/\/home\/nakul\/thesis\/VLSIPlacement\/code/\$$CODE_ROOT/g' tags > _tags
	sed 's/\/home\/mukherts\/thesis\/VLSIPlacement\/code/\$$CODE_ROOT/g' _tags > __tags
	sed 's/\/home\/rameshul/thesis\/VLSIPlacement\/code/\$$CODE_ROOT/g' _tags > __tags
	mv __tags tags
	rm -rf _tags

debug: all-g

clean:
	rm -rf $(OBJDIR)/*.o
	rm -rf $(CODE_ROOT)/the_exec
	rm -rf $(CODE_ROOT)/the_exec-g
