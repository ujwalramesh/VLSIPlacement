CODE_ROOT=/home/rameshul/thesis/VLSIPlacement/code
OBJDIR=$(CODE_ROOT)/objdir
INCLUDE_DIR=$(CODE_ROOT)/include
INCLUDE_WNNLIB=/home/rameshul/thesis/VLSIPlacement/wnlib/acc/h
STATIC_LIB_PATH=/home/rameshul/thesis/VLSIPlacement/wnlib/acc/text.a 

CFLAGS_LOCAL = $(CFLAGS) -c

include src.mk

HFILES:=$(HFILES:=.h)
HFILES:= $(addprefix $(INCLUDE_DIR)/, $(HFILES))

OBJFILES:= $(SRCFILES:=.o)
OBJGFILES:= $(SRCFILES:=-g.o)
OBJS := $(addprefix $(OBJDIR)/,$(OBJFILES))
OBJSG := $(addprefix $(OBJDIR)/,$(OBJGFILES))
SRCCPPFILES:= $(SRCFILES:=.cxx)

$(OBJDIR)/%.o : %.cxx $(HFILES)
	g++ $(CFLAGS_LOCAL) $*.cxx -I $(INCLUDE_DIR) -I $(INCLUDE_WNNLIB) -o $(OBJDIR)/$*.o $(STATIC_LIB_PATH) -lm

$(OBJDIR)/%-g.o : %.cxx $(HFILES)
	g++ -g $(CFLAGS_LOCAL) $*.cxx -I $(INCLUDE_DIR) -I $(INCLUDE_WNNLIB) -o $(OBJDIR)/$*-g.o $(STATIC_LIB_PATH) -lm

all:  $(OBJS)
	ld -r $(OBJS) -o $(OBJDIR)/$(MODULE).o

debug:  $(OBJSG)
	ld -r $(OBJSG) -o $(OBJDIR)/$(MODULE)-g.o
