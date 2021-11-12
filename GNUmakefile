# $Id: GNUmakefile,v 1.2 2008/10/21 14:34:02 sedlak Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := musrSim
G4TARGET := $(name)
G4EXLIB := true
##LDFLAGS := $(shell root-config --glibs)
#ROOTLIBS  := ${shell root-config --glibs}
#LDFLAGS   := -L./ ${ROOTLIBS}

# Root (exlude libNew and libpthread from library list)
ROOTINC       = -I$(ROOTSYS)/include
 
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --glibs) -lMinuit -lHtml
ROOTLIBS      := $(filter-out -lNew,$(ROOTLIBS))
#ROOTLIBS      := $(filter-out -lpthread,$(ROOTLIBS))
#
#ROOTLIBS      := $(filter-out -lThread,$(ROOTLIBS))
##    ROOTLIBS := $(filter-out -lpthread,$(ROOTLIBS))
#ROOTLIBS      := $(filter-out -pthread,$(ROOTLIBS))
 
# Extra flags for G4
CPPFLAGS += $(ROOTINC)
#LDLIBS   += $(ROOTLIBS)
EXTRALIBS += $(ROOTLIBS)
CPPFLAGS += -g
ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
