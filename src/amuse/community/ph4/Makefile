
# Makefile to build the worker code.  First draft was auto-generated;
# subsequent versions are hand-coded, so save a copy.  Note that
# all compiles must use MPICXX, not CXX.

# standard amuse configuration include
# config.mk will be made after ./configure has run
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
endif
-include ${AMUSE_DIR}/config.mk

ALL	 = ph4_worker # ph4_worker_gpu
OPT	  = -I$(AMUSE_DIR)/lib/stopcond

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g -O2 $(OPT)
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -L$(AMUSE_DIR)/lib/stopcond -lstopcond -lm $(MUSE_LD_FLAGS)

CUDA_TK  ?= /usr/local/cuda

SAPPORO_LIBDIRS ?= -L$(AMUSE_DIR)/lib/sapporo_light
SAPPORO_LIBS ?= $(SAPPORO_LIBDIRS) -lsapporo

# Boost is needed by sapporo.

BOOSTLIBS ?= -L/opt/local/lib -lboost_thread-mt

CUDA_LIBDIRS ?= -L$(CUDA_TK)/lib -L$(CUDA_TK)/lib64
CUDA_LIBS ?= $(CUDA_LIBDIRS) -lcudart

LDFLAGS_GPU =	$(LDFLAGS) \
		$(SAPPORO_LIBS) \
		$(CUDA_LIBS) \
		$(BOOSTLIBS)

OBJS = interface.o

CODELIB = src/libph4.a

CODELIB_GPU = src/libph4gpu.a


AMUSE_DIR ?= ../../../..

CODE_GENERATOR ?= $(AMUSE_DIR)/build.py

all: $(ALL)

save:
	@ for f in Makefile interface.py interface.cc src/Makefile; do \
		cp -p $$f $$f.save; done

$(CODELIB):
	make -C src all CXX='$(MPICXX)'

$(CODELIB_GPU):
	make -C src gpulib CXX=$(MPICXX)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py ph4Interface -o $@

interface.h: interface.py
	$(CODE_GENERATOR) --type=h  -i amuse.support.codes.stopping_conditions.StoppingConditionInterface interface.py ph4Interface -o $@

ph4_worker: worker_code.cc interface.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) worker_code.cc $(OBJS) $(CODELIB) -o $@ $(CODELIB) $(LDFLAGS)  $(LIBS)

ph4_worker_gpu: worker_code.cc interface.h $(CODELIB_GPU) $(OBJS)
	$(MPICXX) $(CXXFLAGS) worker_code.cc $(OBJS) $(CODELIB_GPU) -o $@ $(LDFLAGS_GPU) -DGPU  $(LIBS)

%.o: %.cc
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 

wclean:
	$(RM) -r -f *.so *.dSYM *.o *.mpio *.pyc worker_code.cc interface.h 
	$(RM) -f *~ ph4_worker ph4_worker_gpu ph4_worker_sockets ph4_worker_gpu_sockets
	$(RM) -f worker_code.cc worker_code-sockets.cc

clean:
	make wclean
	make -C src clean
