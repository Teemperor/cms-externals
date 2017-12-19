include ../../../../make.inc

LIB_VER = 1.1.0

#Name of libraries to create
LIB_SPIN2_SO = libSpin2.so
LIB_SPIN2_A  = libSpin2.a

SPIN2_OBJECTS      = ME/SPIN2_init.o ME/SPIN2_distr.o ME/SPIN2_functions.o \
                     ME/cc_s2.o  ME/ccx_s2.o ME/cd_s2.o  \
                     ME/cs_s2.o ME/csx_s2.o ME/cu_s2.o \
                     ME/cux_s2.o  ME/dcx_s2.o  ME/dd_s2.o \
                     ME/ddx_s2.o ME/dsx_s2.o ME/gc_s2.o \
                     ME/gcx_s2.o ME/gg_s2.o ME/gu_s2.o \
                     ME/gux_s2.o  ME/scx_s2.o  ME/sd_s2.o \
                     ME/sdx_s2.o  ME/ss_s2.o  ME/ssx_s2.o \
                     ME/su_s2.o  ME/sux_s2.o ME/ucx_s2.o \
                     ME/ud_s2.o ME/cdx_s2.o ME/udx_s2.o \
                     ME/usx_s2.o ME/uu_s2.o \
                     ME/uux_s2.o ME/spin2distr.o

all:
	@make libSpin2


	@echo "##################################################################"	
	@echo "         library for Spin2 created in SPIN2/lib/ directory        "
	@echo "##################################################################"

%.o: %.cxx
	$(CC) $(CFLAGS) -c $< -o $@ -Iinclude
%.o: %.f
	$(F77) $(F77FLAGS) -c $< -o $@


libSpin2: $(SPIN2_OBJECTS)
	ar cr lib/$(LIB_SPIN2_A) $(SPIN2_OBJECTS)
	$(LD) $(LDFLAGS) $(SOFLAGS) $(SPIN2_OBJECTS) -o lib/$(LIB_SPIN2_SO).$(LIB_VER)
	ln -sf $(LIB_SPIN2_SO).$(LIB_VER) lib/$(LIB_SPIN2_SO)
clean:
	rm -f  lib/* ME/*.o 



