include ../make.inc
#######################################################

RCHLCURRENTS  = new-currents/RChL-currents
OTHERCURRENTS = new-currents/other-currents

LOCAL_LIB   =  glib.a 
LIB_OBJECTS =  formf.o tauola.o curr_cleo.o pkorb.o f3pi.o \
               $(RCHLCURRENTS)/rcht_3pi/f3pi_rcht.o        $(RCHLCURRENTS)/rcht_3pi/funct_3pi.o          \
               $(RCHLCURRENTS)/rcht_common/funct_rpt.o     $(RCHLCURRENTS)/rcht_common/value_parameter.o \
               $(RCHLCURRENTS)/rcht_common/FA1RCHL.o       $(RCHLCURRENTS)/rcht_common/ffwid3pi.o        \
               $(RCHLCURRENTS)/rcht_common/initA1TabKKpi.o $(RCHLCURRENTS)/rcht_common/wid_a1_fit.o      \
               $(RCHLCURRENTS)/rcht_common/initA1Tab.o     $(RCHLCURRENTS)/rcht_common/wid_a1_fitKKpi.o  \
               $(RCHLCURRENTS)/rcht_common/gaus_integr.o \
               $(RCHLCURRENTS)/rcht_common/gfact.o         $(OTHERCURRENTS)/frho_pi_belle.o

CPP = g++

##############   What has to be built
all:	$(LOCAL_LIB) 

.c.o:
	$(CPP) $(CFLAGS) -c $< -o $@

.f.o:
	$(F77) $(FFLAGS) -c $< -o $@
$(LOCAL_LIB): $(LIB_OBJECTS)
	$(AR) rv  $(LOCAL_LIB) $(LIB_OBJECTS)
	$(RANLIB) $(LOCAL_LIB)
clean:
	rm -f *.o; rm -f *.a; rm -f *~; rm -f $(RCHLCURRENTS)/*/*.o;  rm -f $(OTHERCURRENTS)/*.o
#######################################################

