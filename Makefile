MODEL = driver_Transform

FC = mpiifort

#F_FLAGS = -extend-source 132 -traceback -O2
F_FLAGS = -extend-source 132 -traceback -check all -warn all -g -check noarg_temp_created 

#LINK_FLAGS = -liomp5 -lpthread -limf
LINK_FLAGS = 

OBJS = \
   Grid_Config.o \
   Grid_to_FFT.o \
   FFT_to_Legendre.o \
   Legendre_to_Spectral.o \
   Spectral_to_Legendre.o \
   Legendre_to_FFT.o \
   FFT_to_Grid.o \
   alloc_status.o \
   Transform_driver.o

.SUFFIXES: .F .f

$(MODEL): $(OBJS)
	$(FC) $(LINK_FLAGS) $(OBJS) -o $@

.F.o:
	$(FC) -c $(F_FLAGS) $<

.f.o:
	$(FC) -c $(F_FLAGS) $<

clean:
	rm -f *.o *genmod* *.mod $(MODEL)
