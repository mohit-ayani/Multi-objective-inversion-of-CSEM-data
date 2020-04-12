# GFortran Compiler (its free!):
# 
# module load swset/2018.05 intel/18.0.1 intel-mpi/2018.1.163
FC = /usr/local/bin/mpif90  
#FC = /apps/u/intel/18.0.1/intel-mpi/2018.1.163-jtnxk7c/compilers_and_libraries_2018.1.163/linux/mpi/intel64/bin/mpif90
#FC = /usr/local/bin/gfortran # use this on your machine
#Use this for debugging source code
FCFLAGS = -g -fbacktrace -fcheck=all # used to debug allocation error
#
# Intel compiler
#FC = /apps/INTEL/2013/impi/4.1.1.036/intel64/bin/mpif90 # use this on mt moran
# 
#Intel compiler for cheyenne
#FC = /glade/u/apps/ch/opt/ncarcompilers/0.4.1/mpi/mpif90 # use on cheyenne

# Use this for debugging source code:
# FCFLAGS =    -Wall  -pedantic -std=f2003
#FCFlAGS = -g -traceback -check=all # used to debug allocation error on mount moran

#----------------------------------------------------
# Build commands:
#----------------------------------------------------

TARGETS= main clean
OBJECTS= FilterModules.o Dipole1D.o CallDipole1D.o modnsga_parallel.o \
		 prognsga.o

#OBJECTS= FilterModules.o Dipole1D.o CallDipole1D.o modnsga_parallel.o prognsga.o
	 
		
#OBJSDP= FilterModules.o Dipole1D.o CallDipole1D.o
				
all:  $(TARGETS)
		
clean:	clean_msg
		rm -f *.o *~ core *.mod
		

main: build_msg_main $(OBJECTS)
		$(FC) $(FCFLAGS) -o $@ $(OBJECTS)
 
	
#		
# Compile rules
#		

# General Fortran compile:
%.o: %.f90 
	$(FC) $(FCFLAGS)    -c -o $@ $^
	
	
#	
# Build Messages:
#	
clean_msg: 
	@printf "#\n# Cleaning files: \n#\n"
	
build_msg_main: 
	@printf "#\n# Building main: \n#\n"
