# define the C++ compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -Wall -g 
#-std=c++11

# define any directories containing header files other than /usr/include
INCLUDES = -I.

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = 

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = 

# define the C++ source files
SRCS = UtilityRun.cpp Hungarian.cpp AuxMathGa.cpp MarquesEnantiomers.cpp EnantiomerIdentification.cpp IsomersToMol.cpp IdentifyIsomers.cpp Geometries.cpp ChangeNames.cpp AllMolecularFormulas.cpp CauchyIndex.cpp Combinations.cpp AdjustSaParameters.cpp AuxMath.cpp BestPermutation.cpp BuildComplex.cpp ComplexCreator.cpp ControlMopac.cpp FindIsomers.cpp Fitness.cpp KabschRmsd.cpp Ligand.cpp main.cpp MyExceptions.cpp PointAnalysis.cpp ReadInput.cpp ReadMopac.cpp ReadQuantumOutput.cpp RootMeanSquareDeviation.cpp WriteQuantumInput.cpp 


# define the C++ object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.cpp=.o)

# define the executable file 
MAIN = lumpacview.exe

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)
	@echo  lumpacview.exe foi compilado com sucesso:

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
