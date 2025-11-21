# executable name
TARGET = viz.exe

# Compiler
F95 = h5pfc
#F95 = ftn

# Include the current directory such that local header files take precedence
INC = -I .

# Libraries and includes
LIBS =

# Optimized compile
#FLAGS = -O2 -p -mcmodel=medium -rdynamic -fno-align-commons -fallow-argument-mismatch
FLAGS = -O2 -mcmodel=medium

# SOURces (modules have to be compiled first)
SOURCES = dtyps.f95 indata.f95 chmech.f95 arrinit.f95 HDF5_write.f95 XDMF_write.f95 post-processing.f95

# Make a list of OBJECT files
OBJECTS = $(SOURCES:.f95=.o)

# Rule for compiling f95s into objects
%.o: %.f95
	$(F95) $(INC) $(FLAGS) -c $<

# Target all - same as target
all: $(TARGET)

# Target - link the executable
$(TARGET): $(OBJECTS)
	$(F95) $(FLAGS) -o $@ $^ $(LIBS)

#! cleanup
clean:
	rm -f *.exe *.o *.mod *genmod*
#list:
#	echo $(VPATH)
