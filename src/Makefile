# Makefile for Kuramoto simulation
# Adjust compiler and flags as needed

# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -O3 -march=native

# FFTW library flags (adjust path if needed)
FFTW_LIBS = -lfftw3 -lm

# Source files
MODULES = mod_constants.f90 mod_fft_utils.f90 mod_kuramoto.f90
MAIN = kuramoto_main.f90

# Object files
MOD_OBJS = $(MODULES:.f90=.o)
MAIN_OBJ = $(MAIN:.f90=.o)

# Module files
MOD_FILES = $(MODULES:.f90=.mod)

# Executable name
EXEC = kuramoto_simulation

# Default target
all: $(EXEC)

# Build executable
$(EXEC): $(MOD_OBJS) $(MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(FFTW_LIBS)

# Compile modules (order matters!)
mod_constants.o: mod_constants.f90
	$(FC) $(FFLAGS) -c $<

mod_fft_utils.o: mod_fft_utils.f90
	$(FC) $(FFLAGS) -c $<

mod_kuramoto.o: mod_kuramoto.f90 mod_constants.o mod_fft_utils.o
	$(FC) $(FFLAGS) -c $<

# Compile main program
$(MAIN_OBJ): $(MAIN) $(MOD_OBJS)
	$(FC) $(FFLAGS) -c $(MAIN)

# Clean up
clean:
	rm -f $(MOD_OBJS) $(MAIN_OBJ) $(MOD_FILES) $(EXEC) *.dat

# Run the simulation
run: $(EXEC)
	./$(EXEC)

.PHONY: all clean run
