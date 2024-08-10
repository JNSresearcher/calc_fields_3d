OS=win

F90=gfortran
# F90=ifort

DO=./
DB=./
DS=./
RES=out
F=fields_3d
EXE=$(F)

ifeq ($(F90),gfortran) 
Csw = -c  -O2 -finit-local-zero  -fbounds-check -ffree-line-length-none 
Lsw = 

# Defining variables  -fbounds-check -std=gnu  
OBJ =              \
$(DO)m_vxc2data.o \
$(DO)vxc2data.o   \
$(DO)solvers.o    \
$(DO)utilites.o   \
$(DO)$(F).o 

EXE: $(OBJ)
	$(F90) $(OBJ) $(Lsw)  -o $(DB)$(EXE) 

$(OBJ): $(DO)%.o: $(DS)%.f90
	$(F90) $(Csw) $< -o $@
endif

ifeq ($(F90),ifort) 
Csw = -c /QO2  /check:bounds /Qzero  /heap-arrays:2000 /assume:buffered_io
# 
Lsw = 

OBJ =                \
$(DO)m_vxc2data.obj \
$(DO)vxc2data.obj   \
$(DO)solvers.obj    \
$(DO)utilites.obj   \
$(DO)$(F).obj 

EXE: $(OBJ)
	$(F90) $(OBJ) $(Lsw)  -o $(DB)$(EXE) 

$(OBJ): $(DO)%.obj: $(DS)%.f90
	$(F90) $(Csw) $< -o $@
endif


# Cleaning everything
clean:
ifeq ($(F90),gfortran) 
	del *.o; del *.mod; del *.exe
endif
ifeq ($(F90),ifort) 
	del *.obj; del *.mod; del *.exe
endif

cleanout:
	rmdir /s /Q out

