EXEPATH = ..
COMPILER = icpc

CPPFLAGS = -Ofast -fp-model fast -ipo -std=c++0x 
IFLAGS = 
LFLAGS = 

OBJS_MD = Table.o Potentials.o System.o Handlers.o Integrators.o NSearch.o Force.o Simulation.o IO.o SimulationSetup.o main.o
OBJS_CG = System.o IO.o CG.o
OBJS_dMD = Table.o Potentials.o System.o Handlers.o Integrators.o NSearch.o Force.o Simulation.o IO.o SimulationSetup.o main_direct.o
OBJS_Conv = Table.o Potentials.o Convert_alloy_tab.o

EXE_MD = $(EXEPATH)/mse627
EXE_CG = $(EXEPATH)/CG
EXE_dMD = $(EXEPATH)/d_mse627
EXE_Conv = $(EXEPATH)/Convert_alloy_tab

# compile and load
default:
	@echo " "
	@echo "Compiling Code MSE6270_MD"
	@echo "C++ version"
	$(MAKE) $(EXE_MD)
	$(MAKE) $(EXE_CG)
#	$(MAKE) $(EXE_dMD)
#	$(MAKE) $(EXE_Conv)

$(EXE_MD): $(OBJS_MD)
	$(COMPILER) $(CPPFLAGS) $(IFLAGS) -o $(EXE_MD) $(OBJS_MD) $(LFLAGS)
	
$(EXE_CG): $(OBJS_CG)
	$(COMPILER) $(CPPFLAGS) $(IFLAGS) -o $(EXE_CG) $(OBJS_CG) $(LFLAGS)
	
$(EXE_dMD): $(OBJS_dMD)
	$(COMPILER) $(CPPFLAGS) $(IFLAGS) -o $(EXE_dMD) $(OBJS_dMD) $(LFLAGS)
	
$(EXE_Conv): $(OBJS_Conv)
	$(COMPILER) $(CPPFLAGS) $(IFLAGS) -o $(EXE_Conv) $(OBJS_Conv) $(LFLAGS)

%.o: %.cpp 
	$(COMPILER) $(CPPFLAGS) $(IFLAGS) -c $< -o $@
	
clean:
	rm -f *.o

