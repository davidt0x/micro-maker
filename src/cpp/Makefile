MEX_SUFFIX = mexa64

MATLAB = /home1/apps/matlab/2013a/
INTEL = /opt/apps/intel/13/composer_xe_2013.2.146/bin/intel64/

CC = icpc
OPTS = -DMATLAB_MEX_FILE -DMKL_ILP64 -fPIC -fno-omit-frame-pointer -DMX_COMPAT_32 -O3 -xHost -fopenmp -std=c++0x
TARGETS = SearchNNBMEX.${MEX_SUFFIX} Extract3DNeighborhoods.${MEX_SUFFIX} GetAllNHoodsMEX.${MEX_SUFFIX} OptimizeSolidMEX.${MEX_SUFFIX} 

INC = -I$(MATLAB)/extern/include -I$(MATLAB)/simulink/include
LIBS = -shared -Wl,--version-script,$(MATLAB)/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -fopenmp
LIBS += -Wl,-rpath-link,$(MATLAB)/bin/glnxa64 -L$(MATLAB)/bin/glnxa64 -lmx -lmex -lmat -lstdc++ -Wl,--start-group -Wl,--end-group -lpthread 

.SUFFIXES: .cpp .o .${MEX_SUFFIX}

.cpp.o:
	$(CC) $(OPTS) $(INC) -c $< -o bin/$@

.o.${MEX_SUFFIX}:
	$(CC) bin/$< $(OPTS) -o bin/$@ $(LIBS)
	@echo Binary created!!

all: ${TARGETS}

clean:
	set nonomatch; rm -f bin/*.$(MEX_SUFFIX) bin/*.o

