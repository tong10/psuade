INCLUDE = -I../../Src/Optimizers -I../../Src/DataIO -I../../Src/Util -I../../Src/Base -I../../Src/FuncApprox -I../../Src/Analysis -I../../Src/Samplings -I../../Src/Comm

all: OptBobyqa OptOUUBobyqa OptCobyla OptOUUCobyla OptLincoa OptOUULincoa OptNewuoa OptOUUNewuoa OptLBFGS OptOUULBFGS OptNomad OptOUUNomad anaGSA psuade_api

OptCobyla: OptCobyla.cpp
	g++ $(INCLUDE) -c -o OptCobyla.o OptCobyla.cpp
	g++ -o OptCobyla OptCobyla.o ../../build/lib/libpsuade.so ../../build/lib/libbobyqa.so -llapack -lblas -Wl,-rpath,../../build/lib:

OptBobyqa: OptBobyqa.cpp
	g++ $(INCLUDE) -c -o OptBobyqa.o OptBobyqa.cpp
	g++ -o OptBobyqa OptBobyqa.o ../../build/lib/libpsuade.so ../../build/lib/libbobyqa.so -llapack -lblas -Wl,-rpath,../../build/lib:

OptNomad: OptNomad.cpp
	g++ $(INCLUDE) -c -o OptNomad.o OptNomad.cpp
	g++ -o OptNomad OptNomad.o ../../build/lib/libpsuade.so ../../build/lib/libbobyqa.so -llapack -lblas -Wl,-rpath,../../build/lib:

OptNewuoa: OptNewuoa.cpp
	g++ $(INCLUDE) -c -o OptNewuoa.o OptNewuoa.cpp
	g++ -o OptNewuoa OptNewuoa.o ../../build/lib/libpsuade.so ../../build/lib/libbobyqa.so -llapack -lblas -Wl,-rpath,../../build/lib:

OptLincoa: OptLincoa.cpp
	g++ $(INCLUDE) -c -o OptLincoa.o OptLincoa.cpp
	g++ -o OptLincoa OptLincoa.o ../../build/lib/libpsuade.so ../../build/lib/libbobyqa.so -llapack -lblas -Wl,-rpath,../../build/lib:

OptLBFGS: OptLBFGS.cpp
	g++ $(INCLUDE) -c -o OptLBFGS.o OptLBFGS.cpp
	g++ -o OptLBFGS OptLBFGS.o ../../build/lib/libpsuade.so ../../build/lib/libbobyqa.so -llapack -lblas -Wl,-rpath,../../build/lib:

OptOUUCobyla: OptOUUCobyla.cpp
	g++ $(INCLUDE) -c -o OptOUUCobyla.o OptOUUCobyla.cpp
	g++ -o OptOUUCobyla OptOUUCobyla.o ../../build/lib/libpsuade.so ../../build/lib/libbobyqa.so -llapack -lblas -Wl,-rpath,../../build/lib:

OptOUUBobyqa: OptOUUBobyqa.cpp
	g++ $(INCLUDE) -c -o OptOUUBobyqa.o OptOUUBobyqa.cpp
	g++ -o OptOUUBobyqa OptOUUBobyqa.o ../../build/lib/libpsuade.so -llapack -lblas -Wl,-rpath,../../build/lib:

OptOUULincoa: OptOUULincoa.cpp
	g++ $(INCLUDE) -c -o OptOUULincoa.o OptOUULincoa.cpp
	g++ -o OptOUULincoa OptOUULincoa.o ../../build/lib/libpsuade.so ../../build/lib/liblincoa.so -llapack -lblas -Wl,-rpath,../../build/lib:

OptOUUNewuoa: OptOUUNewuoa.cpp
	g++ $(INCLUDE) -c -o OptOUUNewuoa.o OptOUUNewuoa.cpp
	g++ -o OptOUUNewuoa OptOUUNewuoa.o ../../build/lib/libpsuade.so -llapack -lblas -Wl,-rpath,../../build/lib:

OptOUULBFGS: OptOUULBFGS.cpp
	g++ $(INCLUDE) -c -o OptOUULBFGS.o OptOUULBFGS.cpp
	g++ -o OptOUULBFGS OptOUULBFGS.o ../../build/lib/libpsuade.so -llapack -lblas -Wl,-rpath,../../build/lib:

OptOUUNomad: OptOUUNomad.cpp
	g++ $(INCLUDE) -c -o OptOUUNomad.o OptOUUNomad.cpp
	g++ -o OptOUUNomad OptOUUNomad.o ../../build/lib/libpsuade.so -llapack -lblas -Wl,-rpath,../../build/lib:

anaGSA: anaGSA.cpp
	g++ $(INCLUDE) -c -o anaGSA.o anaGSA.cpp
	g++ -o anaGSA anaGSA.o ../../build/lib/libpsuade.so -llapack -lblas -Wl,-rpath,../../build/lib:

psuade_api: psuade_api.cpp
	g++ $(INCLUDE) -c -o psuade_api.o psuade_api.cpp
	g++ -o psuade_api psuade_api.o ../../build/lib/libpsuade.so -llapack -lblas -Wl,-rpath,../../build/lib:

run:
	OptBobyqa 
	OptOUUBobyqa 
	OptCobyla 
	OptOUUCobyla 
	OptLincoa 
	OptOUULincoa 
	OptNewuoa 
	OptOUUNewuoa 
	OptLBFGS 
	OptOUULBFGS 
	OptNomad
	OptOUUNomad
	psuade_api

clean:
	rm -f *.o OptBobyqa OptOUUBobyqa OptCobyla OptOUUCobyla OptLincoa OptOUULincoa OptNewuoa OptOUUNewuoa OptLBFGS OptOUULBFGS OptNomad anaGSA OptOUUNomad *.m psuade_api


