#
#
# Makefile
#
#
 SRCS = main.cpp FBMC.cpp FBMC_MakeTables.cpp FBMC_PhononScattering.cpp FBMC_ImpactIonization.cpp EnergyBand.cpp

 CC++ = g++
# C++FLAGS = -Og -Wall -Wextra -std=c++17
 C++FLAGS = -O3

 PROGRAM = a.out
 OBJS = $(SRCS:.cpp=.o)


all	: $(OBJS)
	$(CC++) $(C++FLAGS) $(OBJS) -o $(PROGRAM)

clean	:
	rm -f *.o

print	:
	cat Makefile $(SRCS)

#
main.o	: main.cpp FBMC.h Definitions.h
	$(CC++) $(C++FLAGS) -c main.cpp

EnergyBand.o : EnergyBand.cpp EnergyBand.h Definitions.h
	$(CC++) $(C++FLAGS) -c EnergyBand.cpp

FBMC.o	: FBMC.cpp FBMC.h Definitions.h EnergyBand.h
	$(CC++) $(C++FLAGS) -c FBMC.cpp

FBMC_MakeTables.o	: FBMC_MakeTables.cpp FBMC.h Definitions.h EnergyBand.h
	$(CC++) $(C++FLAGS) -c FBMC_MakeTables.cpp

FBMC_PhononScattering.o	: FBMC_PhononScattering.cpp FBMC.h Definitions.h EnergyBand.h
	$(CC++) $(C++FLAGS) -c FBMC_PhononScattering.cpp

FBMC_ImpactIonization.o	: FBMC_ImpactIonization.cpp FBMC.h Definitions.h EnergyBand.h
	$(CC++) $(C++FLAGS) -c FBMC_ImpactIonization.cpp
