#
#
# Makefile
#
#
 SRCS = FBMC.cpp FBMC_MakeTables.cpp FBMC_PhononScattering.cpp FBMC_ImpactIonization.cpp EnergyBand.cpp

CC++ = g++
#C++FLAGS = -Og -Wall -Wextra -std=c++17
C++FLAGS = -O3
#CC++ = icpx
#C++FLAGS = -ipo -O3 -static -fp-model fast=2 -xHost

OBJS = $(SRCS:.cpp=.o)

all:	avalanche test

avalanche: main_avalanche.o $(OBJS)
	$(CC++) $(C++FLAGS) main_avalanche.o $(OBJS) -o avalanche.out

test:	main_test.o $(OBJS)
	$(CC++) $(C++FLAGS) $(OBJS) main_test.o -o test.out

clean	:
	rm -f *.o *.out

print	:
	cat Makefile $(SRCS)

#
main_avalanche.o: main_avalanche.cpp FBMC.h Definitions.h
	$(CC++) $(C++FLAGS) -c main_avalanche.cpp

main_test.o: main_test.cpp FBMC.h Definitions.h
	$(CC++) $(C++FLAGS) -c main_test.cpp

EnergyBand.o: EnergyBand.cpp EnergyBand.h Definitions.h
	$(CC++) $(C++FLAGS) -c EnergyBand.cpp

FBMC.o: FBMC.cpp FBMC.h Definitions.h EnergyBand.h
	$(CC++) $(C++FLAGS) -c FBMC.cpp

FBMC_MakeTables.o: FBMC_MakeTables.cpp FBMC.h Definitions.h EnergyBand.h
	$(CC++) $(C++FLAGS) -c FBMC_MakeTables.cpp

FBMC_PhononScattering.o: FBMC_PhononScattering.cpp FBMC.h Definitions.h EnergyBand.h
	$(CC++) $(C++FLAGS) -c FBMC_PhononScattering.cpp

FBMC_ImpactIonization.o: FBMC_ImpactIonization.cpp FBMC.h Definitions.h EnergyBand.h
	$(CC++) $(C++FLAGS) -c FBMC_ImpactIonization.cpp
