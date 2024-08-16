# Compile the .cpp file
g++ Code_scalar_activity_Langevin.cpp -lm -fopenmp -o Scalar_Activity_Langevin.o

# Run the compiled executable

#./Scalar_Activity_Langevin.o -rho 0.8 -tau 0.01 -dt 0.001 -ns 5000 -p 50 -Th 80.0 -Tc 2.0 -dumpFreq 10 -dumpFilename 80_2.lammpstrj -r 1 -restartFilename restart.jeet > log.out

./Scalar_Activity_Langevin.o -rho 0.8 -relax 0.01 -dt 0.002 -ns 5000 -p 100 -Th 5.0 -Tc 5.0 -dumpFreq 10 -dumpFilename 5_5.lammpstrj -r 0 > log.out
