# Image processing with OMP and MPI
To execute the MPI code:
```
mpicxx -o parmph.out par-mph.cpp  -lpng
mpirun -np 2 ./parmph d.png out2.png

mpicxx -o parmph.out par-mph.cpp  -lpng && mpirun -np 2 ./parmph.out d.png out.png
```

To compile the OpenMP code:
```
g++ -o cp2.out custom-parallelv2.cpp -lpng -fopenmp
```
Execute:
```
./cp2.out d.png out.png
```