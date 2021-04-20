# Image processing with OMP and MPI
`
mpicxx -o parmph par-mph.cpp  -lpng
mpirun -np 2 ./parmph d.png out2.png

mpicxx -o parmph par-mph.cpp  -lpng && mpirun -np 2 ./parmph d.png out2.png
`
