#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>
#include <png++/png.hpp>
#include "stdio.h"
#include "string.h"
#include <string>
#include <sstream>
#include <chrono>
#include "mpi.h"
#include <iostream>
#include <iterator>

using namespace std;

typedef vector<double> Array;
typedef vector<Array> Matrix;
typedef vector<Matrix> Image;

Matrix lowPass(int height, int width, double sigma)
{
    Matrix kernel(height, Array(width));
    double sum=0.0;
    int i,j;

    for (i=0 ; i<height ; i++) {
        for (j=0 ; j<width ; j++) {
            kernel[i][j] = 1;
            sum += kernel[i][j];
        }
    }

    for (i=0 ; i<height ; i++) {
        for (j=0 ; j<width ; j++) {
            kernel[i][j] /= sum;
        }
    }

    return kernel;
}

Image loadImage(const char *filename)
{
    png::image<png::rgb_pixel> image(filename);
    Image imageMatrix(3, Matrix(image.get_height(), Array(image.get_width())));

    int h,w;
    for (h=0 ; h<image.get_height() ; h++) {
        for (w=0 ; w<image.get_width() ; w++) {
            imageMatrix[0][h][w] = image[h][w].red;
            imageMatrix[1][h][w] = image[h][w].green;
            imageMatrix[2][h][w] = image[h][w].blue;
        }
    }

    return imageMatrix;
}

void saveImage(Image &image, string filename)
{
    assert(image.size()==3);

    int height = image[0].size();
    int width = image[0][0].size();
    int x,y;

    png::image<png::rgb_pixel> imageFile(width, height);

    for (y=0 ; y<height ; y++) {
        for (x=0 ; x<width ; x++) {
            imageFile[y][x].red = image[0][y][x];
            imageFile[y][x].green = image[1][y][x];
            imageFile[y][x].blue = image[2][y][x];
        }
    }
    imageFile.write(filename);
}


Image applyFilter(Image &image, Matrix &filter){
    assert(image.size()==3 && filter.size()!=0);

    int height = image[0].size();
    int width = image[0][0].size();
    int filterHeight = filter.size();
    int filterWidth = filter[0].size();
    int newImageHeight = height-filterHeight+1;
    int newImageWidth = width-filterWidth+1;
    int d,i,j,h,w;
    int rank,size,rc,tag;
    MPI_Status status;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    Image newImage(3, Matrix(newImageHeight, Array(newImageWidth)));
    int slices [size];
    int numOfSlices = size -1;
    int slicesPer = newImageHeight / numOfSlices;
    int remaining = newImageHeight % numOfSlices;

    for (int r=1;r<size;r++){
        if (r <= remaining) {
            slices[r] = slicesPer+1;
        } else {
            slices[r] =slicesPer;
        }
    }



    
    if (rank == 0){
        std::cout << "newImageWidth: " << newImageWidth<< " " << std::endl;
        std::cout << "newImageHeight: " << newImageHeight<< " " << std::endl;
        std::cout << "numOfSlices: " << numOfSlices<< " " << std::endl;
        std::cout << "slicesPer: " << slicesPer<< " " << std::endl;
        std::cout << "remaining: " << remaining<< " " << std::endl;
        for (int i = 0; i<sizeof(slices)/sizeof(slices[0]); ++i)
        {
            std::cout << "slice: " <<slices[i] << std::endl;
        }
    }
    
    tag = 100;
    int finalHeight,firstHeight;
    firstHeight=0;
    finalHeight=0;
    int heightOverlap;
    if (rank==0) {
        for(int n = 1; n < size; n++){
            finalHeight += slices[n];
            if (n == size-1){
                    heightOverlap=0;
                } else {
                    heightOverlap = filterHeight;
                }
            
            for (int d = 0; d < 3; d++){
                for (int i = firstHeight; i < finalHeight + heightOverlap; i++){
                    rc = MPI_Send(&image[d][i][0], newImageWidth, MPI_DOUBLE, n, tag, MPI_COMM_WORLD);
                }
            }
            firstHeight += slices[n];
        }
    } else {
        std::cout << "rank: " << rank<< " " << std::endl;
        Image partImage(3, Matrix(slices[rank] + filterHeight+100, Array(newImageWidth)));
        for (int d = 0; d < 3; d++){
                if (rank == size-1){
                    heightOverlap=0;
                } else {
                    heightOverlap = filterHeight;
                }
            for (int i = 0; i < slices[rank] + heightOverlap; i++){
                rc = MPI_Recv(&partImage[d][i][0], newImageWidth, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
            }
        }
        Matrix filter = lowPass(10, 10, 50.0);
        Image newPartImage(3, Matrix(slices[rank], Array(newImageWidth)));
        int sum = 0;
        for (d=0;d<3;d++) {
            for (i=0 ; i< slices[rank] ; i++) {
                for (j=0 ; j<newImageWidth ; j++) {
                    sum +=partImage[d][i][j];
                }
            }
        }
        std::cout << "sum: " << sum << " " << std::endl;
        int x;
        sum=0;
        for (d=0 ; d<3 ; d++) {
            for (i=0 ; i< slices[rank] ; i++) {
                for (j=0 ; j<newImageWidth ; j++) {
                    x=0;
                    for (h=i ; h<i+filterHeight ; h++) {
                        for (w=j ; w<j+filterWidth ; w++) {
                            x+=filter[h-i][w-j]*partImage[d][h][w];
                        }
                    }
                    newPartImage[d][i][j] = x;
                    sum +=x;
                }
            }
        }
        std::cout << "sum2: " << sum << " " << std::endl;
        std::cout << "rank_finished_loop: " << rank<< " " << std::endl;
        for (int d = 0; d < 3; d++){
            for (int i = 0; i < slices[rank]; i++){
                rc = MPI_Send(&newPartImage[d][i][0], newImageWidth, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
            }
        }
    }
    if (rank == 0){
        firstHeight=0;
        finalHeight=0;
        std::cout << "rank_recv: " << rank << " " << std::endl;
        for(int n = 1; n < size; n++){
            finalHeight += slices[n];

            std::cout << "rank_send: " << n << " " << std::endl;

            for (int d = 0; d < 3; d++){
                for (int i = firstHeight; i < finalHeight; i++){
                    rc = MPI_Recv(&newImage[d][i][0], newImageWidth, MPI_DOUBLE, n, tag, MPI_COMM_WORLD,&status);
                }
            }
            firstHeight += slices[n];
            std::cout << "rank_recv: " << n << " " << std::endl;
        }
    }

    return newImage;
}

Image applyFilter(Image &image, Matrix &filter, int times)
{
    Image newImage = image;
    for(int i=0 ; i<times ; i++) {
        newImage = applyFilter(newImage, filter);
    }
    return newImage;
}

int main(int agrc, char *argv[])
{
    MPI_Init( &agrc, &argv );
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    auto t1 = std::chrono::high_resolution_clock::now();
    Matrix filter = lowPass(10, 10, 50.0);
    if (rank==0)
    cout << "Loading image..." << endl;
    
    Image image = loadImage(argv[1]);
    if (rank==0)
    cout << "Applying filter..." << endl;
    auto t1_1 = std::chrono::high_resolution_clock::now();
    
    Image newImage = applyFilter(image, filter);
    if (rank==0){
    auto t2_1 = std::chrono::high_resolution_clock::now();
    auto duration_1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2_1 - t1_1).count();
    std::cout << "Computation time: " << (float) (duration_1 / 1000.0) << " sec" << std::endl;
    
    cout << "Saving image..." << endl;

    // Generamos el nombre del fichero 
    stringstream ss;
    ss << argv[2];
    string str = ss.str();
    string ficheroGuardar = str;

    std::cout << "" << ss.str() << endl;
    saveImage(newImage, ficheroGuardar);
    cout << "Done!" << endl;
    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Execution time: " << (float) (duration / 1000.0) << " sec" << std::endl;
    }
    MPI_Finalize();
}