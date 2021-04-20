for(int n = 1; n < size; n++){
            finalHeight = newImageHeightNode * n;
            for (int j = 0; j < 3; j++){
                for (int i = firstHeight; i < finalHeight; i++){
                    rc = MPI_Send(&image[j][i][0], userWidth, MPI_DOUBLE, n, tag, MPI_COMM_WORLD);
                }
            }
            firstHeight += newImageHeightNode;
        }