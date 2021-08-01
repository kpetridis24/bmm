#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <iostream>
#include <cstdlib>

namespace util{

    void printArray(int *arr, int len){

        for(int i = 0; i < len; i++)
            std::cout << arr[i] << std::endl;
    }


    void printMatrix(int **mat, int rows, int cols){

        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                std::cout << mat[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
}

#endif