//#ifndef IMG_ENC_HPP
#define IMG_ENC_HPP

#include <opencv2/opencv.hpp>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;
using PixelRow = vector<vector<uint16_t>>;
using PixelMatrix = vector<PixelRow>;

using BigPixelRow = vector<vector<int32_t>>;;
using BigPixelMatrix = vector<BigPixelRow>;

using TransformationMatrix = vector<vector<int32_t>>;


//#endif
void testVec (vector<vector<int32_t>> seg) {
	int rows = seg.size();
    	int cols = seg[0].size();
	cout << "testvec:"<<endl;
	cout << "Rows: " << rows << ", Cols: " << cols  << endl;
	for (int row = 0; row < 8; row++){
	for (int col = 0; col < 8; col++){
	cout << seg[row][col] << " " ;
	}cout << endl;}cout << endl;
}

vector<vector<int32_t>> matrix_multiply(const vector<vector<int32_t>>& matrix1, const vector<vector<int32_t>>& matrix2) {
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int cols2 = matrix2[0].size();
    
    vector<vector<int32_t>> result(rows1, vector<int32_t>(cols2, 0));
    
    for (int i = 0; i < rows1; i++) {
        for (int j = 0; j < cols2; j++) {
            int32_t sum = 0;
            for (int k = 0; k < cols1; k++) {
                sum += matrix1[i][k] * matrix2[k][j];
            }
            result[i][j] = static_cast<int32_t>(sum);
        }
    }
    
    return result;
}

TransformationMatrix get_transposed_transformation_matrix(uint8_t matrix_size) {
	if (matrix_size == 4) {
		return {
			{1, 1, 1, 1},
    		{1, -1, 1, -1},
    		{1, 1, -1, -1},
    		{1, -1, -1, 1}
		};
	}
	
	if (matrix_size == 8) {/*
	return {
    {  8, 11, 10,  9,  8,  7,  5,  3 },
    {  8, 10,  7,  5,  0, -5, -7, -11 },
    {  8,  9,  3,  0, -8, -11, -3,  7 },
    {  8,  7,  0, -3, -8, -9, 10,  5 },
    {  8,  5, -3, -7,  8,  0, 10, -10 },
    {  8,  3, -7, -9,  8,  9, -3,  7 },
    {  8,  1, -10, -11, -8, 10, -7,  0 },
    {  8,  0, -11, -10, -8,  3,  5, -9 }
};
		return {
		    {4, 6, 5, 4, 4, 3, 2, 1},
		    {4, 5, 2, -1, -4, -6, -5, -3},
		    {4, 3, -2, -6, -4, 1, 5, 5},
		    {4, 1, -5, -3, 4, 5, -2, -6},
		    {4, -1, -5, 3, 4, -5, -2, 6},
		    {4, -3, -2, 6, -4, -1, 5, -5},
		    {4, -5, 2, 1, -4, 6, -5, 3},
		    {4, -6, 5, -4, 4, -3, 2, -1}

		};
		return {
			{64,  89,  83,  75,  64,  50,  36,  18},
		    	{64,  75,  36, -18, -64, -89, -83, -50},
		    	{64,  50, -36, -89, -64, -18,  83,  75},
		    	{64,  18, -83, -50,  64,  75, -36, -89},
		    	{64, -18, -83,  50,  64, -75, -36,  89},
		    	{64, -50, -36,  89, -64, -18,  83, -75},
		    	{64, -75,  36,  18, -64,  89, -83,  50},
		    	{64, -89,  83, -75,  64, -50,  36, -18}
		};
		return {
			{32,  48,  48,  48,  40,  56,  56,  56},
		    {48, -53,  47, -78,  60, -41,  59, -66},
		    {48,  47, -53, -78,  60,  59, -41, -66},
		    {48, -78, -78,  72,  60, -66, -66,  84},
		    {40,  60,  60,  60, -100, -80, -80, -80},
		    {56, -41,  59, -66, -80, 123, -77,  98},
		    {56,  59, -41, -66, -80, -77, 123,  98},
		    {56, -66, -66,  84, -80,  98,  98, -52}
		   };*/
		return {
			{1, 1, 1, 1, 1, 1, 1, 1},
			{1, -1, 1, -1, 1, -1, 1, -1},
			{1, 1, -1, -1, 1, 1, -1, -1},
			{1, -1, -1, 1, 1, -1, -1, 1},
			{1, 1, 1, 1, -1, -1, -1, -1},
			{1, -1, 1, -1, -1, 1, -1, 1},
			{1, 1, -1, -1, -1, -1, 1, 1},
			{1, -1, -1, 1, -1, 1, 1, -1}
		};		
	}
	return {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
}
TransformationMatrix get_transformation_matrix(uint8_t matrix_size) {
	if (matrix_size == 4) {
		return {
			{1, 1, 1, 1},
    		{1, -1, 1, -1},
    		{1, 1, -1, -1},
    		{1, -1, -1, 1}
		};
	}
	if (matrix_size == 8) {/*
	return {
    {  8,  8,  8,  8,  8,  8,  8,  8 },
    { 11, 10,  9,  7,  5,  3,  1,  0 },
    { 10,  7,  3,  0, -3, -7, -10, -11 },
    {  9,  5,  0, -3, -7, -9, -11, -10 },
    {  8,  0, -8, -8,  8,  8, -8, -8 },
    {  7, -5, -11, -9,  0,  9,  10,  3 },
    {  5, -7, -3, 10, 10, -3, -7,  5 },
    {  3, -11,  7,  5, -10,  7,  0, -9 }
};
		return {
		    {4, 4, 4, 4, 4, 4, 4, 4},
		    {6, 5, 3, 1, -1, -3, -5, -6},
		    {5, 2, -2, -5, -5, -2, 2, 5},
		    {4, -1, -6, -3, 3, 6, 1, -4},
		    {4, -4, -4, 4, 4, -4, -4, 4},
		    {3, -6, 1, 5, -5, -1, 6, -3},
		    {2, -5, 5, -2, -2, 5, -5, 2},
		    {1, -3, 5, -6, 6, -5, 3, -1}

		};
		return {
			{64,  64,  64,  64,  64,  64,  64,  64},
	 		{89,  75,  50,  18, -18, -50, -75, -89},
	    		{83,  36, -36, -83, -83, -36,  36,  83},
	    		{75, -18, -89, -50,  50,  89,  18, -75},
		    	{64, -64, -64,  64,  64, -64, -64,  64},
		    	{50, -89, -18,  75, -75,  18,  89, -50},
		    	{36, -83,  83, -36, -36,  83, -83,  36},
		    	{18, -50,  75, -89,  89, -75,  50, -18}
		};
		return {
			{2,  2,  2,  2,  2,  1,  1,  1},
		    	{2, -2,  2, -1,  1, -1,  1, -1},
		    	{2,  2, -2, -1,  1,  1, -1, -1},
		    	{2, -1, -1,  1,  1, -1, -1,  1},
		    	{2,  1,  1,  1, -1, -1, -1, -1},
		    	{1, -1,  1, -1, -1,  1, -1,  1},
		    	{1,  1, -1, -1, -1, -1,  1,  1},
		    	{1, -1, -1,  1, -1,  1,  1, -1}
		};*/
		return {
			{1, 1, 1, 1, 1, 1, 1, 1},
			{1, -1, 1, -1, 1, -1, 1, -1},
			{1, 1, -1, -1, 1, 1, -1, -1},
			{1, -1, -1, 1, 1, -1, -1, 1},
			{1, 1, 1, 1, -1, -1, -1, -1},
			{1, -1, 1, -1, -1, 1, -1, 1},
			{1, 1, -1, -1, -1, -1, 1, 1},
			{1, -1, -1, 1, -1, 1, 1, -1}
		};		
	}
	return {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
}
