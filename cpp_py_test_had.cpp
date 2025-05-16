#include "image_encoding_had.hpp"
#include "encoding_library/ArithmeticCoder.hpp"
#include "encoding_library/BitIoStream.hpp"
#include "encoding_library/FrequencyTable.hpp"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>

#define BLOCK_SIZE 8

using namespace std;



std::unordered_map<int32_t, int32_t> computeFrequencyTable(const std::vector<int32_t>& data) {
    std::unordered_map<int32_t, int32_t> freqMap;
    for (const auto& val : data) {
        freqMap[val]++;
    }
    return freqMap;
}

void arithmeticEncodeFromFile(const std::string& inputFile, const std::string& outputFile) {
	time_t timestamp;
    // Open input file for reading
    std::ifstream inFile(inputFile, std::ios::binary);
    if (!inFile.is_open()) {
        throw std::runtime_error("Failed to open input file");
    }
    // Read metadata (rows, cols, channels)
    int32_t rows, cols, channels;
    inFile.read(reinterpret_cast<char*>(&rows), sizeof(rows));
    inFile.read(reinterpret_cast<char*>(&cols), sizeof(cols));
    inFile.read(reinterpret_cast<char*>(&channels), sizeof(channels));

    // Validate metadata
    if (rows <= 0 || cols <= 0 || channels <= 0) {
        throw std::runtime_error("Invalid image dimensions in the input file");
    }
	cout << "Encodin with dimensions (W, H)" << cols << rows <<endl;
    // Read pixel data
    size_t totalPixels = static_cast<size_t>(rows) * cols * channels;
    std::vector<int32_t> imgData(totalPixels);
    inFile.read(reinterpret_cast<char*>(imgData.data()), totalPixels * sizeof(int32_t));
    inFile.close();
    // Compute frequency table
    auto freqMap = computeFrequencyTable(imgData);
    // Prepare frequency table for Nayuki encoder
    std::vector<uint32_t> freqs;
    std::vector<int32_t> symbols;
    for (const auto& pair : freqMap) {
        symbols.push_back(pair.first);
        freqs.push_back(pair.second);
    }
    SimpleFrequencyTable freqTable(freqs);
    // Open output file for writing
    std::ofstream outFile(outputFile, std::ios::binary);
    if (!outFile.is_open()) {
        throw std::runtime_error("Failed to open output file");
    }

    // Write image shape to file
    outFile.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    outFile.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
    outFile.write(reinterpret_cast<const char*>(&channels), sizeof(channels));

    // Write number of unique symbols
    int32_t numSymbols = freqMap.size();
    outFile.write(reinterpret_cast<const char*>(&numSymbols), sizeof(numSymbols));

    // Write symbols and frequencies
    for (size_t i = 0; i < symbols.size(); ++i) {
        int32_t symbol = symbols[i];
        int32_t freq = freqs[i];
        outFile.write(reinterpret_cast<const char*>(&symbol), sizeof(symbol));
        outFile.write(reinterpret_cast<const char*>(&freq), sizeof(freq));
    }

    // Initialize arithmetic encoder
    BitOutputStream bitOut(outFile);
    ArithmeticEncoder encoder(32, bitOut);
    // Encode data
    // Precompute symbol index for fast lookup
std::unordered_map<int32_t, size_t> symbolIndexMap;
size_t index = 0;
for (const auto& pair : freqMap) {
    symbolIndexMap[pair.first] = index++;
}

// Encode data using precomputed indexes (O(1) lookup instead of std::find)
for (const auto& val : imgData) {
    encoder.write(freqTable, symbolIndexMap[val]);
}

    // Finish encoding
    encoder.finish();
    //bitOut.close();
    outFile.close();

    std::cout << "Encoding complete. Output written to " << outputFile << std::endl;
}


BigPixelMatrix get_pixel_matrix(string image_path) {
    // Open the file in binary mode
    std::ifstream inFile(image_path, std::ios::binary);
    if (!inFile.is_open()) {
        throw std::runtime_error("Failed to open input file");
    }

    // Read metadata (rows, cols, channels)
    int32_t rows, cols, channels;
    inFile.read(reinterpret_cast<char*>(&rows), sizeof(rows));
    inFile.read(reinterpret_cast<char*>(&cols), sizeof(cols));
    inFile.read(reinterpret_cast<char*>(&channels), sizeof(channels));

    // Validate metadata
    if (rows <= 0 || cols <= 0 || channels <= 0) {
        throw std::runtime_error("Invalid image dimensions in the input file");
    }

    // Compute total number of pixels
    size_t totalPixels = static_cast<size_t>(rows) * cols * channels;
    
    // Create a buffer to read the pixel data
    vector<int32_t> buffer(totalPixels);

    // Read pixel data into the buffer
    inFile.read(reinterpret_cast<char*>(buffer.data()), totalPixels * sizeof(int32_t));
    inFile.close();

    // Initialize the BigPixelMatrix to store the 3D image data
    BigPixelMatrix imgData(rows, vector<vector<int32_t>>(cols, vector<int32_t>(channels)));

    // Populate imgData from the flat buffer
    size_t index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < channels; ++k) {
                imgData[i][j][k] = buffer[index++];
            }
        }
    }

    return imgData;
}

void write_pixel_matrix(const BigPixelMatrix& imgData, const string& output_path) {
    // Get dimensions from the matrix
    int32_t rows = imgData.size();
    int32_t cols = imgData[0].size();
    int32_t channels = imgData[0][0].size();
    
    // Open output file for writing in binary mode
    std::ofstream outFile(output_path, std::ios::binary);
    if (!outFile.is_open()) {
        throw std::runtime_error("Failed to open output file");
    }
    
    // Write metadata (rows, cols, channels)
    outFile.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    outFile.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
    outFile.write(reinterpret_cast<const char*>(&channels), sizeof(channels));
    
    // Compute total number of pixels
    size_t totalPixels = static_cast<size_t>(rows) * cols * channels;
    
    // Create a buffer to store the flattened data
    vector<int32_t> buffer(totalPixels);
    
    // Flatten the 3D matrix into the buffer
    size_t index = 0;
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < channels; ++k) {
                buffer[index++] = imgData[i][j][k];
            }
        }
    }
    
    // Write the pixel data
    outFile.write(reinterpret_cast<const char*>(buffer.data()), totalPixels * sizeof(int32_t));
    
    outFile.close();
    
    if (!outFile) {
        throw std::runtime_error("Error occurred while writing to file");
    }
}

BigPixelMatrix transform_pixels(BigPixelMatrix pixelMatrix, TransformationMatrix transformationMatrix, TransformationMatrix transposedTransformationMatrix) {
    int orig_height = pixelMatrix.size();
    int orig_width = pixelMatrix[0].size();
    int channels = pixelMatrix[0][0].size();
    int block_size = BLOCK_SIZE;
    
    // Calculate padded dimensions
    int padded_height = ((orig_height + block_size - 1) / block_size) * block_size;
    int padded_width = ((orig_width + block_size - 1) / block_size) * block_size;
    
    cout << "Original dimensions: " << orig_height << "x" << orig_width << endl;
    cout << "Padded dimensions: " << padded_height << "x" << padded_width << endl;
    
    // Create padded matrix
    BigPixelMatrix transformed(padded_height, BigPixelRow(padded_width, vector<int32_t>(channels, 0)));
    int block_counter = 0;
    
    for (int channel = 0; channel < channels; ++channel) {
        for (int row = 0; row < padded_height; row += block_size) {
            for (int col = 0; col < padded_width; col += block_size) {
                block_counter++;
                vector<vector<int32_t>> segment(block_size, vector<int32_t>(block_size, 0)); // Initialize with zeros
                
                // Extract the segment for the current channel, using original data where available
                for (int i = 0; i < block_size; ++i) {
                    for (int j = 0; j < block_size; ++j) {
                        if (row + i < orig_height && col + j < orig_width) {
                            segment[i][j] = pixelMatrix[row + i][col + j][channel];
                        }
                        // else leave as 0 (padding)
                    }
                }
                // Apply transformation
                if (channel==0 && row == 0 && col == 0) {
                	testVec(segment);
                }
                segment = matrix_multiply(transformationMatrix, segment);
                if (channel==0 && row == 0 && col == 0) {
                	testVec(segment);
                }
                segment = matrix_multiply(segment, transposedTransformationMatrix);
                if (channel==0 && row == 0 && col == 0) {
                	testVec(segment);
                }
                
                // Store the transformed segment
                for (int i = 0; i < block_size; ++i) {
                    for (int j = 0; j < block_size; ++j) {
                        transformed[row + i][col + j][channel] = (segment[i][j]);
                    }
                }
            }
        }
    }
    
    cout << "Number of blocks created: " << block_counter << endl;
    return transformed;
}

TransformationMatrix get_quantization_matrix_luminance() {
	
	return {
		{4, 3, 4, 4, 4, 6, 11, 15},
		{3, 3, 3, 4, 5, 8, 14, 19},
		{3, 4, 4, 5, 8, 12, 16, 20},
		{4, 5, 6, 7, 12, 14, 18, 20},
		{6, 6, 9, 11, 14, 17, 21, 23},
		{9, 12, 12, 18, 23, 22, 25, 21},
		{11, 13, 15, 17, 21, 23, 25, 21},
		{13, 12, 12, 13, 16, 19, 21, 21},
	};
}

TransformationMatrix get_quantization_matrix_chrominance() {
	return {
		{12, 12, 12, 12, 12, 12, 12, 12},
		{12, 12, 12, 12, 12, 12, 12, 12},
		{12, 12, 12, 12, 12, 12, 12, 12},
		{12, 12, 12, 12, 12, 12, 12, 12},
		{12, 12, 12, 12, 12, 12, 12, 12},
		{12, 12, 12, 12, 12, 12, 12, 12},
		{12, 12, 12, 12, 12, 12, 12, 12},
		{12, 12, 12, 12, 12, 12, 12, 12},
	};
	
	return {
		{4, 4, 6, 10, 21, 21, 21, 21},
		{4, 5, 6, 21, 21, 21, 21, 21},
		{6, 6, 12, 21, 21, 21, 21, 21},
		{10, 14, 21, 21, 21, 21, 21, 21},
		{21, 21, 21, 21, 21, 21, 21, 21},
		{21, 21, 21, 21, 21, 21, 21, 21},
		{21, 21, 21, 21, 21, 21, 21, 21},
		{21, 21, 21, 21, 21, 21, 21, 21}
	};
}

BigPixelMatrix quantize_pixel_matrix (BigPixelMatrix pixelMatrix, int quantization_factor) {
	int row_c = 0;
	int col_c = 0;
	int zero_count=0;
	
	TransformationMatrix cqm = get_quantization_matrix_chrominance();
	TransformationMatrix lqm = get_quantization_matrix_luminance();
	
	if (quantization_factor == 0) {
		return pixelMatrix;
	}
	
	cout << "quantize" << endl;
	for (auto& row : pixelMatrix) {
		col_c = 0;		
        for (auto& pixel : row) {
             	
        	for (int i = 0; i < pixel.size(); i++) {
        		pixel[i] = (pixel[i] / (cqm[row_c % 8][col_c % 8]*quantization_factor))*(cqm[row_c % 8][col_c % 8]*quantization_factor);
        	}
        	col_c++;
        	
        	
        }
        row_c++;
    }
    
    cout << "zeros: " << zero_count << endl;
    return pixelMatrix;
}


void test (BigPixelMatrix pixelMatrix) {
	int rows = pixelMatrix.size();
    	int cols = pixelMatrix[0].size();
    	int ch = pixelMatrix[0][0].size();

	cout << "Rows: " << rows << ", Cols: " << cols << ", Ch: " << ch << endl;
	for (int row = 0; row < 8; row++){
	for (int col = 0; col < 8; col++){
	cout << pixelMatrix[row][col][0] << " " ;
	}cout << endl;}cout << endl;
}

int main(int argc, char *argv[]) {
	string src_image_path = argv[1];
	int quantization_level = atoi(argv[2]);
	
	BigPixelMatrix pixelMatrix = get_pixel_matrix(src_image_path);
	cout<<"start conversion"<<endl;
	TransformationMatrix transformationMatrix = get_transformation_matrix(BLOCK_SIZE);
	TransformationMatrix transposedTransformationMatrix = get_transposed_transformation_matrix(BLOCK_SIZE);
	
	if (quantization_level == 0) {
		write_pixel_matrix(pixelMatrix, src_image_path);
		cout<< "finished writing" <<endl;
		arithmeticEncodeFromFile(src_image_path, src_image_path);
		return 0;
	}
	BigPixelMatrix pixelMatrixConverted = transform_pixels(pixelMatrix, transformationMatrix, transposedTransformationMatrix);
	
	cout<<"finished conversion"<<endl;
	
	BigPixelMatrix pixelMatrixQuantized = quantize_pixel_matrix (pixelMatrixConverted, quantization_level);
	
	write_pixel_matrix(pixelMatrixQuantized, src_image_path);
	cout<< "finished writing" <<endl;
	//test(pixelMatrix);
	//test(pixelMatrixConverted);
	//test(pixelMatrixQuantized);
	//arithmeticMatrixEncode(pixelMatrix, src_image_path);
	arithmeticEncodeFromFile(src_image_path, src_image_path);
	return 0;
}



