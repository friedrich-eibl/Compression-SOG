#include "image_encoding_had.hpp"
#include "encoding_library/ArithmeticCoder.hpp"
#include "encoding_library/BitIoStream.hpp"
#include "encoding_library/FrequencyTable.hpp"
#include <fstream>
#include <vector>
#include <iostream>
#include <stdexcept>

#define BLOCK_SIZE 8

using namespace std;

void arithmeticDecodeToFile(const std::string& inputFile, const std::string& outputFile) {
    // Open input file for reading
    std::ifstream inFile(inputFile, std::ios::binary);
    if (!inFile.is_open()) {
        throw std::runtime_error("Failed to open input file.");
    }

    // Read metadata (rows, cols, channels)
    int32_t rows, cols, channels;
    inFile.read(reinterpret_cast<char*>(&rows), sizeof(rows));
    inFile.read(reinterpret_cast<char*>(&cols), sizeof(cols));
    inFile.read(reinterpret_cast<char*>(&channels), sizeof(channels));

    // Validate metadata
    if (rows <= 0 || cols <= 0 || channels <= 0) {
        throw std::runtime_error("Invalid image dimensions.");
    }

    // Read number of unique symbols
    uint32_t numSymbols;
    inFile.read(reinterpret_cast<char*>(&numSymbols), sizeof(numSymbols));

    // Read symbols and frequencies
    std::vector<uint16_t> symbols(numSymbols);
    std::vector<uint32_t> freqs(numSymbols);
    for (size_t i = 0; i < numSymbols; ++i) {
        uint32_t symbol, freq;
        inFile.read(reinterpret_cast<char*>(&symbol), sizeof(uint32_t));
        inFile.read(reinterpret_cast<char*>(&freq), sizeof(uint32_t));
        symbols[i] = static_cast<uint16_t>(symbol);
        freqs[i] = freq;
    }

    // Create a frequency table for decoding
    SimpleFrequencyTable freqTable(freqs);

    // Initialize arithmetic decoder
    BitInputStream bitIn(inFile);
    ArithmeticDecoder decoder(32, bitIn);

    // Decode pixel data
    size_t totalPixels = static_cast<size_t>(rows) * cols * channels;
    std::vector<int32_t> imgData(totalPixels);
    for (size_t i = 0; i < totalPixels; ++i) {
        uint32_t idx = decoder.read(freqTable);
        imgData[i] = symbols[idx];
    }

    // Close input file
    inFile.close();

    // Write decoded data to output binary file
    std::ofstream outFile(outputFile, std::ios::binary);
    if (!outFile.is_open()) {
        throw std::runtime_error("Failed to open output file.");
    }

    // Write metadata and decoded pixel data
    outFile.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    outFile.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
    outFile.write(reinterpret_cast<const char*>(&channels), sizeof(channels));
    outFile.write(reinterpret_cast<const char*>(imgData.data()), imgData.size() * sizeof(uint16_t));
    outFile.close();

    std::cout << "Decoding complete. Output written to " << outputFile << std::endl;
}





BigPixelMatrix vectorToMatrix(const vector<int32_t>& data, int32_t rows, int32_t cols, int32_t channels) {
    BigPixelMatrix matrix(rows, vector<vector<int32_t>>(cols, vector<int32_t>(channels)));
    
    size_t index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < channels; ++k) {
                matrix[i][j][k] = data[index++];
            }
        }
    }
    return matrix;
}

BigPixelMatrix arithmeticDecodeToMatrix(const std::string& inputFile) {
    // Open input file for reading
    std::ifstream inFile(inputFile, std::ios::binary);
    if (!inFile.is_open()) {
        throw std::runtime_error("Failed to open input file.");
    }

    // Read metadata (rows, cols, channels)
    int32_t rows, cols, channels;
    inFile.read(reinterpret_cast<char*>(&rows), sizeof(rows));
    inFile.read(reinterpret_cast<char*>(&cols), sizeof(cols));
    inFile.read(reinterpret_cast<char*>(&channels), sizeof(channels));

    // Validate metadata
    if (rows <= 0 || cols <= 0 || channels <= 0) {
        throw std::runtime_error("Invalid image dimensions.");
    }

    // Read number of unique symbols
    int32_t numSymbols;
    inFile.read(reinterpret_cast<char*>(&numSymbols), sizeof(numSymbols));

    // Read symbols and frequencies
    std::vector<int32_t> symbols(numSymbols);
    std::vector<uint32_t> freqs(numSymbols);
    for (size_t i = 0; i < numSymbols; ++i) {
        uint32_t symbol, freq;
        inFile.read(reinterpret_cast<char*>(&symbol), sizeof(int32_t));
        inFile.read(reinterpret_cast<char*>(&freq), sizeof(int32_t));
        symbols[i] = static_cast<int32_t>(symbol);
        freqs[i] = freq;
    }

    // Create a frequency table for decoding
    SimpleFrequencyTable freqTable(freqs);

    // Initialize arithmetic decoder
    BitInputStream bitIn(inFile);
    ArithmeticDecoder decoder(32, bitIn);

    // Decode pixel data
    size_t totalPixels = static_cast<size_t>(rows) * cols * channels;
    std::vector<int32_t> imgData(totalPixels);
    for (size_t i = 0; i < totalPixels; ++i) {
        uint32_t idx = decoder.read(freqTable);
        imgData[i] = symbols[idx];
    }

    // Close input file
    inFile.close();
	
	
	return vectorToMatrix(imgData, rows, cols, channels);
    
}

void removeLastRowAndColumn(BigPixelMatrix& matrix) {
    if (matrix.empty() || matrix[0].empty()) {
        return;
    }
    
    // Remove last row
    matrix.pop_back();
    
    // Remove last column from each remaining row
    for (auto& row : matrix) {
        row.pop_back();
    }
}



void writeToFile(BigPixelMatrix imgData, const std::string& outputFile) {
	// Write decoded data to output binary file
	std::ofstream outFile(outputFile, std::ios::binary);
	if (!outFile.is_open()) {
		throw std::runtime_error("Failed to open output file.");
	}
	//removeLastRowAndColumn(imgData);
	int32_t rows = imgData.size();
    	int32_t cols = imgData[0].size();
    	int32_t channels = imgData[0][0].size();
    	cout << "C++ Size:" << rows << cols << endl;
    	// Create a flat buffer
    	vector<int32_t> flatData;
    	flatData.reserve(rows * cols * channels);
    
	    // Flatten the 3D matrix
	    for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
		    for (int k = 0; k < channels; ++k) {
		        flatData.push_back(static_cast<int32_t>(imgData[i][j][k]));
		    }
		}
	    }
    	
	// Write metadata and decoded pixel data
	outFile.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
	outFile.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
	outFile.write(reinterpret_cast<const char*>(&channels), sizeof(channels));
	outFile.write(reinterpret_cast<const char*>(flatData.data()), flatData.size() * sizeof(int32_t));
	outFile.close();

	std::cout << "Decoding complete. Output written to " << outputFile << std::endl;
}



BigPixelMatrix reverse_transform_matrix(const BigPixelMatrix& transformedMatrix, 
                                   const TransformationMatrix& transformationMatrix, 
                                   const TransformationMatrix& transposedTransformationMatrix) {
    int padded_height = transformedMatrix.size();
    int padded_width = transformedMatrix[0].size();
    int channels = transformedMatrix[0][0].size();
    int block_size = BLOCK_SIZE;
    int orig_height = 940;
    int orig_width = 940;
    
    // Create output matrix with original dimensions
    BigPixelMatrix pixelMatrix(orig_height, BigPixelRow(orig_width, vector<int32_t>(channels, 0)));
    
    for (int channel = 0; channel < channels; ++channel) {
        for (int row = 0; row < padded_height; row += block_size) {
            for (int col = 0; col < padded_width; col += block_size) {
                vector<vector<int32_t>> segment(block_size, vector<int32_t>(block_size));
                
                // Extract the segment for the current channel
                for (int i = 0; i < block_size; ++i) {
                    for (int j = 0; j < block_size; ++j) {
                        segment[i][j] = transformedMatrix[row + i][col + j][channel];
                    }
                }

                // Apply inverse transformation
                segment = matrix_multiply(transposedTransformationMatrix, segment);
                //test
                if (channel == 0 && row == 0 && col == 0) testVec(segment);
                segment = matrix_multiply(segment, transformationMatrix);
                //test
                if (channel == 0 && row == 0 && col == 0) testVec(segment);
                
                // Store only within original dimensions
                for (int i = 0; i < block_size; ++i) {
                    for (int j = 0; j < block_size; ++j) {
                        if (row + i < orig_height && col + j < orig_width) {
                            int32_t value = segment[i][j]/(BLOCK_SIZE*BLOCK_SIZE);360000;
                            // Clamp to uint16_t range (0 to 65535)
                            value = std::max(0, std::min(2147483647, value));
                            pixelMatrix[row + i][col + j][channel] = static_cast<int32_t>(value);
                        }
                    }
                }
            }
        }
    }
    
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
	//arithmeticDecodeToFile(src_image_path, "/root/tcu/data/results_cluster/tmp_truck/compression/iteration_30000/cmdc_compressed/temp_t.bin");
	TransformationMatrix transformationMatrix = get_transformation_matrix(8);
	TransformationMatrix transposedTransformationMatrix = get_transposed_transformation_matrix(8);
	BigPixelMatrix decodedMatrix = arithmeticDecodeToMatrix(src_image_path);
	test(decodedMatrix);
	cout << "finished decoding" <<endl;
	if(decodedMatrix.size() % BLOCK_SIZE != 0) {
		writeToFile(decodedMatrix, "/root/tcu/data/results_cluster/tmp_truck/compression/iteration_30000/cmdc_compressed/temp_t.bin");
		return 0;
	}
	BigPixelMatrix revTransMatrix = reverse_transform_matrix(decodedMatrix, transformationMatrix, transposedTransformationMatrix);
	//test(revTransMatrix);
	writeToFile(revTransMatrix, "/root/tcu/data/results_cluster/tmp_truck/compression/iteration_30000/cmdc_compressed/temp_t.bin");
	return 0;
}
