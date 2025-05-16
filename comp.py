import os
import cv2
import numpy as np
import struct
import sys

def get_frequency_table(data):
    """Calculate frequency table for the data using numpy"""
    hist, _ = np.histogram(data, bins=np.arange(np.max(data) + 2))
    return [int(count) for count in hist if count > 0]

def float_to_int(img_arr, value_range):
    return (img_arr * (value_range - 1)).astype(np.int32)

def int_to_float(img_arr, value_range):
    return (img_arr.astype(np.float32) / (value_range - 1))


def encode_cpp(img_arr, outfile, value_range, quantization_level):
    minval = np.amin(img_arr)
    maxval = np.amax(img_arr)
    zero_count = np.count_nonzero(img_arr == 0)
    print('\n Zero count before' ,  zero_count , '\n MaxVal', maxval, '\n')
    print(len(img_arr), " ", len(img_arr[0]), " ", len(img_arr[0][0]))

    save_as_binary_new(img_arr, outfile, value_range)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    encode_path = os.path.join(script_dir, "cpyth")
    command = f"{encode_path} {outfile} {quantization_level}"
    exit_code = os.system(command)
    # print(f"cpp encoded {outfile}: {minval} - {maxval}")
    print(f"quantization_level: {quantization_level}")
    return

def decode_cpp(filename):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    encode_path = os.path.join(script_dir, "decpyth")
    tmp_file = "/root/tcu/data/results_cluster/tmp_truck/compression/iteration_30000/cmdc_compressed/temp_t.bin"
    command = f"{encode_path} {filename}"
    exit_code = os.system(command)
    
    with open(tmp_file, 'rb') as f:
        rows, cols, channels = struct.unpack('iii', f.read(12))
        pixel_data = np.frombuffer(f.read(), dtype=np.int32)
        
        img_arr = pixel_data.reshape((rows, cols, channels))
        minval = np.amin(img_arr)
        maxval = np.amax(img_arr)
        zero_count = np.count_nonzero(img_arr == 0)
        print('\n Zero count after' , zero_count , '\n MaxVal', maxval, '\n')

        value_range = 50000
        if '_features_' in filename: #or 'opacity' in filename: ## no recognizable difference in size but bigger error
            value_range = 256

        ##if 'scale' in filename or 'rotation' in filename:
        ##    value_range = 256

        img_arr = int_to_float(img_arr, value_range)
        minval = np.amin(img_arr)
        maxval = np.amax(img_arr)
    
    print("Decoded image shape:", img_arr.shape, " Min:", minval, " Max:", maxval)
    return img_arr

def save_as_binary_new(img_arr, filename="temp_image.bin", value_range=50000):
    img_arr = np.ascontiguousarray(img_arr)  # Ensure contiguous memory
    rows, cols, channels = img_arr.shape
    with open(filename, 'wb') as f:
        f.write(np.array([rows, cols, channels], dtype=np.int32).tobytes())
        img_arr = (img_arr * (value_range - 1)).astype(np.int32)
        f.write(img_arr.tobytes())
