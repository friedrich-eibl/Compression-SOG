from compression.codec import Codec

import numpy as np
from PIL import Image
import imagecodecs

import comp


class CMDCCodec(Codec):
    def encode_image(self, image, out_file, **kwargs):
        value_range = kwargs["value_range"]
        quantization_level = kwargs["quantization_level"]
        return comp.encode_cpp(image, out_file, value_range, quantization_level)

    def decode_image(self, file_name):
        return comp.decode_cpp(file_name)

    def file_ending(self):
        return "bin"
