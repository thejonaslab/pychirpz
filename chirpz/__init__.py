import os
from pychirpz import chirpz, chirpz2d, zoom_fft, zoom_fft2

def get_includes():
    import eigency

    SOURCE_DIR = os.path.dirname(os.path.abspath(__file__)) 
    inc = eigency.get_includes(include_eigen=True)
    inc.append(SOURCE_DIR)
    return inc

def get_sources():
    inc.append(SOURCE_DIR)
    return inc

