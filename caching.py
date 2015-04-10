
from joblib import Memory
joblib_memory = Memory(cachedir='./_cache/', verbose=0,  mmap_mode='r')