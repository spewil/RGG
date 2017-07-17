# profile test

import pstats, cProfile
import pyximport
pyximport.install()
import generation 

cProfile.runctx("Sparse_RGG_Sampler()", globals(), locals(), "Profile.prof")

s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()