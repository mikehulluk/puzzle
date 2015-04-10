import pickle
from matplotlib.pyplot import plot, show, grid, axvline
with open("weights.pckl","r") as f:
    weights = pickle.load( f)
    
weights =  weights.tolist()
weights.sort()

sz = (8,8)
total_sides = sz[0]*sz[1] *4
flats = 2*sz[0] + 2*sz[1]
ins = (total_sides - flats) / 2
outs = ins

print flats, ins, outs

flat_start = ins
flat_end = ins + flats

plot(weights,'x')
grid()
axvline(flat_start)
axvline(flat_end)

show()
    