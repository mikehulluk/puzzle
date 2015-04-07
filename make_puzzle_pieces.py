
import numpy as np
from scipy import misc
import pylab
from matplotlib.pyplot import hist, figure
from scipy.ndimage.measurements import label
from matplotlib.image import imsave

from scipy.ndimage.interpolation import rotate
import random

random.seed(1)


def build_puzzle_pieces(puzzle_skeleton_fname, puzzle_img_fname):
    im = misc.imread(puzzle_skeleton_fname)

    # Turn to black and white
    im = np.mean(im, axis=2)
    im = im > 200
    pylab.imshow(im)

    

    # Build the image pieces:
    color_img = im.copy() * 0.0

    labels, numL = label(im)
    label_indices = [(labels == i).nonzero() for i in xrange(1, numL+1)]
    n_areas = len(label_indices)
    indiv_pieces = []
    for i,li in enumerate(label_indices):
        print li
        color_img[li] = 200./n_areas * i
        
        im_piece = im * 0
        im_piece[li] = 1
        
        
        indiv_pieces.append(im_piece)
    
    
    if False:    
        pylab.figure()
        pylab.imshow(color_img)
        pylab.show()
        
    
    for i, piece in enumerate(indiv_pieces):
        print "shape", piece.shape
        
        #Crop the image
        valid_ycoords = np.nonzero( piece.sum(axis=0) )[0]
        ymin, ymax = min(valid_ycoords), max(valid_ycoords)
        valid_xcoords = np.nonzero( piece.sum(axis=1) )[0]
        xmin, xmax = min(valid_xcoords), max(valid_xcoords)
        
        
        
        #new_image_size_x = (xmax - xmin) * 3
        #new_image_size_y = (xmax - xmin) * 3
        #new_image = 
        
        
        print "%d:%d, %d:%d" % (xmin,xmax, ymin, ymax) 
        
        
        padding = 3
        xmin,xmax,ymin,ymax = xmin-padding, xmax+padding, ymin-padding, ymax+padding
        xmin,xmax,ymin,ymax = max(0,xmin),min(piece.shape[0],xmax),max(0,ymin),min(ymax,piece.shape[1])
        only_piece = piece[xmin:xmax,ymin:ymax]
        
        
        # Rotate:
        piece_rot = rotate(only_piece, random.randint(0,360))
        
        imsave("build/piece%02d.png"%i, piece)    
        imsave("build/piece%02d_single.png"%i, only_piece )
        imsave("build/piece%02d_single_rot.png"%i, piece_rot)
    
    print "Len: li", len(label_indices) 
    print label_indices
    


    #figure()
    #hist(im.flatten())
    #pylab.show()
    #for (x,y),val in np.ndenumerate(im):
    #    if not handled_pixels[x,y]:
    #        pass
    #    print x,y, val
    
    
build_puzzle_pieces("src_skeletons/48ec8e583591a4f17e05349b0f42be85.png", None)
