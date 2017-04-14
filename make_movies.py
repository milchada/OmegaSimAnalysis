import yt
import glob
import os
import matplotlib.pylab as plt
import numpy as np 
from cv2 import VideoWriter, VideoWriter_fourcc, imread, resize

sim_sub_type = os.getcwd().split('_')[-1]
snaps = sorted(glob.glob(os.getcwd()+"/DAT/*.art"))
a_plus = [snap.split('a0')[-1] for snap in snaps][2:]
a_sp = [aplus.split('.')[-2] for aplus in a_plus]
a_s = np.asarray([float(asp)/pow(10,len(asp)) for asp in a_sp])

print a_s
field = "entropy"
snapsdone = len(glob.glob(os.getcwd()+'/yt/'+field+'/*'))

def find_center(ds):
	v, c = ds.find_max("density")
	return c

def make_images():
	for a in a_s[snapsdone:-1]: #because the last snapshot won't start with a0., instead with a1.
		name = str(a)
		if len(name) < 6:
			for i in xrange(6-len(name)):
				name += '0'
		ds = yt.load("DAT/L100_a%s.art" %  name)
		p = yt.ProjectionPlot(ds, "y", field)
		p.save()

def make_video(images, outimg=None, fps=5, size=None,
               is_color=True, format="XVID"):
    """
    Create a video from a list of images.
 
    @param      outvid      output video
    @param      images      list of images to use in the video
    @param      fps         frame per second
    @param      size        size of each frame
    @param      is_color    color
    @param      format      see http://www.fourcc.org/codecs.php
    @return                 see http://opencv-python-tutroals.readthedocs.org/en/latest/py_tutorials/py_gui/py_video_display/py_video_display.html
 
    The function relies on http://opencv-python-tutroals.readthedocs.org/en/latest/.
    By default, the video will have the size of the first image.
    It will resize every image to this size before adding them to the video.
    """
    fourcc = VideoWriter_fourcc(*format)
    vid = None
    for image in images:
        if not os.path.exists(image):
            raise FileNotFoundError(image)
        img = imread(image)
        if vid is None:
            if size is None:
                size = img.shape[1], img.shape[0]
            vid = VideoWriter(outvid, fourcc, float(fps), size, is_color)
        if size[0] != img.shape[1] and size[1] != img.shape[0]:
            img = resize(img, size)
        vid.write(img)
    vid.release()
    return vid

if __name__=='__main__':
	make_images()