Usage: fast and dirty way of making sure your tracked ultrasound captures structures you want, based off the edge detect tutorial on vtk and Tom Pheiffer's code 
This assumes you have a set of ultrasound framegrabs from an Acuson Antares system along with a set of calibrated tracked 
transformations i.e. the tracked transform is calibrated such the head of the ultrasound probe is lined up with the origin. The code does the following

1. Grab the framegrab bitmaps and crops them to obtain only the ultrasound images.
2. Segment the images using edge detection and attempt to link points to form linear structures e.g. vessels. If the images are poor quality, you can batch preprocess the images
using something like ImageJ or FIJI.I found that the difference of Gaussians filter to be very useful for detecting edges. 
3. Find the corresponding transformation file. 
4. Transform the 2D segmentations into 3D point polydata. 


needs python and vtk, written for vtk6, but can be tweaked to work with 5 or below, see vtk documentation for changes

windows cmd commands

vtkpython.exe segment.py basedirectory

basedirectory format 

<basedirectory>
		bmpfiles
		<tracking_calibration>
			#tc.xform
	

C:\Users\yife>C:\Users\yife\Desktop\VTK-6.1.0_build\bin\Release\vtkpython.exe C:
\Users\yife\Desktop\segment.py C:\Users\yife\Desktop\1210_stuff\20141210_yifei
C:\Users\yife\Desktop\1210_stuff\20141210_yifei/tracking_calibration/


if no contours seen, change gradient parameters. Last section of code can be uncommented for interactive frame by frame mode.