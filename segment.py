#!/usr/bin/env python

import os
import sys
import glob
import re

from vtk import *
from math import *

vtk.vtkObject.GlobalWarningDisplayOff() #disable warnings


# sanity check for the number of input arguments
if len(sys.argv)<1:
    print 'arg1 bmp directory'    
    raise RuntimeError('Program does not have reqd. inputs')
    
# base directory
path = sys.argv[1]
tcpath = path+ '/tracking_calibration/'
print tcpath


# pull the bmps
for bmpfile in glob.glob(os.path.join(path, '*.bmp')):
    #get corresponding tc.xform file
	
	txt=os.path.basename(bmpfile)
	re1='.*?'	# Non-greedy match on filler
	re2='(\\d+)';	# Integer Number 1
	rg = re.compile(re1+re2,re.IGNORECASE|re.DOTALL)
	m = rg.search(txt)
	if m:
		int1=m.group(1)
		#print bmpfile
		#print "("+int1+")"+"\n"
	tcfile=tcpath + int1 + 'tc.xform'
	
	
	i=0
	tcmatrix=vtk.vtkTransform()
	with open(tcfile) as f: # open the file for reading
		for line in f: # iterate over each line
			m1,m2,m3,m4 = line.split() # split it by whitespace
			tcmatrix.GetMatrix().SetElement(i,0,float(m1))
			tcmatrix.GetMatrix().SetElement(i,1,float(m2))
			tcmatrix.GetMatrix().SetElement(i,2,float(m3))
			tcmatrix.GetMatrix().SetElement(i,3,float(m4))
			i=i+1

	
	

	



	imageIn = vtk.vtkBMPReader() 
	imageIn.SetFileName(bmpfile) 
	imageIn.Update(); 


	 
		# Define viewport ranges
	  # (xmin, ymin, xmax, ymax)
	originalViewport = [0.0, 0.0, 0.5, 1.0];
	edgeViewport = [0.5, 0.0, 1.0, 1.0];
	 
	originalRenderer =vtkRenderer()
	originalRenderer.SetViewport(originalViewport);
	edgeRenderer =vtkRenderer()
	edgeRenderer.SetViewport(edgeViewport);

	renderWindow =vtkRenderWindow()
	renderWindow.SetSize(800,800);
	renderWindow.SetMultiSamples(0);
	renderWindow.AddRenderer(originalRenderer);
	renderWindow.AddRenderer(edgeRenderer);

	interactor =vtkRenderWindowInteractor()
	interactor.SetRenderWindow(renderWindow);



	imageActor =vtkImageActor()
	imageActor.SetInputData(imageIn.GetOutput());
	originalRenderer.AddActor(imageActor);




	il =vtkImageLuminance()
	il.SetInputConnection(imageIn.GetOutputPort());


#clip the roi
	bmpclipper=vtkImageClip()
	bmpclipper.SetInputConnection(il.GetOutputPort())
	bmpclipper.ClipDataOff()
	minX = 150;
	maxX = 338;
	minY = 126;
	maxY = 429; 
	minZ = 0;
	maxZ = 0;
	bmpclipper.SetOutputWholeExtent(minX,maxX,minY,maxY,minZ,maxZ)
	bmpclipper.Update()
	
	pad2 =vtkImageConstantPad()
	pad2.SetInputConnection(bmpclipper.GetOutputPort());
	pad2.SetOutputWholeExtent(0,640,0,480,0,0)

	
	


	ic =vtkImageCast()
	ic.SetOutputScalarTypeToFloat();
	ic.SetInputConnection(pad2.GetOutputPort());
	ic.Update()

	# Smooth the image
	
	
	# diffusion=vtkImageAnisotropicDiffusion2D()
	# diffusion.SetInputConnection(ic.GetOutputPort())
	# diffusion.SetDiffusionFactor(1.0)
	# diffusion.SetDiffusionThreshold(70.0)
	# diffusion.SetNumberOfIterations(5)
	
	
	
	#change deviation of smoothing #parameter
	gs =vtkImageGaussianSmooth()
	gs.SetInputConnection(ic.GetOutputPort());
	gs.SetDimensionality(2);
	gs.SetStandardDeviation(6)
	gs.SetRadiusFactors(1, 1, 0);

	# Gradient the image
	imgGradient =vtkImageGradient()
	imgGradient.SetInputConnection(gs.GetOutputPort());
	imgGradient.SetDimensionality(2);
	
	# range=imgGradient.GetOutput().GetScalarRange()
	# gradInvert=vtkImageShiftScale()
	# gradInvert.SetShift( -1.0*range[ 1 ] );
	# gradInvert.SetScale( 1.0 /( range[ 0 ] - range[ 1 ] ) );
	# gradInvert.SetOutputScalarTypeToFloat();
	# gradInvert.SetInputConnection( imgGradient.GetOutputPort() );
	# gradInvert.Update();
	
	

	imgMagnitude =vtkImageMagnitude()
	imgMagnitude.SetInputConnection(imgGradient.GetOutputPort());

	# Non maximum suppression
	nonMax =vtkImageNonMaximumSuppression()
	imgMagnitude.Update();
	nonMax.SetMagnitudeInputData(imgMagnitude.GetOutput());
	imgGradient.Update();
	nonMax.SetVectorInputData(imgGradient.GetOutput());
	#endif
	nonMax.SetDimensionality(2);

	pad =vtkImageConstantPad()
	pad.SetInputConnection(imgGradient.GetOutputPort());
	pad.SetOutputNumberOfScalarComponents(3);
	pad.SetConstant(0);

	
	i2sp1 =vtkImageToStructuredPoints()
	i2sp1.SetInputConnection(nonMax.GetOutputPort());
	#else
	pad.Update();
	i2sp1.SetVectorInputData(pad.GetOutput());
	#endif

	# Link edgles #parameter to determine threshold before linkage
	imgLink =vtkLinkEdgels()
	imgLink.SetInputConnection(i2sp1.GetOutputPort());
	imgLink.SetGradientThreshold(4);

	# Threshold links #parameter define size of links i.e largest/smallest segment length
	thresholdEdgels =vtkThreshold()
	thresholdEdgels.SetInputConnection(imgLink.GetOutputPort());
	thresholdEdgels.ThresholdBetween(8,55);
	thresholdEdgels.AllScalarsOff();

	gf =vtkGeometryFilter()
	gf.SetInputConnection(thresholdEdgels.GetOutputPort());

	i2sp =vtkImageToStructuredPoints()
	i2sp.SetInputConnection(imgMagnitude.GetOutputPort());
	#else
	pad.Update();
	i2sp.SetVectorInputData(pad.GetOutput());
	#endif

	# Subpixel them
	spe =vtkSubPixelPositionEdgels()
	spe.SetInputConnection(gf.GetOutputPort());
	#else
	i2sp.Update();
	spe.SetGradMapsData(i2sp.GetStructuredPointsOutput());
	#endif

	strip =vtkStripper()
	strip.SetInputConnection(spe.GetOutputPort());
	strip.Update()

	#get non-edge points #parameter if you have a lot of surface/backside noise. point[1] is the y coord of the points.
	k=strip.GetOutput().GetNumberOfPoints()
	realpoints = vtk.vtkPoints()
	realptids=vtkIdList()

	for i in range(0,k):
		point = [0,0,0]
		strip.GetOutput().GetPoint(i, point)
		#print point
		#print point[1]
		if point[1]<350 and point[1]>150:
				realptids.InsertUniqueId(i)
	#print realptids
			
	for i in range(0,realptids.GetNumberOfIds()):
		realpoints.InsertPoint(i,strip.GetOutput().GetPoint(realptids.GetId(i)))
		# display of picked points and cells
	#pickedPoints.SetPoints(points)
		
		
	myids=vtk.vtkIdTypeArray()
	myids.SetNumberOfComponents(1)

	for i in range(0,realptids.GetNumberOfIds()):
		myids.InsertNextValue(realptids.GetId(i))



	selectionNode=vtk.vtkSelectionNode()
	selectionNode.SetFieldType(vtk.vtkSelectionNode.POINT)
	selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES)
	selectionNode.SetSelectionList(myids);
	selectionNode.GetProperties().Set(vtk.vtkSelectionNode.CONTAINING_CELLS(), 1)

	selection=vtk.vtkSelection()
	selection.AddNode(selectionNode)

	extractSelection=vtk.vtkExtractSelection()

	extractSelection.SetInputData(0,strip.GetOutput())
	extractSelection.SetInputData(1,selection)
	extractSelection.Update()


	selected=vtk.vtkUnstructuredGrid()
	selected.ShallowCopy(extractSelection.GetOutput())

	surfaceFilter=vtk.vtkDataSetSurfaceFilter()
	surfaceFilter.SetInputData(selected)
	surfaceFilter.Update()

	mypoly=vtk.vtkPolyData()
	mypoly=surfaceFilter.GetOutput()
		
	xform_filter=vtk.vtkTransformPolyDataFilter()
	xform_filter.SetTransform(tcmatrix)
	xform_filter.SetInputData(mypoly)
	xform_filter.Update()
		

	outpath=bmpfile +".vtk"
	
	writer = vtk.vtkPolyDataWriter() 
	writer.SetFileName(outpath) 
	writer.SetInputData(xform_filter.GetOutput()) 
	writer.Update() 
	writer.Write() 
		
		
		
		
		
#uncomment below to enable interactive mode

	# dsm =vtkPolyDataMapper()
	# dsm.SetInputConnection(strip.GetOutputPort());
	# dsm.ScalarVisibilityOff();

	# planeActor =vtkActor()
	# planeActor.SetMapper(dsm);
	# planeActor.GetProperty().SetAmbient(1.0);
	# planeActor.GetProperty().SetDiffuse(0.0);

	# #Add the actors to the renderer, set the background and size
	# edgeRenderer.AddActor(planeActor);

	# #Render the image
	# interactor.Initialize();
	# renderWindow.Render();
	# renderWindow.Render();

	# interactor.Start();



