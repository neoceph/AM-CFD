# This a preprocessing code for thermosolutal flow code (c++) 
# Python Version

#### Import python modules 
import vtk
import math
import numpy as np
import os, sys, struct

#### Load the data from the input file 
inputfile = open(sys.argv[1], 'r')
inputdata = inputfile.readlines()

indx = [3, 12, 21]
space = [2, 4, 6]

# adjust for 0-based indexing
indx = indx - np.ones(3, dtype=int)

# grab the x grid information 
nzx = int(inputdata[indx[0]].split()[0])
xzone = np.array(map(float, inputdata[indx[0]+space[0]].split()))
ncvx = np.array(map(int, inputdata[indx[0]+space[1]].split()))
powrx = np.array(map(float, inputdata[indx[0]+space[2]].split()))

# grab the y grid information
nzy = int(inputdata[indx[1]].split()[0])
yzone = np.array(map(float, inputdata[indx[1]+space[0]].split()))
ncvy = np.array(map(int, inputdata[indx[1]+space[1]].split()))
powry = np.array(map(float, inputdata[indx[1]+space[2]].split()))

# grab the z grid information
nzz = int(inputdata[indx[2]].split()[0])
zzone = np.array(map(float, inputdata[indx[2]+space[0]].split()))
ncvz = np.array(map(int, inputdata[indx[2]+space[1]].split()))
powrz = np.array(map(float, inputdata[indx[2]+space[2]].split()))

## calculate total domain
ni = np.sum([ncvx])+2
nj = np.sum([ncvy])+2
nk = np.sum([ncvz])+2

## declare some variables (coordinates and faces)
x = np.zeros(ni)
y = np.zeros(nj)
z = np.zeros(nk)

xu = np.zeros(ni)
yv = np.zeros(nj)
zw = np.zeros(nk)

## define some varibles (boundary/internal nodes)
nim1 = ni-1
njm1 = nj-1
nkm1 = nk-1

nim2 = ni-2
njm2 = nj-2
nkm2 = nk-2

## x-grid-----------------------------------------
ist = 2
statloc = 0.0
for i in range(nzx):
    for j in range(ncvx[i]):
        if(powrx[i] >= 0.0):
            ## grids transit from fine to coarse
            term = pow((float)(j+1)/ncvx[i], powrx[i])
        else:
            ## grids transit from coarse to fine
            term = 1.0 - pow(1.0-(float)(j+1)/ncvx[i], -powrx[i])
        xu[j+ist] = statloc + xzone[i]*term
    ist += int(ncvx[i])
    statloc += xzone[i]
## coordinate value of nim1 scalar nodes can be interpolately got from that of ni uVel nodes

## central interpolation
for i in range(nim1):
    x[i] = (xu[i+1] + xu[i])*0.5

## last coordinate value of scalar node equals to that of uVel nodes
x[nim1] = xu[nim1]
 
## y-grid-----------------------------------------
ist = 2
statloc = 0.0

## loop in sub-region
for i in range(nzy):
    for j in range(ncvy[i]):
        if(powry[i] >= 0.0):
            term = pow(float(j+1) / ncvy[i], powry[i])
        else:
            term = 1.0 - pow(1.0-float(j+1) / ncvy[i], -powry[i])
        yv[j+ist] = statloc + yzone[i]*term
    ist += int(ncvy[i])
    statloc += yzone[i]

for i in range(njm1):
    y[i] = (yv[i+1] + yv[i])*0.5

## last coordinate value of scalar node equals to that of vVel nodes
y[njm1] = yv[njm1]

## z-grid-----------------------------------------
ist = 2
statloc = 0.0

## loop in sub-region
for i in range(nzz):
    for j in range(ncvz[i]):
        if(powrz[i] >= 0.0):
            term = pow(float(j+1) / ncvz[i], powrz[i])
        else:
            term = 1.0 - pow(1.0-float(j+1) / ncvz[i], -powrz[i])
        zw[j+ist] = statloc + zzone[i]*term
    ist += int(ncvz[i])
    statloc += zzone[i]

for i in range(nkm1):
    z[i] = (zw[i+1] + zw[i])*0.5

## last coordinate value of scalar node equals to that of wVel nodes
z[nkm1] = zw[nkm1]

## write out the binary paraview file
offSetCtr = 0
byteCtr = 0
npoints = nim2*njm2*nkm2

print('Total number of control volumes: {}'. format(npoints))

filename = 'preprocess_trial_mesh.vts'
fid = open(filename, 'wb')

## define internal node indecies
x1 = 1 
y1 = 1 
z1 = 1 
x2 = nim2 
y2 = njm2 
z2 = nkm2 

## general header
fid.write("<?xml version=\"1.0\"?> ")
fid.write("<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"> ")
fid.write("<StructuredGrid WholeExtent=\"%i %i %i %i %i %i\"> " % (x1, x2, y1, y2, z1, z2))
fid.write("<Piece Extent=\"%i %i %i %i %i %i\"> " % (x1, x2, y1, y2, z1, z2))

## coordinate header
fid.write("<Points> ")
fid.write("<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset = \"%i\" /> " % offSetCtr) 
offSetCtr += struct.calcsize('d') * npoints * 3 + struct.calcsize('i')
fid.write("</Points> ")

## point data header
fid.write("<PointData> ")
fid.write("</PointData> ")

## cell data header
fid.write("<CellData> ")
fid.write("</CellData> ")

## end header
fid.write("</Piece> ")
fid.write("</StructuredGrid> ")

## append header
fid.write("<AppendedData encoding=\"raw\">_")

## write out the corrdinates (in binary)
byteCtr = struct.calcsize('d') * npoints * 3;
fid.write(struct.pack('<i', byteCtr))
for k in range(1,nkm1):
    for j in range(1,njm1):
        for i in range (1,nim1):
            fid.write(struct.pack('<d', x[i]))
            fid.write(struct.pack('<d', y[j]))
            fid.write(struct.pack('<d', z[k]))


## end of appended data
fid.write("\n")
fid.write("</AppendedData> </VTKFile>")

## close the file
fid.close()

## read in the data
reader = vtk.vtkXMLStructuredGridReader()
reader.SetFileName(filename)
reader.Update()

## save off the structred grid
structGrid = reader.GetOutput()
#for i in range(530604,820000,1):
#    structGrid.BlankPoint(i)

# detele the file
#os.remove(filename)

## set up mapper for rendering the data
mapper = vtk.vtkDataSetMapper()
mapper.SetInputData(structGrid)


## set up actors
    ## main object
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().EdgeVisibilityOn()
    ## axes
transformAxes = vtk.vtkTransform()
transformGrid = vtk.vtkTransform()
transformGrid.Scale(1000*np.ones(3))
actor.SetUserTransform(transformGrid)

bounds = actor.GetBounds()
d = math.sqrt((bounds[0]-bounds[1])**2 +\
              (bounds[2]-bounds[3])**2 +\
              (bounds[4]-bounds[5])**2)

transformAxes.Scale(d*0.05*np.ones(3))
transformAxes.Translate(-d*0.03*np.ones(3))
axes = vtk.vtkAxesActor()
axes.SetUserTransform(transformAxes)

## set up the writer
writer = vtk.vtkXMLStructuredGridWriter()
writer.SetInputData(structGrid)
writer.SetFileName("trialOut.vts")
writer.SetCompressorTypeToNone()
writer.Write()

## set up the renderer
renderer = vtk.vtkRenderer()
renderer.SetBackground(1.0, 1.0, 1.0)
renderer.AddActor(actor)
renderer.AddActor(axes)
renderer.ResetCamera()

# Create the render window
render_window = vtk.vtkRenderWindow()
render_window.SetWindowName("Pre-processed Mesh")
render_window.SetSize(1800, 950)
render_window.AddRenderer(renderer)

# Create an interactor
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)

# Initialize the interactor and start the rendering loop
interactor.Initialize()
render_window.Render()
interactor.Start()
