# *************************************************************
# Unaxial Tensile Test
# *************************************************************
from abaqus import *
from abaqusConstants import *
import regionToolset
import os
import numpy as np
import sys
import json

# Change work direction
simulation_name = sys.argv[-1]
working_directory = os.path.join(os.getcwd(),simulation_name)
os.chdir(working_directory)

# Read the initial data
data = json.load(open('data.json'))

print('--------------GENERATE NEW TENSILE TEST MODEL--------------')
# Create a new model data base
Mdb()
session.viewports['Viewport: 1'].setValues(displayedObject=None)

# ---------------------------------------------------------------
# create the model
mdb.models.changeKey(fromName='Model-1', toName='Specimen')
specimenModel = mdb.models['Specimen']

# ----------------------------------------------------------------
# Create the part
import sketch
import part

# a) Sketch the beam cross section using rectangle tool
width = 25
mid_width = 12.5
length = 180
gauge_length = 60
r_round = 25
# a = 43.4641
e = (width - mid_width) / 2
f = (length - gauge_length)/2
''' a is the solution of equation system: (x-x_0)^2 + (y-y_0)^2 =  R^2 and y = width'''
a = -np.sqrt(r_round**2 - (width - ((width - e) + r_round ))**2) + f
b = f - a
c = 5
d = gauge_length - c * 2

# Save geometry dimensions
data['specimen_length'] = length
data['gauge_length'] = gauge_length
json.dump(data, open('data.json', 'w'))

SpecimenProfileSketch = specimenModel.ConstrainedSketch(name='Specimen Profile', sheetSize=length)
SpecimenProfileSketch.Line(point1=(0, 0), point2=(a, 0))
SpecimenProfileSketch.Line(point1=(a, 0), point2=(a, e))
SpecimenProfileSketch.Line(point1=(a, e), point2=(a + b , e))
SpecimenProfileSketch.Line(point1=(a + b , e), point2=(a + b + c, e))
SpecimenProfileSketch.Line(point1=(a + b + c , e), point2=(a + b + c + d, e))
SpecimenProfileSketch.Line(point1=(a + b + c + d, e), point2=(a + b + 2 * c + d, e))
SpecimenProfileSketch.Line(point1=(a + b + 2 * c + d, e), point2=(a + 2 * b + 2 * c + d, e))
SpecimenProfileSketch.Line(point1=(a + 2 * b + 2 * c + d, e), point2=(a + 2 * b + 2 * c + d, 0))
SpecimenProfileSketch.Line(point1=(a + 2 * b + 2 * c + d, 0), point2=(length, 0))
SpecimenProfileSketch.Line(point1=(length, 0), point2=(length, width))
SpecimenProfileSketch.Line(point1=(length, width), point2=(length - a, width))
SpecimenProfileSketch.Line(point1=(length - a, width), point2=(length - a, width - e))
SpecimenProfileSketch.Line(point1=(length - a, width - e), point2=(length - a - b, width - e))
SpecimenProfileSketch.Line(point1=(length - a - b , width - e), point2=(length - a - b - c, width - e))
SpecimenProfileSketch.Line(point1=(length - a - b - c , width - e), point2=(length - a - b - c - d, width - e))
SpecimenProfileSketch.Line(point1=(length - a - b - c - d, width - e), point2=(length - a - b - 2 * c - d, width - e))
SpecimenProfileSketch.Line(point1=(length - a - b - 2 * c - d, width - e), point2=(length - a - 2 * b - 2 * c - d, width - e))
SpecimenProfileSketch.Line(point1=(length - a - 2 * b - 2 * c - d, width - e), point2=(length - a - 2 * b - 2 * c - d, width))
SpecimenProfileSketch.Line(point1=(length - a - 2 * b - 2 * c - d, width), point2=(0, width))
SpecimenProfileSketch.Line(point1=(0, width), point2=(0, 0))

# b) Create a 3D deformable part named "Beam" by extruding the sketch !
thickness = data['thickness']
specimenPart = specimenModel.Part(name='Specimen', dimensionality=THREE_D, type=DEFORMABLE_BODY)
specimenPart.BaseSolidExtrude(sketch=SpecimenProfileSketch, depth=thickness)
edge1 = specimenPart.edges.findAt((a, e,0.5*thickness,))
specimenPart.Round(radius=r_round, edgeList=(edge1,))
edge2 = specimenPart.edges.findAt((length - a, e,0.5*thickness,))
specimenPart.Round(radius=r_round, edgeList=(edge2,))
edge3 = specimenPart.edges.findAt((length - a, width - e,0.5*thickness,))
specimenPart.Round(radius=r_round, edgeList=(edge3,))
edge4 = specimenPart.edges.findAt((a, width - e,0.5*thickness,))
specimenPart.Round(radius=r_round, edgeList=(edge4,))

# Partition the specimen
# 1
selected_cell = specimenPart.cells.findAt((length/2, width/2,0.5*thickness,))
point1 = (a,e,0)
point2 = (a,e,thickness)
point3 = (a,width-e,0)
specimenPart.PartitionCellByPlaneThreePoints(point1=point1, point2=point2, point3=point3,
        cells=selected_cell)

# 2 check
selected_cell = specimenPart.cells.findAt((length/2, width/2,0.5*thickness,))
point1 = (a + b,e,0)
point2 = (a + b,e,thickness)
point3 = (a + b,width-e,0)
specimenPart.PartitionCellByPlaneThreePoints(point1=point1, point2=point2, point3=point3,
        cells=selected_cell)

# 3
selected_cell = specimenPart.cells.findAt((length/2, width/2,0.5*thickness,))
point1 = (a+b+c,e,0)
point2 = (a+b+c,e,thickness)
point3 = (a+b+c,width-e,0)
specimenPart.PartitionCellByPlaneThreePoints(point1=point1, point2=point2, point3=point3,
        cells=selected_cell)

# 4
selected_cell = specimenPart.cells.findAt((length/2, width/2,0.5*thickness,))
point1 = (a+2*b+2*c + d,e,0)
point2 = (a+2*b+2*c + d,e,thickness)
point3 = (a+2*b+2*c + d,width-e,0)
specimenPart.PartitionCellByPlaneThreePoints(point1=point1, point2=point2, point3=point3,
        cells=selected_cell)

# 5
selected_cell = specimenPart.cells.findAt((length/2, width/2,0.5*thickness,))
point1 = (a+b+2*c + d,e,0)
point2 = (a+b+2*c + d,e,thickness)
point3 = (a+b+2*c + d,width-e,0)
specimenPart.PartitionCellByPlaneThreePoints(point1=point1, point2=point2, point3=point3,
        cells=selected_cell)

# 6
selected_cell = specimenPart.cells.findAt((length/2, width/2,0.5*thickness,))
point1 = (a+b+c+d,e,0)
point2 = (a+b+c+d,e,thickness)
point3 = (a+b+c+d,width-e,0)
specimenPart.PartitionCellByPlaneThreePoints(point1=point1, point2=point2, point3=point3,
        cells=selected_cell)
# # ----------------------------------------------------------------
# Create material
import material

# Create material AISI tees Steel by assigning mass density, youngs
# modulus and poissons ratio
specimenMaterial = specimenModel.Material(name='AISI 1005 Steel')
density = data['density']
young_modulus = data['young_modulus']
nuy = data['nuy']
plastic_table = []
strain = 0
s0, K, t, h = data['s0'],data['K'],data['t'],data['h']
while strain <= 1:
    strain = round(strain, 4)
    sigma = s0 + K * (1 - np.exp(-t * strain)) * (strain + 0.002) ** h
    sigma = round(sigma, 4)
    plastic_table.append((sigma, strain))
    strain += 0.005
plastic_table = tuple(plastic_table)
density_table = ((density,),)
elastic_table = ((young_modulus, nuy),)
specimenMaterial.Density(table=density_table)
specimenMaterial.Elastic(table=elastic_table)
specimenMaterial.Plastic(table=plastic_table)

# ----------------------------------------------------------------
# Create solid section and assign the beam to it
import section

# create a section to assign to the beam
specimenSection = specimenModel.HomogeneousSolidSection(name='Specimen Section', material='AISI 1005 Steel')
# Assign the beam to this section
specimen_region = (specimenPart.cells,)
specimenPart.SectionAssignment(region=specimen_region, sectionName='Specimen Section')

# ----------------------------------------------------------------
# Create the assembly
import assembly

# Create the part instance
specimenAssembly = specimenModel.rootAssembly
specimenInstance = specimenAssembly.Instance(name='Specimen Instance', part=specimenPart, dependent=ON)

# ------------------------------------------------------------------------
# Create the step
import step
# Create a static general step
displacement_step = specimenModel.StaticStep(name='Apply Displacement Load', previous='Initial', description='Displacement load is applied during this step')
displacement_step.setValues(initialInc=0.1)
# ------------------------------------------------------------------------
# Create the field output request
# Change the name of field output request 'F-Output-1' to 'Selected Field Outputs'
specimenModel.fieldOutputRequests.changeKey(fromName='F-Output-1', toName='Selected Field Outputs')
# Since F-Output-1 is applied at the 'Apply Load' step by default, 'Selected Field
# Outputs' will be too
# We only need to set the required variables
specimenModel.fieldOutputRequests['Selected Field Outputs'].setValues(variables=('S', 'E', 'PEMAG', 'U', 'RF' , 'CF'))

# ------------------------------------------------------------------------
# Create the history output request
# We try a slightly different method from that used in field output request
#Create a new history output request called 'Default History Outputs' and assign
# both the step and the variables
specimenModel.HistoryOutputRequest(name='Default History Outputs', createStepName='Apply Displacement Load', variables=PRESELECT)
#Now delete the original history output request 'H-Output-1'
del specimenModel.historyOutputRequests['H-Output-1']

# ------------------------------------------------------------------------
# Apply encastre (fixed) boundary condition to one end to make it cantilever
# First we need to locate and select the top surface
# We place a point somewhere on the top surface based on our knowledge of the
# geometry
fixed_face_pt_x = 0
fixed_face_pt_y = 0.5 * width
fixed_face_pt_z = 0.5 * thickness
fixed_face_pt = (fixed_face_pt_x, fixed_face_pt_y, fixed_face_pt_z)
# The face on which that point lies is the face we are looking for
fixed_face = specimenInstance.faces.findAt((fixed_face_pt,))
# we extract the region of the face choosing which direction its normal points in
fixed_face_region=regionToolset.Region(faces=fixed_face)
specimenModel.EncastreBC(name='Encaster one end', createStepName='Initial', region=fixed_face_region)

# ------------------------------------------------------------------------
# Apply displacement load to top surface
# First we need to locate and select the top surface
# We place a point somewhere on the top surface based on our knowledge of the
# geometry
displacement_table = []
time = 0
amplitude = 0
while time <= 1:
    displacement_table.append([time, amplitude])
    time += 0.1
    amplitude += 0.1
displacement_table = tuple(displacement_table)
specimenModel.TabularAmplitude(name='Displacement Amp', timeSpan=STEP,
        smooth=SOLVER_DEFAULT, data=displacement_table)

end_face_pt_x = length
end_face_pt_y = 0.5 * (width)
end_face_pt_z = 0.5 * (thickness)
end_face_pt = (end_face_pt_x, end_face_pt_y, end_face_pt_z)
# The face on which that point lies is the face we are looking for
end_face = specimenInstance.faces.findAt((end_face_pt,))
# We extract the region of the face choosing which direction its normal points in
end_face_region=regionToolset.Region(faces=end_face)
#Apply the pressure load on this region in the 'Apply Load' step
specimenModel.DisplacementBC(name='displacementBC',
                             createStepName='Apply Displacement Load', region=end_face_region, u1=5.0,
                             u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='Displacement Amp',
                             fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

# # ------------------------------------------------------------------------
# Create the mesh
import mesh
# Assign element type to all cells
elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD,
                          kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF,
                          hourglassControl=DEFAULT, distortionControl=DEFAULT)
specimenCells=specimenPart.cells
specimenMeshRegion=(specimenCells,)
specimenPart.setElementType(regions=specimenMeshRegion, elemTypes=(elemType1,))
edges = specimenPart.edges
# seed cell 1
x_min,y_min,z_min,x_max,y_max,z_max = 0,0,0,a,width,thickness
edges_1 = edges.getByBoundingBox(x_min,y_min,z_min,x_max,y_max,z_max )
specimenPart.seedEdgeBySize(edges=edges_1, size=3.0, deviationFactor=0.1,
        constraint=FINER)
# seed cell 2
x_min,y_min,z_min,x_max,y_max,z_max = a,0,0,a + b,width,thickness
edges_1 = edges.getByBoundingBox(x_min,y_min,z_min,x_max,y_max,z_max )
specimenPart.seedEdgeBySize(edges=edges_1, size=2.5, deviationFactor=0.1,
        constraint=FINER)
# seed cell 3
x_min,y_min,z_min,x_max,y_max,z_max = a+b,0,0,a + b + c,width,thickness
edges_1 = edges.getByBoundingBox(x_min,y_min,z_min,x_max,y_max,z_max )
specimenPart.seedEdgeBySize(edges=edges_1, size=2.0, deviationFactor=0.1,
        constraint=FINER)
# seed cell 4 (gauge cell)
x_min,y_min,z_min,x_max,y_max,z_max = a + b + c,0,0,a + b + c + d,width,thickness
edges_1 = edges.getByBoundingBox(x_min,y_min,z_min,x_max,y_max,z_max )
specimenPart.seedEdgeBySize(edges=edges_1, size=1.0, deviationFactor=0.1,
        constraint=FINER)
# seed cell 5
x_min,y_min,z_min,x_max,y_max,z_max = a + b + c + d,0,0,a + b + 2*c + d,width,thickness
edges_1 = edges.getByBoundingBox(x_min,y_min,z_min,x_max,y_max,z_max )
specimenPart.seedEdgeBySize(edges=edges_1, size=2.0, deviationFactor=0.1,
        constraint=FINER)
# seed cell 6
x_min,y_min,z_min,x_max,y_max,z_max = a + b + 2*c + d,0,0,a + 2*b + 2*c + d,width,thickness
edges_1 = edges.getByBoundingBox(x_min,y_min,z_min,x_max,y_max,z_max )
specimenPart.seedEdgeBySize(edges=edges_1, size=2.5, deviationFactor=0.1,
        constraint=FINER)
# seed cell 7
x_min,y_min,z_min,x_max,y_max,z_max = a + 2*b + 2*c + d,0,0,length,width,thickness
edges_1 = edges.getByBoundingBox(x_min,y_min,z_min,x_max,y_max,z_max )
specimenPart.seedEdgeBySize(edges=edges_1, size=3.0, deviationFactor=0.1,
        constraint=FINER)

specimenPart.generateMesh()

# ------------------------------------------------------------------------
# Create and run the job
import job
# Create the job
mdb.Job(name='SpecimenJob', model='Specimen', type=ANALYSIS,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE,
        description='Job simulates a loaded cantilever beam',
        parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT,
        numDomains=1, userSubroutine='', numCpus=1, memory=50,
        memoryUnits=PERCENTAGE, scratch='', echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF)
# Run the job
mdb.jobs['SpecimenJob'].submit(consistencyChecking=OFF)
# Do not return control till job is finished running
mdb.jobs['SpecimenJob'].waitForCompletion()
# End of run job
# Create the .txt file to show that simulation is finished
finish_file_name = 'simulation_finish.txt'
dispFile = open(finish_file_name, 'w')
dispFile.close()

# ------------------------------------------------------------------------
# # Post processing
# import visualization
# specimen_viewport = session.Viewport(name='Specimen Results Viewport')
# specimen_Odb_Path = 'SpecimenJob.odb'
# an_odb_object = session.openOdb(name=specimen_Odb_Path)
# specimen_viewport.setValues(displayedObject=an_odb_object)
# specimen_viewport.odbDisplay.display.setValues(plotState=(DEFORMED,))