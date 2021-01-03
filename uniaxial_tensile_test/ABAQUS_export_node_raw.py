# *************************************************************
# Unaxial Tensile Test
# *************************************************************
from abaqus import *
from abaqusConstants import *
import regionToolset
import os
import numpy as np
import sys
from caeModules import *
from driverUtils import executeOnCaeStartup
import json

# Change work direction, pls comment these line if run by ABAQUS GUI
simulation_name = sys.argv[-1]
working_directory = os.path.join(os.getcwd(),simulation_name)
os.chdir(working_directory)

# These lines need to be added to run abaqus no GUI
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=240,
                         height=143)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(referenceRepresentation=ON)

# Read odb file
odb_file = session.openOdb(name='SpecimenJob.odb')

# Export node result
item = 'node_raw'
node_solution_folder_path = os.path.join(os.getcwd(), item)
if os.path.isdir(node_solution_folder_path):
    print('CONTINUE')
else:
    os.makedirs(node_solution_folder_path)

step_name = odb_file.steps.keys()[0]
num_of_step = len(odb_file.steps[step_name].frames)
for step in range(num_of_step):
    # These 2 lines need to be added to run abaqus no GUI
    ########################################################
    session.viewports['Viewport: 1'].setValues(displayedObject=odb_file)
    session.fieldReportOptions.setValues(printTotal=OFF, printMinMax=OFF)
    ########################################################
    file_name = '{}_output_step_{}.txt'.format(item, step)
    file_path = os.path.join(node_solution_folder_path, file_name)
    session.writeFieldReport(fileName=file_path, append=OFF,
                             sortItem='Node Label', odb=odb_file, step=0, frame=step,
                             outputPosition=NODAL, variable=(('U', NODAL), ('RF', NODAL)))

# Export node coordinates
# Create a file
x_y_z_coordinate_file_name = 'label_x_y_z.txt'
x_y_z_coordinate_file_path = os.path.join(node_solution_folder_path, x_y_z_coordinate_file_name)
# Open the file
dispFile = open(x_y_z_coordinate_file_path, 'w')
# Identify the instance you need the coordinates
instance_name = odb_file.rootAssembly.instances.keys()[0]
myInstance = odb_file.rootAssembly.instances[instance_name]
# Identify the total number of nodes
numNodes = len(myInstance.nodes)
# Write the information: node, coord1, coord2, coord3
for curNode in myInstance.nodes:
    dispFile.write('%s %10.10E %10.10E %10.10E\n' % (
    curNode.label, curNode.coordinates[0], curNode.coordinates[1], curNode.coordinates[2]))
dispFile.close()

# Write num of step to data
data = json.load(open('data.json'))
data['num_step'] = num_of_step
data['num_node'] = numNodes
json.dump(data, open('data.json', 'w'))

