import os
import numpy as np
import json
import pandas as pd
import matplotlib.pyplot as plt


class uniaxial_tensile_test():
    def __init__(self, simulation_name):
        self.name = simulation_name
        self.simulation_path = None
        self.create_simulation_path()
        # Geometry
        self.thickness = 2.0
        # Elasticity
        self.density = 7.8E-09
        self.young_modulus = 208000
        self.nuy = 0.3
        # Plasticity
        self.s0 = 490
        self.K = 928.4
        self.t = 147.28
        self.h = 0.345
        self.data_path = None
        self.save_initial_data()

    def create_simulation_path(self):
        current_dir = os.getcwd()
        if os.path.isdir(os.path.join(current_dir, self.name)):
            print('WARNING: SIMULATION FOLDER EXISTED, OLD FOLDER WILL BE OVERWRITTEN')
            input_command = str(input('Continue y/n: '))
            if input_command == 'y':
                self.simulation_path = os.path.join(current_dir, self.name)
            else:
                print('QUIT')
                quit()
        else:
            self.simulation_path = os.path.join(current_dir, self.name)
            os.makedirs(self.simulation_path)

    def save_initial_data(self):
        data = {'name': self.name,
                'thickness': self.thickness,
                'density': self.density,
                'young_modulus': self.young_modulus,
                'nuy': self.nuy,
                's0': self.s0,
                'K': self.K,
                't': self.t,
                'h': self.h}
        self.data_path = os.path.join(self.simulation_path, 'data.json')
        json.dump(data, open(self.data_path, 'w'))

    def process(self):
        finish_file_path = os.path.join(os.getcwd(), self.name, 'simulation_finish.txt')
        if os.path.exists(finish_file_path):
            os.remove(finish_file_path)
        abaqus_script = 'ABAQUS_uniaxial_tensile_test.py -- {}'.format(self.name)
        print('RUN SIMULATION')
        os.system('abaqus cae noGUI={}'.format(abaqus_script))
        print('FINISH SIMULATION')

    def export_node_raw(self):
        finish_file_path = os.path.join(os.getcwd(), self.name, 'simulation_finish.txt')
        while not os.path.exists(finish_file_path):
            print('ERROR: THE SIMULATION IS NOT FINISHED')
        abaqus_script = 'ABAQUS_export_node_raw.py -- {}'.format(self.name)
        print('EXPORTING RAW DATA')
        os.system('abaqus cae noGUI={}'.format(abaqus_script))
        print('FINISH EXPORTING RAW DATA')

    def convert_node_raw2csv(self):
        print('CONVERTING RAW DATA TO CSV')
        # set up some paths
        current_dir = os.getcwd()
        raw_data_folder_path = os.path.join(current_dir, self.name, "node_raw")
        csv_data_path = os.path.join(current_dir, self.name, "node_csv")
        if os.path.isdir(csv_data_path):
            print('CONTINUE')
        else:
            os.makedirs(csv_data_path)

        # get num step
        data = json.load(open(self.data_path))
        num_step = data['num_step']

        # Load xyz coordinate data
        xyz_header = ['label', 'x_coord', 'y_coord', 'z_coord']
        x_y_coord_file_path = os.path.join(raw_data_folder_path, 'label_x_y_z.txt')
        x_y_z_csv = pd.read_csv(x_y_coord_file_path, sep=" ", header=None)
        x_y_z_csv.columns = xyz_header

        # remove unnecessary information in text result data of abaqus file
        for step in range(num_step):
            txt_file_path = os.path.join(raw_data_folder_path, 'Node_raw_output_step_%s.txt' % step)
            txt_file = open(txt_file_path, 'r')
            lines = txt_file.readlines()
            txt_file.close()
            new_txt_file_path = os.path.join(raw_data_folder_path, 'new_Node_output_step_%s.txt' % step)
            new_txt_file = open(new_txt_file_path, 'w')
            for idx, line in enumerate(lines):
                if idx < 19:
                    continue
                elif idx >= 19 + len(x_y_z_csv[['label']]) + 1:
                    break
                else:
                    new_txt_file.write(line)
            new_txt_file.close()

        # Create a blank frame of full data
        full_data = pd.DataFrame(
            columns=['step', 'label', 'x_coord', 'y_coord', 'z_coord',
                     'Umagnitude', 'U1', 'U2', 'U3',
                     'RFmagnitude', 'RF1', 'RF2', 'RF3'])
        full_node_output_file_path = os.path.join(csv_data_path, 'full_node_output.csv')
        # Loop for full data = raw data + xyz_coordinate
        for step in range(num_step):
            # Load raw data
            txt_file_path = os.path.join(raw_data_folder_path, 'new_Node_output_step_%s.txt' % step)
            raw_header = ['label',
                          'Umagnitude', 'U1', 'U2', 'U3',
                          'RFmagnitude', 'RF1', 'RF2', 'RF3']
            # Here \s+ means any one or more white space character.
            csv_file = pd.read_csv(txt_file_path, sep='\s+', header=None)
            csv_file = csv_file.dropna()
            csv_file.columns = raw_header
            # Append column xyz coordinate to raw data
            csv_file['x_coord'] = x_y_z_csv['x_coord']
            csv_file['y_coord'] = x_y_z_csv['y_coord']
            csv_file['z_coord'] = x_y_z_csv['z_coord']
            csv_file['step'] = step
            csv_file['label'] = x_y_z_csv['label']
            # Append new raw data to full data and save it
            csv_file = csv_file[['step', 'label', 'x_coord', 'y_coord', 'z_coord',
                                 'Umagnitude', 'U1', 'U2', 'U3',
                                 'RFmagnitude', 'RF1', 'RF2', 'RF3']]
            full_data = pd.concat((full_data, csv_file))
            full_data.to_csv(full_node_output_file_path, index=False, header=True)
        print('FINISH CONVERTING RAW DATA TO CSV')

    def visualization(self):
        print('VISUALING')
        current_dir = os.getcwd()
        node_data_path = os.path.join(os.path.join(current_dir, self.name,
                                                   "node_csv",
                                                   'full_node_output.csv'))
        node_data = pd.read_csv(node_data_path)
        header = ['step', 'label', 'x_coord', 'y_coord', 'z_coord',
                  'Umagnitude', 'U1', 'U2', 'U3',
                  'RFmagnitude', 'RF1', 'RF2', 'RF3']
        initial_data = json.load(open(self.data_path))
        num_step = initial_data['num_step']

        # Visualize the Reaction Force at the end face of the specimen
        specimen_length = initial_data['specimen_length']
        gauge_length = initial_data['gauge_length']
        gauge_point1_x = (specimen_length - gauge_length) * 0.5 + 5
        gauge_point2_x = (specimen_length - gauge_length) * 0.5 + 55
        end_face_data = node_data[node_data['x_coord'] == specimen_length]
        gauge_point1_data = node_data[node_data['x_coord'] == gauge_point1_x]
        gauge_point2_data = node_data[node_data['x_coord'] == gauge_point2_x]
        RF1_list = []
        gauge_length_displacement_list = []
        for step in range(num_step):
            # Collect displacement data at end face
            i_step_end_face_data = end_face_data[end_face_data['step'].between(step - 0.5, step + 0.5)]
            i_step_RF1_list = i_step_end_face_data['RF1'].values
            i_step_sum_RF1 = sum(i_step_RF1_list)
            RF1_list.append(i_step_sum_RF1)

            # Collect gauge data
            i_step_gauge_point1_data = gauge_point1_data[gauge_point1_data['step'].between(step - 0.5, step + 0.5)]
            i_step_gauge_point2_data = gauge_point2_data[gauge_point2_data['step'].between(step - 0.5, step + 0.5)]
            U1_gauge_point1 = np.mean(i_step_gauge_point1_data['U1'].values)
            U1_gauge_point2 = np.mean(i_step_gauge_point2_data['U1'].values[0])
            gauge_length_displacement = U1_gauge_point2 - U1_gauge_point1
            gauge_length_displacement_list.append(gauge_length_displacement)

        fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))
        ax1.plot(gauge_length_displacement_list, RF1_list)
        # ax1.set_yscale('log')
        ax1.set_xlabel("Gauge Length Displacement U1 (mm)")
        ax1.set_ylabel("Load (N)")
        ax1.set_title('Load according to gauge length displacement U1')
        fig_path = os.path.join(self.simulation_path,'Load_displacement.png')
        fig.savefig(fig_path)
        print('FINISH VISUALING')


if __name__ == '__main__':
    test1 = uniaxial_tensile_test(simulation_name='test1')
    test1.process()
    test1.export_node_raw()
    test1.convert_node_raw2csv()
    test1.visualization()
