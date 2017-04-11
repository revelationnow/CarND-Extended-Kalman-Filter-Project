import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', dest='filename')
args = parser.parse_args()

data = open(args.filename,'r')
gt_points = []
est_points = []


for line in data:
    data_arr = line.split('\t')
    est_points.append([float(data_arr[0]), float(data_arr[1])])
    gt_points.append([float(data_arr[6]), float(data_arr[7])])

data.close()


est_points = np.array(est_points)
gt_points = np.array(gt_points)

fig, ax = plt.subplots()
ax.scatter(est_points[:,0], est_points[:,1], color='r')
ax.scatter(gt_points[:,0], gt_points[:,1], color='b')

plt.show()
