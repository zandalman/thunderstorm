import struct
from types import SimpleNamespace
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable
import const
import os

def read(data_path='../data'):
  ''' 
  Read the simulation data. 

  Args
  data_path: Path of the simulation data directory.
  
  Return
  id_list: List of particle IDs.
  data: List of particle data.
  '''
  format_string = '5i9d'
  event_size = struct.calcsize(format_string)

  data, data_temp = [], {}
  col_name_list = ['id', 'nstep', 'Zelem', 'interaction', 'ion', 'time', 'x', 'y', 'z', 'ener', 'cos_th', 'ener_loss', 'ener_sec', 'ener_loss_sync']
  for col_name in col_name_list: data_temp[col_name] = []

  for filename in os.listdir(data_path):

    with open(os.path.join(data_path, filename), 'rb') as file:
      data_raw = file.read()
          
    num_event = len(data_raw) // event_size
    for i in range(num_event):
      event_data = data_raw[i*event_size:(i+1)*event_size]
      event = struct.unpack(format_string, event_data)
      for j, col_name in enumerate(col_name_list):
        data_temp[col_name_list[j]].append(event[j])
          
  for col_name in col_name_list:
    data_temp[col_name] = np.array(data_temp[col_name])

  id_list = np.sort(list(set(data_temp["id"])))
  for id in id_list:
      cond = data_temp["id"] == id
      data.append(SimpleNamespace(**{col_name: data_temp[col_name][cond] for col_name in col_name_list}))
  
  return id_list, data

def plot_traj(data, var_c=None, unit_l=const.AU, unit_c=const.day, unit_l_label='AU', cbar_label=r'$t$ [${\rm day}$]', cmap='jet', cval=0., do_top=False):
  '''
  Plot the trajectory of a particle.

  Args
  data: Particle data.
  var_c: Color variable.
  unit_l: Length unit [cm].
  unit_c: Color unit.
  unit_l_label: Length unit label.
  cbar_label: Colorbar label.
  cmap: Colormap.
  cval: Color variable bound.
  do_top: Plot the color variable as a function of time on the top.
  '''
  ax = plt.figure().add_subplot(projection='3d')
    
  nstep = np.max(data.nstep+1)
  nseg = min(nstep, 256)
  skip = int(nstep/nseg)

  if type(var_c) == type(None): var_c = data.time
  cmin = min(cval/unit_c, var_c[-1]/unit_c)
  cmax = max(cval/unit_c, var_c[-1]/unit_c)
  norm = Normalize(cmin, cmax)
  color_list = mpl.colormaps[cmap](norm(var_c[::skip]/unit_c))
    
  ax.plot([0], [0], [0], color='black', marker='o')
  ax.plot([0, data.x[0]/unit_l], [0, data.y[0]/unit_l], [0, data.z[0]/unit_l], color=color_list[0], lw=2)
  for i in range(nseg):
    ax.plot(data.x[i*skip:(i+1)*skip+1]/unit_l, data.y[i*skip:(i+1)*skip+1]/unit_l, data.z[i*skip:(i+1)*skip+1]/unit_l, color=color_list[i], lw=2)

  sm = ScalarMappable(norm, cmap=cmap)
  cax = ax.inset_axes([1.1, 0, 0.075, 1])
  cbar = plt.colorbar(sm, cax=cax)
  cbar.set_label(cbar_label, fontsize=14)
    
  ax.set_aspect('equal')
  ax.set_xlabel(r'$x$ [%s]' % unit_l_label, fontsize=14)
  ax.set_ylabel(r'$y$ [%s]' % unit_l_label, fontsize=14)
  ax.set_zlabel(r'$z$ [%s]' % unit_l_label, fontsize=14)
    
  if do_top:
    axtop = ax.inset_axes([0, 1.1, 1.175, 0.5])
    axtop.plot(data.time, var_c/unit_c, lw=2)
    axtop.set_xlabel(r'$t$ [${\rm s}$]', fontsize=14)
    axtop.set_ylabel(cbar_label, fontsize=14)
    