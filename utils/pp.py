import struct
from types import SimpleNamespace
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable
from datetime import datetime
import const
import os

col_event = SimpleNamespace(id=0, nstep=1, Zelem=2, interaction=3, ion=4, time=5, x=6, y=7, z=8, ener=9, cos_th=10, ener_loss=11, ener_sec=12, ener_loss_sync=13, ener_loss_cher=14, ener_loss_moller=15)
col_interact = SimpleNamespace(death=0, scat=1, brem=2, exc=3, ion=4, bturb=5, moller=6)
format_string = '5i11d'
event_size = struct.calcsize(format_string)

def clear_fig():
    for filename in os.listdir(os.path.join('..', 'figures', 'current')): 
        os.remove(os.path.join('..', 'figures', 'current', filename))

def save_fig(fig_name, filetype="png", dpi=256):
  '''
  Save the current matplotlib figure.

  Args
  name (string): figure name
  filetype (string): file type
  dpi (int): dots per inch
  '''
  datetime_string = datetime.now().strftime("%m%d%Y%H%M")
  filename = "%s-%s.%s" % (fig_name, datetime_string, filetype)
  plt.savefig(os.path.join('..', 'figures', 'current', filename), bbox_inches="tight", dpi=dpi)
  print("Saved figure as '%s'" % filename)

def read_chunk(file, chunk_size, lim=None):
  '''
  Read a file chunk by chunk.

  Args
  file: File object.
  chunk_size: Chunk size in bytes.
  lim: Number of chunks to read.
  '''
  count = 0
  while type(lim)==type(None) or count < lim:
    data = file.read(chunk_size)
    if not data: break # stop if there is no data
    count += 1
    yield data

def processing(func):
  '''
  Converts a function on one chunk into a function on the whole data set.

  Args
  shape_data: Shape of the data returned by the function.
  '''
  def processing_func(data_path='../data', chunk_size=128, chunk_num=128, **kwargs):
    '''
    Function to process the data.

    Args
    data_path:  Path to the data.
    chunk_size: Number of events per chunk.
    chunk_num:  Number of chunks.
    kwargs:     Function arguments.

    Return
    id_list: List of particle IDs.
    data_list: List of particle data.
    '''
    id_list, data_list = [], []  
    for filename in os.listdir(data_path):
      if filename == 'info.txt': continue
      file = open(os.path.join(data_path, filename), 'rb') # open file
      for data_raw in read_chunk(file, chunk_size*event_size, lim=chunk_num):          
        data_1chunk = []
        for j in range(len(data_raw) // event_size):
          data_raw_1event = data_raw[j*event_size:(j+1)*event_size]
          data_1event = struct.unpack(format_string, data_raw_1event)
          id = data_1event[0]     
          if id in id_list:
            data_1chunk.append(data_1event)
          else:
            if len(data_1chunk)>0:
              if type(data_list[-1])==type(None):
                data_list[-1] = func(np.array(data_1chunk), **kwargs)
              else:
                data_list[-1] += func(np.array(data_1chunk), **kwargs)
            id_list.append(id)
            data_list.append(None)         
        if len(data_1chunk)>0:
          if type(data_list[-1])==type(None):
            data_list[-1] = func(np.array(data_1chunk), **kwargs)
          else:
            data_list[-1] += func(np.array(data_1chunk), **kwargs)
        if type(chunk_num) != type(None): chunk_num = chunk_num - 1
      file.close()
    return np.array(id_list), np.array(data_list)
  return processing_func

@processing
def ener_loss_mech(data_1chunk):
  data = np.zeros((6))
  for i, interaction in enumerate([col_interact.brem, col_interact.exc, col_interact.ion, col_interact.moller]):
    cond = data_1chunk[:, col_event.interaction] == interaction
    data[i] = np.sum(data_1chunk[:, col_event.ener_loss][cond])
  data[3] += np.sum(data_1chunk[:, col_event.ener_loss_moller])
  data[4] = np.sum(data_1chunk[:, col_event.ener_loss_sync])
  data[5] = np.sum(data_1chunk[:, col_event.ener_loss_cher])
  return data

@processing
def num_ion(data_1chunk, Zelem_bins=np.arange(0.5, 99.5)):
  data, _ = np.histogram(data_1chunk[:, col_event.Zelem], bins=Zelem_bins)
  return data

@processing
def num_sec(data_1chunk, ener_bins=np.logspace(0, 5, 64)):
  interaction = data_1chunk[:, col_event.interaction]
  ener_sec = data_1chunk[:, col_event.ener_sec][interaction==col_interact.ion]
  data, _ = np.histogram(ener_sec, bins=ener_bins)
  return data

@processing
def ener_time(data_1chunk, time_bins=np.linspace(0, 3600, 64)):
  ener = data_1chunk[:, col_event.ener]
  nevent_bin, _ = np.histogram(data_1chunk[:, col_event.time], bins=time_bins)
  ener_bin, _ = np.histogram(data_1chunk[:, col_event.time], bins=time_bins, weights=ener)
  data = np.array([nevent_bin, ener_bin])
  return data

@processing
def ener_loss_time(data_1chunk, time_bins=np.linspace(0, 3600, 64)):
  ener_loss = data_1chunk[:, col_event.ener_loss] + data_1chunk[:, col_event.ener_loss_sync] + data_1chunk[:, col_event.ener_loss_cher] + data_1chunk[:, col_event.ener_loss_moller]
  data, _ = np.histogram(data_1chunk[:, col_event.time], bins=time_bins, weights=ener_loss)
  return data

@processing
def ener_loss_dis(data_1chunk, dis_bins=np.linspace(0, 0.01*const.AU, 64)):
  ener_loss = data_1chunk[:, col_event.ener_loss] + data_1chunk[:, col_event.ener_loss_sync] + data_1chunk[:, col_event.ener_loss_cher] + data_1chunk[:, col_event.ener_loss_moller]
  dis = np.sqrt(data_1chunk[:, col_event.x]**2 + data_1chunk[:, col_event.y]**2 + data_1chunk[:, col_event.z]**2)
  data, _ = np.histogram(dis, bins=dis_bins, weights=ener_loss)
  return data

def read_1part(data_path='../data', chunk_size=128, idx_file=0):
  ''' 
  Read event data for one particle.

  Args
  data_path: Path to the data.
  
  Return
  data: Particle data.
  '''
  id = None
  col_name_list = ['id', 'nstep', 'Zelem', 'interaction', 'ion', 'time', 'x', 'y', 'z', 'ener', 'cos_th', 'ener_loss', 'ener_sec', 'ener_loss_sync', 'ener_loss_cher', 'ener_loss_moller']
  data = {col_name: [] for col_name in col_name_list}

  filename = [filename for filename in os.listdir(data_path) if filename != 'info.txt'][idx_file]
  file = open(os.path.join(data_path, filename), 'rb')
  for data_raw in read_chunk(file, chunk_size*event_size, lim=None):
    for i in range(len(data_raw) // event_size):
      data_raw_1event = data_raw[i*event_size:(i+1)*event_size]
      data_1event = struct.unpack(format_string, data_raw_1event)
      if type(id)==type(None): id = data_1event[col_event.id]
      if data_1event[col_event.id] != id: break
      for j, col_name in enumerate(col_name_list):
        data[col_name_list[j]].append(data_1event[j])
    if data_1event[col_event.id] != id: break
          
  for col_name in col_name_list:
    data[col_name] = np.array(data[col_name])
  data = SimpleNamespace(**data)
  
  return data

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
    