import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def draw_out(binslog, binsper, GD_pho, \
    GD_mean_velpar, GD_mean_velper, \
    GD_std_velpar,  GD_std_velper, \
    LOG=True, Nlevel=20):

  if(LOG == True):
    GD_pho           = np.mean(GD_pho,axis=0)
    pho_min          = np.min(GD_pho[GD_pho>0])
    GD_pho[GD_pho<0] = pho_min
    GD_pho           = np.log10(GD_pho)
  
  GD_mean_velpar = np.mean(GD_mean_velpar,axis=0)
  GD_mean_velper = np.mean(GD_mean_velper,axis=0)
  GD_std_velpar  = np.mean(GD_std_velpar,axis=0)
  GD_std_velper  = np.mean(GD_std_velper,axis=0)
  
  ###############  DELTA  #########################################
  fig = plt.figure()
  fig.subplots_adjust(hspace=0.3)
  plt.ion()
  
  plt.subplot(1, 1, 1)
  
  cset  = plt.contour(binslog,binsper,GD_pho,Nlevel,colors='k')
  cset = plt.contourf(binslog,binsper,GD_pho,Nlevel,cmap=cm.coolwarm)
  cbar = plt.colorbar(cset)
  plt.title('DELTA')
  
  plt.show()
  
  #################################################################
  
  fig = plt.figure()
  fig.subplots_adjust(hspace=0.3)
  
  ###################### MEDIAS ###################################
  plt.subplot(2, 2, 1)
  
  cset  = plt.contour(binslog,binsper,GD_mean_velper,Nlevel, colors='k')
  cset = plt.contourf(binslog,binsper,GD_mean_velper,Nlevel, cmap=cm.coolwarm)
  cbar = plt.colorbar(cset)
  plt.title('MEDIA PERPENDICULAR')
  
  plt.subplot(2, 2, 3)
  
  cset  = plt.contour(binslog,binsper,GD_mean_velpar,Nlevel, colors='k')
  cset = plt.contourf(binslog,binsper,GD_mean_velpar,Nlevel, cmap=cm.coolwarm)
  cbar = plt.colorbar(cset)
  
  plt.title('MEDIA PARALELA')
  #################################################################
  
  ######################  STD  ####################################
  min_std = np.min([GD_std_velper.min(),GD_std_velpar.min()]) 
  max_std = np.max([GD_std_velper.max(),GD_std_velpar.max()]) 
  
  plt.subplot(2, 2, 2)
  
  cset  = plt.contour(binslog,binsper,GD_std_velper,Nlevel, colors='k')
  cset = plt.contourf(binslog,binsper,GD_std_velper,Nlevel, cmap=cm.coolwarm)
  cbar = plt.colorbar(cset)
  plt.title('STD PERPENDICULAR')
  
  
  plt.subplot(2, 2, 4)
  
  cset  = plt.contour(binslog,binsper,GD_std_velpar,Nlevel, colors='k')
  cset = plt.contourf(binslog,binsper,GD_std_velpar,Nlevel, cmap=cm.coolwarm)
  cbar = plt.colorbar(cset)
  
  plt.title('STD PARALELA')
  
  #################################################################
  
  plt.show()


#################################################################
def draw_hist(data_x,x_label='x_data',num_bins=50,log=False):

  plt.ion()
 
  x_min = data_x.min()-1e-10
  x_max = data_x.max()+1e-10
  print "min %g, max %g" % (x_min, x_max)
  if log==False:
    bins = np.linspace(x_min,x_max,num=num_bins)
  else:
    bins = np.logspace(np.log10(x_min),np.log10(x_max),num=num_bins)
  
  f, axs = plt.subplots(1)
  
  axs.hist(data_x,bins,color="r",normed=True)
  axs.set_xlabel(x_label)
  axs.set_xlim([x_min,x_max])
 
  plt.show()


def draw_scatter(data_x,data_y,x_label='x_data',y_label='y_data',x_log=False,y_log = False):
 
  plt.ion()

  fig, axs = plt.subplots(ncols=2, sharey=True)

  assert(len(data_x)==len(data_y))

  if(x_log): data_x = np.log10(data_x)
  if(y_log): data_y = np.log10(data_y)

  axs[0].set_ylabel(y_label)

  axs[0].set_title("Hexagon binning")
  hb = axs[0].hexbin(data_x, data_y, gridsize=50, cmap=cm.coolwarm)
  cb = fig.colorbar(hb, ax=axs[0])
  cb.set_label('N')
  axs[0].set_xlabel(x_label)

  axs[1].set_title("Log color scale")
  hb = axs[1].hexbin(data_x, data_y, gridsize=50, bins='log', cmap=cm.coolwarm)

  cb = fig.colorbar(hb, ax=axs[1])
  cb.set_label('log10(N)')
  axs[1].set_xlabel(x_label)

  fig.subplots_adjust(hspace=0.5)

  plt.show()
