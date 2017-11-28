import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def draw_out(binslog, binsper, GD_pho, \
    GD_mean_velpar, GD_mean_velper, \
    GD_std_velpar,  GD_std_velper, \
    title='',LOG=True, Nlevel=20):

  if(LOG == True):
    GD_pho           = np.mean(GD_pho,axis=0)
    pho_min          = np.min(GD_pho[GD_pho>0])
    GD_pho[GD_pho<0] = pho_min
    GD_pho           = np.log10(GD_pho)
  
  GD_mean_velpar = np.mean(GD_mean_velpar,axis=0)
  GD_mean_velper = np.mean(GD_mean_velper,axis=0)
  GD_std_velpar  = np.mean(GD_std_velpar,axis=0)
  GD_std_velper  = np.mean(GD_std_velper,axis=0)
  
  #################################################################
  ###############  DELTA  #########################################
  fig = plt.figure()
  fig.subplots_adjust(hspace=0.3)
  
  plt.subplot(1, 1, 1)
  if(LOG == True):
    min_max = [np.log10(0.01), np.log10(5000)]
  else:
    min_max = [0.00, 5000]

  Niveles = np.linspace(min_max[0],min_max[1],Nlevel)

  cset =  plt.contour(binslog,binsper,GD_pho,levels=Niveles,colors='k',vmin=min_max[0],vmax=min_max[1])
  cset = plt.contourf(binslog,binsper,GD_pho,levels=Niveles,cmap=cm.coolwarm,vmin=min_max[0],vmax=min_max[1])
  cbar = plt.colorbar(cset)
  plt.title('DELTA\n %s' % title)

  #################################################################
  ###################### MEDIAS ###################################
  fig = plt.figure()
  fig.subplots_adjust(hspace=0.3)
  
  plt.subplot(2, 2, 1)
  min_max = [-400, 50]
  Niveles = np.linspace(min_max[0],min_max[1],Nlevel)
 
  cset =  plt.contour(binslog,binsper,GD_mean_velper,levels=Niveles,colors='k')
  cset = plt.contourf(binslog,binsper,GD_mean_velper,levels=Niveles,cmap=cm.coolwarm)
  cbar = plt.colorbar(cset)
  plt.title('MEDIA PERPENDICULAR\n %s' % title)
  
  plt.subplot(2, 2, 3)
  min_max = [-200, 300]
  Niveles = np.linspace(min_max[0],min_max[1],Nlevel)

  cset =  plt.contour(binslog,binsper,GD_mean_velpar,levels=Niveles,colors='k')
  cset = plt.contourf(binslog,binsper,GD_mean_velpar,levels=Niveles,cmap=cm.coolwarm)
  cbar = plt.colorbar(cset); plt.clim(-150.,200)
  plt.title('MEDIA PARALELA\n %s' % title)

  #################################################################
  ######################  STD  ####################################
  plt.subplot(2, 2, 2)
  min_max = [0, 500]
  Niveles = np.linspace(min_max[0],min_max[1],Nlevel)
 
  cset =  plt.contour(binslog,binsper,GD_std_velper,levels=Niveles,colors='k')
  cset = plt.contourf(binslog,binsper,GD_std_velper,levels=Niveles,cmap=cm.coolwarm)
  cbar = plt.colorbar(cset)
  plt.title('STD PERPENDICULAR\n %s' % title)
  
  plt.subplot(2, 2, 4)
  min_max = [0, 500]
  Niveles = np.linspace(min_max[0],min_max[1],Nlevel)
  
  cset =  plt.contour(binslog,binsper,GD_std_velpar,levels=Niveles,colors='k')
  cset = plt.contourf(binslog,binsper,GD_std_velpar,levels=Niveles,cmap=cm.coolwarm)
  cbar = plt.colorbar(cset)
  plt.title('STD PARALELA\n %s' % title)
  #################################################################

  plt.show()

#################################################################
def draw_hist(data_x,num_bins=50,log=False):

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
  axs.set_xlim([x_min,x_max])
 
  return axs

#################################################################
def draw_scatter(data_x,data_y,x_label='x_data',y_label='y_data',x_log=False,y_log = False):
 
  fig, axs = plt.subplots(ncols=2, sharey=True)

  assert(len(data_x)==len(data_y))

  if(x_log): data_x = np.log10(data_x)
  if(y_log): data_y = np.log10(data_y)

  hb = axs[0].hexbin(data_x, data_y, gridsize=50, cmap=cm.coolwarm)
  cb = fig.colorbar(hb, ax=axs[0])
  cb.set_label('N')

  hb = axs[1].hexbin(data_x, data_y, gridsize=50, bins='log', cmap=cm.coolwarm)

  cb = fig.colorbar(hb, ax=axs[1])
  cb.set_label('log10(N)')
  axs[1].set_xlabel(x_label)

  fig.subplots_adjust(hspace=0.5)

  return axs
