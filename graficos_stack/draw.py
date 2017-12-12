import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages

def draw_out(binslog, binsper,         \
    GD_pho,         min_max_pho,       \
    GD_mean_velpar, min_max_vmean_par, \
    GD_mean_velper, min_max_vmean_per, \
    GD_std_velpar,  min_max_vrms_par,  \
    GD_std_velper,  min_max_vrms_per,  \
    lista_pho, lista_vel,              \
    title = '', LOG = False, Nlevel = 20):

  GD_pho  = np.mean(GD_pho,axis=0)

  if(LOG == True):
    GD_pho = np.log10(GD_pho)
  
  GD_mean_velpar = np.mean(GD_mean_velpar,axis=0)
  GD_mean_velper = np.mean(GD_mean_velper,axis=0)
  GD_std_velpar  = np.mean(GD_std_velpar,axis=0)
  GD_std_velper  = np.mean(GD_std_velper,axis=0)

  #################################################################
  ###############  DELTA  #########################################
  fig = plt.figure(figsize=(9, 7)) 
  fig.suptitle(title.replace('\n','  '),x=0.47)

  ax = plt.subplot(1, 1, 1)
  if(LOG == True):
    min_max_pho = np.log10(min_max_pho)

  Niveles = np.linspace(min_max_pho[0],min_max_pho[1],Nlevel)

  cset = ax.contour(binslog,binsper,GD_pho,levels=Niveles,colors='k')
  cset = ax.contourf(binslog,binsper,GD_pho,levels=Niveles,cmap=cm.seismic)
  plt.quiver(binslog,binsper,GD_mean_velpar, GD_mean_velper,units='xy',pivot='mid')
  cbar = plt.colorbar(cset)
  ax.set_title('DELTA')

  lista_pho.append(fig)
  plt.close(fig)
  #################################################################

  ###################### MEDIAS ###################################
  fig = plt.figure(figsize=(11,8))
  fig.suptitle(title.replace('\n','  '),x=0.47)
  fig.subplots_adjust(wspace=0.5,hspace=0.3)
  
  ax = plt.subplot(2, 2, 1)
  Niveles = np.linspace(min_max_vmean_per[0],min_max_vmean_per[1],Nlevel)
 
  cset =  ax.contour(binslog,binsper,GD_mean_velper,levels=Niveles,colors='k')
  cset = ax.contourf(binslog,binsper,GD_mean_velper,levels=Niveles,cmap=cm.seismic)
  cbar = plt.colorbar(cset)
  ax.set_title('MEDIA PERPENDICULAR')
  
  ax = plt.subplot(2, 2, 3)
  Niveles = np.linspace(min_max_vmean_par[0],min_max_vmean_par[1],Nlevel)

  cset =  ax.contour(binslog,binsper,GD_mean_velpar,levels=Niveles,colors='k')
  cset = ax.contourf(binslog,binsper,GD_mean_velpar,levels=Niveles,cmap=cm.seismic)
  cbar = plt.colorbar(cset)
  ax.set_title('MEDIA PARALELA')

  #################################################################
  ######################  STD  ####################################
  ax = plt.subplot(2, 2, 2)
  Niveles = np.linspace(min_max_vrms_per[0],min_max_vrms_per[1],Nlevel)
 
  cset =  ax.contour(binslog,binsper,GD_std_velper,levels=Niveles,colors='k')
  cset = ax.contourf(binslog,binsper,GD_std_velper,levels=Niveles,cmap=cm.seismic)
  cbar = plt.colorbar(cset)
  ax.set_title('STD PERPENDICULAR')
  
  ax = plt.subplot(2, 2, 4)
  Niveles = np.linspace(min_max_vrms_par[0],min_max_vrms_par[1],Nlevel)
  
  cset =  ax.contour(binslog,binsper,GD_std_velpar,levels=Niveles,colors='k')
  cset = ax.contourf(binslog,binsper,GD_std_velpar,levels=Niveles,cmap=cm.seismic)
  cbar = plt.colorbar(cset)
  ax.set_title('STD PARALELA')

  lista_vel.append(fig)
  plt.close(fig)

  #################################################################
  
  return lista_pho, lista_vel

#################################################################
def draw_hist(data_x,xlabel,num_bins=50,log=False):

  x_min = data_x.min()-1e-10
  x_max = data_x.max()+1e-10
  print "min %g, max %g" % (x_min, x_max)
  if log==False:
    bins = np.linspace(x_min,x_max,num=num_bins)
  else:
    bins = np.logspace(np.log10(x_min),np.log10(x_max),num=num_bins)
  
  f, ax = plt.subplots(1)
  
  ax.hist(data_x,bins,color="r",normed=True)
  ax.set_xlim([x_min,x_max])
  ax.set_xlabel(xlabel)

  plt.show()

  return

#################################################################
def draw_scatter(data_x,data_y,x_label,y_label,x_log=False,y_log=False):
 
  fig, ax = plt.subplots(ncols=2, sharey=True)

  assert(len(data_x)==len(data_y))

  if(x_log): data_x = np.log10(data_x)
  if(y_log): data_y = np.log10(data_y)

  hb = ax[0].hexbin(data_x, data_y, gridsize=50, cmap=cm.seismic)
  cb = fig.colorbar(hb, ax=ax[0])
  cb.set_label('N')
  ax[0].set_title("Hexagon binning")
  ax[0].set_ylabel(y_label)
  ax[0].set_xlabel(x_label)

  hb = ax[1].hexbin(data_x, data_y, gridsize=50, bins='log', cmap=cm.seismic)
  cb = fig.colorbar(hb, ax=ax[1])
  cb.set_label('log10(N)')
  ax[1].set_title("Log color scale")
  ax[1].set_xlabel(x_label)

  fig.subplots_adjust(hspace=0.5)

  plt.show()

  return

def save_pdf(lista_pho,lista_vel,label,num_bins,mean,percent):

  bdim = len(percent)
  m    = len(mean)

  lista_pho = np.asarray(lista_pho).reshape(m,bdim,order='C')
  lista_vel = np.asarray(lista_vel).reshape(m,bdim,order='C')

  if 'q' in label:

    for i in range(bdim):

      pdf_pho = PdfPages('%.2d_%.2d_%.2d_pho_%s.pdf' % (m,percent[i],num_bins,label))
      pdf_vel = PdfPages('%.2d_%.2d_%.2d_vel_%s.pdf' % (m,percent[i],num_bins,label))

      for j in range(m):

        pdf_pho.savefig(lista_pho[j,i])
        pdf_vel.savefig(lista_vel[j,i])

      pdf_pho.close()
      pdf_vel.close()

  elif 'long' in label:

    for i in range(m):

      pdf_pho = PdfPages('%.2f_%.2d_%.2d_pho_%s.pdf' % (mean[i],m,num_bins,label))
      pdf_vel = PdfPages('%.2f_%.2d_%.2d_vel_%s.pdf' % (mean[i],m,num_bins,label))

      for j in range(bdim):

        pdf_pho.savefig(lista_pho[i,j])
        pdf_vel.savefig(lista_vel[i,j])

      pdf_pho.close()
      pdf_vel.close()

  else:

    print 'label error',label

