import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap

dict_WB = {
         'red': ((0.0, 0.0, 0.0),
                 (1.0, 1.0, 1.0)),

       'green': ((0.0, 0.0, 0.0),
                 (1.0, 1.0, 1.0)),

       'blue':  ((0.0, 0.0, 0.3),
                 (0.0, 0.0, 0.4),
                 (0.0, 0.0, 0.5),
                 (0.0, 0.0, 0.6),
                 (0.0, 0.0, 0.7),
                 (0.0, 0.0, 0.8),
                 (0.0, 0.0, 1.0),
                 (0.0, 0.0, 0.8),
                 (0.0, 0.0, 0.7),
                 (0.0, 0.0, 0.6),
                 (0.0, 0.0, 0.5),
                 (0.0, 0.0, 0.6),
                 (0.0, 0.0, 0.7),
                 (0.0, 0.0, 0.8),
                 (1.0, 1.0, 1.0)),
       }

def draw_out(binslog, binsper,         \
    GD_pho,         min_max_pho,       \
    GD_mean_velpar, min_max_vmean_par, \
    GD_mean_velper, min_max_vmean_per, \
    GD_std_velpar,  min_max_vrms_par,  \
    GD_std_velper,  min_max_vrms_per,  \
    lista_pho, lista_vel, LOG,         \
    title, Nlevel = 20):

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
    print "DELTA LOGARITMICA"
    min_max_pho = np.log10(min_max_pho)

  print min_max_pho
  Niveles = np.linspace(min_max_pho[0],min_max_pho[1],Nlevel)

  #forest = 0.5*(binslog[1:]+binslog[:-1])
  #tree   = 0.5*(GD_pho[:,1:]+GD_pho[:,:-1])

  cset =  ax.contour(binslog,binsper,GD_pho,levels=Niveles,extend='both',colors='k')
  cset = ax.contourf(binslog,binsper,GD_pho,levels=Niveles,extend='both',cmap=cm.seismic)

  #plt.quiver(binslog,binsper,GD_mean_velpar, GD_mean_velper,units='xy',pivot='middle')

  lw = np.sqrt(GD_mean_velpar**2+GD_mean_velper**2)

  plt.streamplot(binslog,binsper,GD_mean_velpar,GD_mean_velper,density=[2.0, 2.0],
  color=GD_pho,cmap=cm.cool,linewidth=3.0*lw/lw.max())
  #color='cyan'
  #pivot : [ 'tail' | 'mid' | 'middle' | 'tip' ]

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

  cset =  ax.contour(binslog,binsper,GD_mean_velper,levels=Niveles,extend='both',colors='k')
  cset = ax.contourf(binslog,binsper,GD_mean_velper,levels=Niveles,extend='both',cmap=cm.seismic)

  #cmapa = LinearSegmentedColormap('White_to_Blue', dict_WB, 256)
  #cset = ax.contourf(binslog,binsper,GD_mean_velper,levels=Niveles,extend='both',cmap=cmapa)

  cbar = plt.colorbar(cset)
  ax.set_title('MEDIA PERPENDICULAR')
  
  ax = plt.subplot(2, 2, 3)
  Niveles = np.linspace(min_max_vmean_par[0],min_max_vmean_par[1],Nlevel)

  cset =  ax.contour(binslog,binsper,GD_mean_velpar,levels=Niveles,extend='both',colors='k')
  cset = ax.contourf(binslog,binsper,GD_mean_velpar,levels=Niveles,extend='both',cmap=cm.seismic)
  cbar = plt.colorbar(cset)
  ax.set_title('MEDIA PARALELA')

  #################################################################
  ######################  STD  ####################################
  ax = plt.subplot(2, 2, 2)
  Niveles = np.linspace(min_max_vrms_per[0],min_max_vrms_per[1],Nlevel)
 
  cset =  ax.contour(binslog,binsper,GD_std_velper,levels=Niveles,extend='both',colors='k')
  cset = ax.contourf(binslog,binsper,GD_std_velper,levels=Niveles,extend='both',cmap=cm.seismic)
  cbar = plt.colorbar(cset)
  ax.set_title('STD PERPENDICULAR')
  
  ax = plt.subplot(2, 2, 4)
  Niveles = np.linspace(min_max_vrms_par[0],min_max_vrms_par[1],Nlevel)
  
  cset =  ax.contour(binslog,binsper,GD_std_velpar,levels=Niveles,extend='both',colors='k')
  cset = ax.contourf(binslog,binsper,GD_std_velpar,levels=Niveles,extend='both',cmap=cm.seismic)
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
  
  ax.hist(data_x,bins,color="r",histtype='step',normed=True)
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

def save_pdf(lista_pho,lista_vel,label,num_bins,mean,percent,fixed):

  bdim = len(percent)
  m    = len(mean)

  lista_pho = np.asarray(lista_pho).reshape(m,bdim,order='C')
  lista_vel = np.asarray(lista_vel).reshape(m,bdim,order='C')

  if 'q' in label:

    for i in range(bdim):
      
      if fixed == True:
        pdf_pho = PdfPages('%.2f_%.2d_%.2d_fixed_pho_%s.pdf' % (m,percent[i],num_bins,label))
        pdf_vel = PdfPages('%.2f_%.2d_%.2d_fixed_vel_%s.pdf' % (m,percent[i],num_bins,label))
      else:
        pdf_pho = PdfPages('%.2f_%.2d_%.2d_pho_%s.pdf' % (m,percent[i],num_bins,label))
        pdf_vel = PdfPages('%.2f_%.2d_%.2d_vel_%s.pdf' % (m,percent[i],num_bins,label))


        
      for j in range(m):

        pdf_pho.savefig(lista_pho[j,i])
        pdf_vel.savefig(lista_vel[j,i])

      pdf_pho.close()
      pdf_vel.close()

  elif 'long' in label:

    for i in range(m):

      if fixed == True:
        pdf_pho = PdfPages('%.2f_%.2d_fixed_pho_%s.pdf' % (mean[i],num_bins,label))
        pdf_vel = PdfPages('%.2f_%.2d_fixed_vel_%s.pdf' % (mean[i],num_bins,label))
      else:
        pdf_pho = PdfPages('%.2f_%.2d_pho_%s.pdf' % (mean[i],num_bins,label))
        pdf_vel = PdfPages('%.2f_%.2d_vel_%s.pdf' % (mean[i],num_bins,label))

      for j in range(bdim):

        pdf_pho.savefig(lista_pho[i,j])
        pdf_vel.savefig(lista_vel[i,j])

      pdf_pho.close()
      pdf_vel.close()

  else:

    print 'label error',label

def realizo_graficos(num_bins, mean, leng, flag, rper, data, MNod, longitud, name_input = '', LOG = False):

    mask     = (longitud>mean[0]-1.0)*(longitud<mean[-1]+1.0)
    leng     = leng[mask]
    flag     = flag[mask]
    rper     = rper[mask]
    data     = data[mask] 
    MNod     = MNod[mask]
    longitud = longitud[mask]


    data["pho"] += 1.0
    data["pho"][data["pho"]<1.0]   = 1.0
    data["pho"][data["pho"]>1000.] = 1000.

    print data["pho"].shape

    min_max_pho    = [np.inf, -np.inf]
    min_max_velpar = [np.inf, -np.inf]
    min_max_velper = [np.inf, -np.inf]
    min_max_rmspar = [np.inf, -np.inf]
    min_max_rmsper = [np.inf, -np.inf]

    for m in mean:
 
      mask_len = (longitud>m-1.0)*(longitud<m+1.0)

      data_mask = np.mean(data["pho"][mask_len],axis=0)
      min_max_pho[0] = min_max_pho[0] if np.min(data_mask)>min_max_pho[0] else np.min(data_mask)
      min_max_pho[1] = min_max_pho[1] if np.max(data_mask)<min_max_pho[1] else np.max(data_mask)

    #  data_mask = np.mean(data["vmean_par"][mask_len],axis=0)
    #  min_max_velpar[0] = min_max_velpar[0] if np.min(data_mask)>min_max_velpar[0] else np.min(data_mask)
    #  min_max_velpar[1] = min_max_velpar[1] if np.max(data_mask)<min_max_velpar[1] else np.max(data_mask)

    #  data_mask = np.mean(data["vmean_per"][mask_len],axis=0)
    #  min_max_velper[0] = min_max_velper[0] if np.min(data_mask)>min_max_velper[0] else np.min(data_mask)
    #  min_max_velper[1] = min_max_velper[1] if np.max(data_mask)<min_max_velper[1] else np.max(data_mask)

    #  data_mask = np.mean(data["vrms_par"][mask_len],axis=0)
    #  min_max_rmspar[0] = min_max_rmspar[0] if np.min(data_mask)>min_max_rmspar[0] else np.min(data_mask)
    #  min_max_rmspar[1] = min_max_rmspar[1] if np.max(data_mask)<min_max_rmspar[1] else np.max(data_mask)

    #  data_mask = np.mean(data["vrms_per"][mask_len],axis=0)
    #  min_max_rmsper[0] = min_max_rmsper[0] if np.min(data_mask)>min_max_rmsper[0] else np.min(data_mask)
    #  min_max_rmsper[1] = min_max_rmsper[1] if np.max(data_mask)<min_max_rmsper[1] else np.max(data_mask)

    if LOG == True:
      min_max_pho    = [1., 1000.]
    else:
      min_max_pho    = [0, 10.]

    min_max_velper = [-300.0,   0.]
    min_max_velpar = [-250.0, 250.]
    min_max_rmsper = [0.0, 500.]
    min_max_rmspar = [0.0, 500.]
    #  min_max_rmsper[1] = min_max_rmsper[1] if np.max(data_mask)<min_max_rmsper[1] else np.max(data_mask)

    #print min_max_pho 
    #print min_max_velpar
    #print min_max_velper
    #print min_max_rmspar
    #print min_max_rmsper

    del flag
    del mask
    ##############################################################

    ##############################################################
    lista_pho, lista_vel = [], []

    for m in mean:
 
      mask_len = (longitud>m-1.0)*(longitud<m+1.0)
      Qq = MNod[mask_len][:,0]/MNod[mask_len][:,1]

      #draw_hist(np.concatenate((MNod[mask_len][:,0],MNod[mask_len][:,1])),"Masa")

      q = np.linspace(0.0,100.0,num_bins+1)
      q = np.percentile(Qq,q); q[0] -= 1e-6; q[-1] += 1e-6

      print '----------------------\n',m,len(Qq),'\n----------------------'
 
      for i in xrange(num_bins):

        q_mask = ((Qq>q[i])*(Qq<=q[i+1]))

        print '%d-quartil %f %f Num %d' % (i,q[i],q[i+1],len(data[mask_len]["pho"][q_mask]))
        print 'Media Min ',np.mean(MNod[mask_len][q_mask][:,0]),' Media Max ',np.mean(MNod[mask_len][q_mask][:,1])

        binslog = np.mean(leng[mask_len][q_mask],axis=0).ravel()
        binsper = np.mean(rper[mask_len][q_mask],axis=0).ravel()  

        name  = 'long > %.2f && long < %.2f\n' % (m-1.0, m+1.0) \
        + name_input + '\nq > %.2f && q <= %.2f' % (q[i], q[i+1])

        draw_out(binslog, binsper, \
        data[mask_len]["pho"][q_mask],           min_max_pho, \
        data[mask_len]["vmean_par"][q_mask],  min_max_velpar, \
        data[mask_len]["vmean_per"][q_mask], min_max_velper, \
        data[mask_len]["vrms_par"][q_mask],   min_max_rmspar, \
        data[mask_len]["vrms_per"][q_mask],  min_max_rmsper, \
        lista_pho, lista_vel, LOG, name, Nlevel=25)

        del q_mask

      leng     = leng[~mask_len]
      rper     = rper[~mask_len]
      data     = data[~mask_len] 
      MNod     = MNod[~mask_len]
      longitud = longitud[~mask_len]
      
      del mask_len,Qq

    del leng,rper,data,MNod,longitud
    ##############################################################

    return lista_pho,lista_vel
