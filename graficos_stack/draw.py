import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.colorbar as bar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
from scipy import interpolate
from scipy import stats
#import matplotlib.gridspec as gridspec

class nf(float):
  def __repr__(self):
    str = '%.1f' % (self.__float__(),)
    if str[-1] == '0':
      return '%.0f' % self.__float__()
    else:
      return '%.1f' % self.__float__()

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def label_pow(x,pos):
  return "{:4.2f}".format(10.0**x)

def jackknife_resampling(data):

    n = data.shape[0]
    assert n > 0, "data must contain at least one measurement"

    resamples = np.empty([n, n-1])

    for i in range(n):
        resamples[i] = np.delete(data, i)

    return resamples


def jackknife_stats(data, statistic, conf_lvl=0.95):
    """
    This function requires `SciPy <http://www.scipy.org>`_ to be installed.

    Parameters
    ----------
    data : numpy.ndarray
        Original sample (1-D array).
    statistic : function
        Any function (or vector of functions) on the basis of the measured
        data, e.g, sample mean, sample variance, etc. The jackknife estimate of
        this statistic will be returned.
    conf_lvl : float, optional
        Confidence level for the confidence interval of the Jackknife estimate.
        Must be a real-valued number in (0,1). Default value is 0.95.

    Returns
    -------
    estimate : numpy.float64 or numpy.ndarray
        The i-th element is the bias-corrected "jackknifed" estimate.

    bias : numpy.float64 or numpy.ndarray
        The i-th element is the jackknife bias.

    std_err : numpy.float64 or numpy.ndarray
        The i-th element is the jackknife standard error.

    conf_interval : numpy.ndarray
        If ``statistic`` is single-valued, the first and second elements are
        the lower and upper bounds, respectively. If ``statistic`` is
        vector-valued, each column corresponds to the confidence interval for
        each component of ``statistic``. The first and second rows contain the
        lower and upper bounds, respectively.
    """

    from scipy.special import erfinv

    # make sure original data is proper
    n = data.shape[0]
    assert n > 0, "data must contain at least one measurement"

    resamples = jackknife_resampling(data)

    stat_data = statistic(data)
    jack_stat = np.apply_along_axis(statistic, 1, resamples)
    mean_jack_stat = np.mean(jack_stat, axis=0)

    # jackknife bias
    bias = (n-1)*(mean_jack_stat - stat_data)

    # jackknife standard error
    std_err = np.sqrt((n-1)*np.mean((jack_stat - mean_jack_stat)*(jack_stat -
                                    mean_jack_stat), axis=0))

    # bias-corrected "jackknifed estimate"
    estimate = stat_data - bias

    # jackknife confidence interval
    assert (conf_lvl > 0 and conf_lvl < 1), "confidence level must be in (0,1)."

    z_score = np.sqrt(2.0)*erfinv(conf_lvl)
    conf_interval = estimate + z_score*np.array((-std_err, std_err))

    return estimate, bias, std_err, conf_interval

def draw_out(binslog, binsper,         \
    GD_pho,         min_max_pho,       \
    GD_mean_velpar, min_max_vmean_par, \
    GD_mean_velper, min_max_vmean_per, \
    GD_std_velpar,  min_max_vrms_par,  \
    GD_std_velper,  min_max_vrms_per,  \
    lista_pho, lista_vel, lista_std,   \
    lin, LOG, flag_vel, title, m, Nlevel = 20):

  #tiene que estar definido dentro de la 
  #funcion para que pueda "ver" la media m
  def label_norm(x,pos):
    return "{:4.2f}".format(x/m)

  name_pdf = title.replace('\in\ ','')
  name_pdf = name_pdf.replace('and\ ','and')
  name_pdf = name_pdf.replace('.',' ')
  name_pdf = name_pdf.replace('$','')
  name_pdf = name_pdf.replace('[','')
  name_pdf = name_pdf.replace(']','')
  name_pdf = name_pdf.replace('(','')
  name_pdf = name_pdf.replace(',','')
  name_pdf = name_pdf.replace('\ ',' ')
  name_pdf = name_pdf.replace(' ','_')

  GD_pho  = np.mean(GD_pho,axis=0)

  if(LOG == True):
    GD_pho = np.log10(GD_pho)
  
  GD_mean_velpar = np.mean(GD_mean_velpar,axis=0)
  GD_mean_velper = np.mean(GD_mean_velper,axis=0)
  GD_std_velpar  = np.mean(GD_std_velpar,axis=0)
  GD_std_velper  = np.mean(GD_std_velper,axis=0)

  # INVIERTO LA MATRIZ
  GD_pho         =         GD_pho[::-1]
  GD_mean_velpar = GD_mean_velpar[::-1]
  GD_mean_velper = GD_mean_velper[::-1]
  GD_std_velpar  =  GD_std_velpar[::-1]
  GD_std_velper  =  GD_std_velper[::-1]

  #####xorg = GD_pho.shape[0]
  #####yorg = GD_pho.shape[1]

  #####assert(yorg%2==0)

  #####mm = yorg/2
  #####binslog = binslog[:mm]
  #####GD_pho = -GD_pho[:,:mm]+np.fliplr(GD_pho[:,mm:])
  #####GD_mean_velpar =  GD_mean_velpar[:,:mm]+np.fliplr(-1.0*GD_mean_velpar[:,mm:])
  #####GD_mean_velper = -GD_mean_velper[:,:mm]+np.fliplr(GD_mean_velper[:,mm:])
  #####GD_std_velpar = -GD_std_velpar [:,:mm]+np.fliplr(GD_std_velpar[:,mm:])
  #####GD_std_velper = -GD_std_velper [:,:mm]+np.fliplr(GD_std_velper[:,mm:])
  #####GD_mean_velperinv = GD_mean_velper

  binsper        = np.hstack((-1.0*np.flipud(binsper  ),binsper))
  GD_mean_velperinv = np.vstack((-1.0*np.flipud(GD_mean_velper),GD_mean_velper))
  GD_pho         = np.vstack((np.flipud(GD_pho        ),GD_pho))
  GD_mean_velpar = np.vstack((np.flipud(GD_mean_velpar),GD_mean_velpar))
  GD_mean_velper = np.vstack((np.flipud(GD_mean_velper),GD_mean_velper))
  GD_std_velpar  = np.vstack((np.flipud(GD_std_velpar ),GD_std_velpar))
  GD_std_velper  = np.vstack((np.flipud(GD_std_velper ),GD_std_velper))
  
  #################################################################
  ###############  DELTA  #########################################
  fig = plt.figure(figsize=(9, 7)) 
  #fig.suptitle(title.replace('\n','  '),x=0.47)

  ax = plt.subplot(1, 1, 1)
  if(LOG == True):
    print "DELTA LOGARITMICA"
    min_max_pho = np.log10(min_max_pho)

  print min_max_pho
  Niveles = np.linspace(min_max_pho[0],min_max_pho[1],Nlevel)

  print binslog.shape
  print binsper.shape
  print GD_pho.shape

  rho_levels  = np.array([1.0, 5.0, 10.0, 100.0, 200.0])       #min_max_pho    = [1., 1000.]
  vper_levels = np.array([-400.0,-300.0, -200.0, -100.0, 0.0]) #min_max_velper = [-500.0, 100.]
  vpar_levels = np.array([-200.0, -100.0, 0.0, 100.0, 200.0])  #min_max_velpar = [-250.0, 250.]
  sper_levels = np.array([50.0, 100.0, 150.0, 200.0, 250.0])   #min_max_rmsper = [0.0, 500.]
  spar_levels = np.array([50.0, 100.0, 150.0, 200.0, 250.0])   #min_max_rmspar = [0.0, 500.]
  plt.rcParams['contour.negative_linestyle'] = 'solid'

  cset = ax.contour(binslog,binsper,GD_pho,levels=np.log10(rho_levels),extend='both',colors='k')
  cset.levels = [nf(10**val) for val in cset.levels]

  if plt.rcParams["text.usetex"]:
    fmt = r'%r'
  else:
    fmt = '%r'

  ax.clabel(cset, cset.levels, inline=True, inline_spacing=6, rightside_up=True, colors='k', fontsize=10, fmt=fmt)
  cset = ax.contourf(binslog,binsper,GD_pho,levels=Niveles,extend='both',cmap=cm.rainbow)
  ax.set_aspect(1.0)

  # create an axes on the right side of ax. The width of cax will be 5%
  # of ax and the padding between cax and ax will be fixed at 0.10 inch
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.10)
  cbar = plt.colorbar(cset, cax=cax)
  cbar.formatter = ticker.FuncFormatter(label_pow)
  cbar.update_ticks()
  #cbar.set_label(r'$\delta + 1$', fontsize=14, rotation=-90, labelpad=20)

  #ax.set_title(title + '\n\n', fontsize=14)
  ax.set_title(title + '\n\n' + '$\delta + 1$', fontsize=14)
  ax.set_xlabel(r'$z\ /\ l$', fontsize=14)
  ax.set_ylabel(r'$r\ [h^{-1} Mpc]$', fontsize=14)
  ax.xaxis.set_major_formatter(ticker.FuncFormatter(label_norm))

  if lin==True:
    plt.savefig("%s_LIN_sin_flechas_%s_rho.pdf" % (name_pdf, flag_vel), format='pdf',bbox_inches='tight')
  else:
    plt.savefig("%s_LOG_sin_flechas_%s_rho.pdf" % (name_pdf, flag_vel), format='pdf',bbox_inches='tight')

  # bicubic interpolation
  xi = np.linspace(binslog.min(), binslog.max(), 50)
  yi = np.linspace(binsper.min(), binsper.max(), 50)

  uCi = interpolate.interp2d(binslog, binsper, GD_mean_velpar   )(xi, yi)
  vCi = interpolate.interp2d(binslog, binsper, GD_mean_velperinv)(xi, yi)
  lw = np.sqrt(uCi**2+vCi**2)
  lw = ax.streamplot(xi,yi,uCi,vCi,density=[1.0, 1.0],color='k',cmap=cm.cool,linewidth=2.0*lw/lw.max())
 
  if lin==True:
    plt.savefig("%s_LIN_%s_rho.pdf" % (name_pdf, flag_vel), format='pdf',bbox_inches='tight')
  else:
    plt.savefig("%s_LOG_%s_rho.pdf" % (name_pdf, flag_vel), format='pdf',bbox_inches='tight')

  lista_pho.append(fig)
  plt.close(fig)
  #################################################################

  ###################### MEDIAS ###################################

  #gs = gridspec.GridSpec(2, 2)
  #ax1 = plt.subplot(gs[0, :])
  #ax2 = plt.subplot(gs[1,:-1])
  #ax3 = plt.subplot(gs[1:, -1])
  #ax4 = plt.subplot(gs[-1,0])
  #ax5 = plt.subplot(gs[-1,-2])

  fig = plt.figure(figsize=(11,8))
  #fig.suptitle(title.replace('\n','  '),x=0.47)
  fig.subplots_adjust(wspace=0.5,hspace=0.5)
  
  ax = plt.subplot(2, 1, 1)
  Niveles = np.linspace(min_max_vmean_per[0],min_max_vmean_per[1],Nlevel)

  norme = MidpointNormalize(midpoint=0.00)
  #cset =  ax.contour(binslog,binsper,GD_mean_velper,levels=Niveles,extend='both',colors='k')
  cset =  ax.contour(binslog,binsper,GD_mean_velper,levels=vper_levels,extend='both',colors='k')
  cset.levels = [val for val in cset.levels]
  if plt.rcParams["text.usetex"]:
    fmt = r'%r'
  else:
    fmt = '%r'
  ax.clabel(cset, cset.levels, inline=True, inline_spacing=10, rightside_up=True, colors='k', fontsize=10, fmt=fmt)
  cset = ax.contourf(binslog,binsper,GD_mean_velper,extend='both',levels=Niveles,norm=norme,cmap=cm.bwr) 
  ax.set_aspect(1.0)
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.10)
  #cbar = plt.colorbar(cset,cax=cax,extend='both')
  cbar = bar.ColorbarBase(cax,cmap=cm.bwr,norm=norme,extend='both',boundaries=Niveles.tolist(),extendfrac='auto',spacing='uniform',orientation='vertical')
  cbar.ax.set_title(r'$s^{-1}\ km$')
  #cbar.set_label(r'$\overline{v_{\perp}}\ [s^{-1} km]$', fontsize=14, rotation=-90, labelpad=20)

  #ax.set_title(title + '\n\n', fontsize=14)
  ax.set_title(title + '\n\n' + r'$\overline{v_{\perp}}$', fontsize=14)
  ax.set_xlabel(r'$z\ /\ l$', fontsize=14)
  ax.set_ylabel(r'$r\ [h^{-1} Mpc]$', fontsize=14)
  ax.xaxis.set_major_formatter(ticker.FuncFormatter(label_norm))

  ax = plt.subplot(2, 1, 2)
  Niveles = np.linspace(min_max_vmean_par[0],min_max_vmean_par[1],Nlevel)

  #cset =  ax.contour(binslog,binsper,GD_mean_velpar,levels=Niveles,extend='both',colors='k')
  cset =  ax.contour(binslog,binsper,GD_mean_velpar,levels=vpar_levels,extend='both',colors='k')
  cset.levels = [val for val in cset.levels]
  if plt.rcParams["text.usetex"]:
    fmt = r'%r'
  else:
    fmt = '%r'
  ax.clabel(cset, cset.levels, inline=True, inline_spacing=5, rightside_up=True, colors='k', fontsize=10, fmt=fmt)
  cset = ax.contourf(binslog,binsper,GD_mean_velpar,levels=Niveles,extend='both',cmap=cm.bwr)
  #cbar = plt.colorbar(cset)
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.10)
  cbar = plt.colorbar(cset, cax=cax)
  cbar.ax.set_title(r'$s^{-1}\ km$')
  #cbar.set_label(r'$\overline{v_{\parallel}}\ [s^{-1}\ km]$', fontsize=14, rotation=-90, labelpad=20)

  ax.set_aspect(1.0)
  ax.set_title('$\overline{v_{\parallel}}$', fontsize=14)
  ax.set_xlabel(r'$z\ /\ l$', fontsize=14)
  ax.set_ylabel(r'$r\ [h^{-1}\ Mpc]$', fontsize=14)
  ax.xaxis.set_major_formatter(ticker.FuncFormatter(label_norm))

  if lin==True:
    plt.savefig("%s_LIN_%s_vel.pdf" % (name_pdf, flag_vel), format='pdf',bbox_inches='tight')
  else:
    plt.savefig("%s_LOG_%s_vel.pdf" % (name_pdf, flag_vel), format='pdf',bbox_inches='tight')

  lista_vel.append(fig)
  plt.close(fig)

  #################################################################
  ######################  STD  ####################################

  fig = plt.figure(figsize=(11,8))
  #fig.suptitle(title.replace('\n','  '),x=0.47)
  fig.subplots_adjust(wspace=0.5,hspace=0.5)

  ax = plt.subplot(2, 1, 1)
  Niveles = np.linspace(min_max_vrms_per[0],min_max_vrms_per[1],Nlevel)
 
  #cset =  ax.contour(binslog,binsper,GD_std_velper,levels=Niveles,extend='both',colors='k')
  cset =  ax.contour(binslog,binsper,GD_std_velper,levels=sper_levels,extend='both',colors='k')
  cset.levels = [val for val in cset.levels]
  if plt.rcParams["text.usetex"]:
    fmt = r'%r'
  else:
    fmt = '%r'
  ax.clabel(cset, cset.levels, inline=True, inline_spacing=7, rightside_up=True, colors='k', fontsize=10, fmt=fmt)
  cset = ax.contourf(binslog,binsper,GD_std_velper,levels=Niveles,extend='both',cmap=cm.bwr)
  #cbar = plt.colorbar(cset)
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.10)
  cbar = plt.colorbar(cset, cax=cax)
  cbar.ax.set_title(r'$s^{-1}\ km$')
  #cbar.set_label(r'$\overline{\sigma_{\perp}}\ [s^{-1} km]$', fontsize=14, rotation=-90, labelpad=20)

  ax.set_aspect(1.0)
  ax.set_title(title + '\n\n' + r'$\overline{\sigma_{\perp}}$', fontsize=14)
  ax.set_xlabel(r'$z\ /\ l $', fontsize=14)
  ax.set_ylabel(r'$r\ [h^{-1} Mpc]$', fontsize=14)
  ax.xaxis.set_major_formatter(ticker.FuncFormatter(label_norm))
 
  ax = plt.subplot(2, 1, 2)
  Niveles = np.linspace(min_max_vrms_par[0],min_max_vrms_par[1],Nlevel)
  
  #cset =  ax.contour(binslog,binsper,GD_std_velpar,levels=Niveles,extend='both',colors='k')
  cset =  ax.contour(binslog,binsper,GD_std_velpar,levels=spar_levels,extend='both',colors='k')
  cset.levels = [val for val in cset.levels]
  if plt.rcParams["text.usetex"]:
    fmt = r'%r'
  else:
    fmt = '%r'
  ax.clabel(cset, cset.levels, inline=True, inline_spacing=7, rightside_up=True, colors='k', fontsize=10, fmt=fmt)
  cset = ax.contourf(binslog,binsper,GD_std_velpar,levels=Niveles,extend='both',cmap=cm.bwr)

  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.10)
  cbar = plt.colorbar(cset, cax=cax)
  cbar.ax.set_title(r'$s^{-1}\ km$')
  #cbar.set_label(r'$\overline{\sigma_{\parallel}}\ (km/s)$', fontsize=14, rotation=-90, labelpad=20)

  ax.set_aspect(1.0)
  ax.set_title(r'$\overline{\sigma_{\parallel}}$', fontsize=14)
  ax.set_xlabel(r'$z\ /\ l$', fontsize=14)
  ax.set_ylabel(r'$r\ [h^{-1} Mpc]$', fontsize=14)
  ax.xaxis.set_major_formatter(ticker.FuncFormatter(label_norm))
  
  if lin==True:
    plt.savefig("%s_LIN_%s_std.pdf" % (name_pdf, flag_vel), format='pdf',bbox_inches='tight')
  else:
    plt.savefig("%s_LOG_%s_std.pdf" % (name_pdf, flag_vel), format='pdf',bbox_inches='tight')

  lista_std.append(fig)
  plt.close(fig)

  #################################################################
  
  return lista_pho, lista_vel, lista_std

#################################################################
def draw_hist(ax,data_x,xlabel,num_bins=20,log=False):

  #x_min = data_x.min()-1e-10
  #x_max = data_x.max()+1e-10
  x_min = 0.0
  x_max = 30.0

  print "min %g, max %g" % (x_min, x_max)
  if log==False:
    bins = np.linspace(x_min,x_max,num=num_bins)
  else:
    bins = np.logspace(np.log10(x_min),np.log10(x_max),num=num_bins)
  
  ax.hist(data_x,bins,color="r",histtype='step',normed=True)
  ax.set_xlim([x_min,x_max])
  ax.set_xlabel(xlabel)

  return

#################################################################
def draw_scatter(data_x,data_y,x_label,y_label,x_log=False,y_log=False):
 
  fig, ax = plt.subplots(ncols=2, sharey=True)

  assert(len(data_x)==len(data_y))

  if(x_log): data_x = np.log10(data_x)
  if(y_log): data_y = np.log10(data_y)

  hb = ax[0].hexbin(data_x, data_y, gridsize=50, cmap=cm.bwr)
  cb = fig.colorbar(hb, ax=ax[0])
  cb.set_label('N')
  ax[0].set_title("Hexagon binning")
  ax[0].set_ylabel(y_label)
  ax[0].set_xlabel(x_label)

  hb = ax[1].hexbin(data_x, data_y, gridsize=50, bins='log', cmap=cm.bwr)
  cb = fig.colorbar(hb, ax=ax[1])
  cb.set_label('log10(N)')
  ax[1].set_title("Log color scale")
  ax[1].set_xlabel(x_label)

  fig.subplots_adjust(hspace=0.5)

  plt.show()

  return

def save_pdf(lista_pho,lista_vel,lista_std,lista_gr,label,num_bins,mean,percent,smooth,extend,fixed,vc,lin,flag_vel):

  bdim = len(percent)
  m    = len(mean)

  lista_pho = np.asarray(lista_pho).reshape(m,bdim,order='C')
  lista_vel = np.asarray(lista_vel).reshape(m,bdim,order='C')
  lista_std = np.asarray(lista_std).reshape(m,bdim,order='C')
  lista_gr  = np.asarray(lista_gr).reshape(m,bdim,order='C')

  if 'long' in label:

    for i in range(m):

      prefix = '%.2f_%.2d' % (mean[i],num_bins)
      if smooth == True: prefix += '_smooth'
      if extend == True: prefix += '_extend'
      if fixed == True: prefix += '_fixed'
      if vc == True: prefix += '_vc'
      if lin == True: 
        prefix += '_lin'
      else:
        prefix += '_log'

      pdf_pho = PdfPages('%s_%s_rho_%s.pdf' % (prefix,flag_vel,label))
      pdf_vel = PdfPages('%s_%s_vel_%s.pdf' % (prefix,flag_vel,label))
      pdf_std = PdfPages('%s_%s_std_%s.pdf' % (prefix,flag_vel,label))
      pdf_gr  = PdfPages('%s_%s_int_%s.pdf' % (prefix,flag_vel,label))

      pdf_gr.savefig(lista_gr[i,0])

      for j in range(bdim):

        pdf_pho.savefig(lista_pho[i,j])
        pdf_vel.savefig(lista_vel[i,j])
        pdf_std.savefig(lista_std[i,j])

      size  = len(lista_gr[i,1])
      error = np.concatenate(lista_gr[i,1:])
      error = np.reshape(error,(len(error)/size,size),order='C')

      x      = error[:,0] 
      dM     = error[:,1] 
      err_dM = error[:,2] 
      M0     = error[:,3] 
      err_M0 = error[:,4] 
      M1     = error[:,5] 
      err_M1 = error[:,6]

      fig = plt.figure(figsize=(9, 7)) 
      ax = plt.subplot(1, 1, 1)
      ax.errorbar(x, dM, yerr=err_dM, fmt='ok', lw=3)
      ax.plot(x, dM,'--k', lw=1)
      #ax.set_ylim([10.0,20.0])
      pdf_gr.savefig(fig)

      fig = plt.figure(figsize=(9, 7)) 
      ax = plt.subplot(1, 1, 1)
      ax.errorbar(x, M1, yerr=err_M1, fmt='ob', lw=3)
      ax.errorbar(x, M0, yerr=err_M0, fmt='og', lw=3)
      ax.plot(x, M1,'--b', lw=1)
      ax.plot(x, M0,'--g', lw=1)
      pdf_gr.savefig(fig)

      pdf_pho.close()
      pdf_vel.close()
      pdf_std.close()
      pdf_gr.close()

  elif 'q' in label:

    print "Nada"
    print "Plots comentados"

    #for i in range(bdim):
    #  
    #  if fixed == True:
    #    pdf_pho = PdfPages('%.2f_%.2d_%.2d_fixed_pho_%s.pdf' % (m,percent[i],num_bins,label))
    #    pdf_vel = PdfPages('%.2f_%.2d_%.2d_fixed_vel_%s.pdf' % (m,percent[i],num_bins,label))
    #    pdf_gr  = PdfPages('%.2f_%.2d_%.2d_fixed_int_%s.pdf' % (m,percent[i],num_bins,label))
    #  else:
    #    pdf_pho = PdfPages('%.2f_%.2d_%.2d_pho_%s.pdf' % (m,percent[i],num_bins,label))
    #    pdf_vel = PdfPages('%.2f_%.2d_%.2d_vel_%s.pdf' % (m,percent[i],num_bins,label))
    #    pdf_gr  = PdfPages('%.2f_%.2d_%.2d_int_%s.pdf' % (m,percent[i],num_bins,label))

    #  for j in range(m):

    #    pdf_pho.savefig(lista_pho[j,i])
    #    pdf_vel.savefig(lista_vel[j,i])
    #    pdf_gr.savefig(lista_gr[j,i])

    #  pdf_pho.close()
    #  pdf_vel.close()

  else:

    print 'label error',label

def realizo_graficos(num_bins, mean, delta_m, rperp, data, MNod, longitud, flag, set_fil, extend, lin, flag_vel, LOG = False, name_input = ''):

    #mask     = (flag==2)
    ##leng     = leng[mask]
    #flag     = flag[mask]
    ##rper     = rper[mask]
    #data     = data[mask] 
    #MNod     = MNod[mask]
    #longitud = longitud[mask]

    data["pho"] += 1.0
    data["pho"][data["pho"]<1.0]   = 1.0
    data["pho"][data["pho"]>1000.] = 1000.
    print data["pho"].shape

    if LOG == True:
      min_max_pho    = [1., 1000.]
    else:
      min_max_pho    = [0, 10.]

    min_max_velper = [-500.0, 100.]
    min_max_velpar = [-250.0, 250.]
    min_max_rmsper = [0.0, 500.]
    min_max_rmspar = [0.0, 500.]

    ##############################################################
    lista_gr, lista_pho, lista_vel, lista_std = [], [], [], []

    for m in mean:
 
      mask_len = (longitud>(m-delta_m))*(longitud<=(m+delta_m))
      Qq = MNod[mask_len][:,0]/MNod[mask_len][:,1]

      set_aux = set_fil[mask_len]
      set_aux_MNod = MNod[mask_len]

      q = np.linspace(0.0,100.0,num_bins+1)
      q = np.percentile(Qq,q); q[0] -= 1e-6; q[-1] += 1e-6

      if(extend == True):
        binslog = np.linspace(-0.5*m,1.5*m,data.shape[2]+1)
      else:
        binslog = np.linspace(0.0,m,data.shape[2]+1)

      if(lin == True):
        binsper = np.linspace(0.0,rperp,data.shape[1]+1)
        print "rint = 0.00"
        print binsper
      else:
        binsper = np.logspace(np.log10(rperp/20.0),np.log10(rperp),data.shape[1]+1)
        print "rint = ",rperp/20.0
        print binsper
      
      binslog = 0.5*(binslog[1:]+binslog[:-1])
      binsper = 0.5*(binsper[1:]+binsper[:-1])

      name  = r'$%s\ Lenght\ \in\ [%.2f, %.2f]$' % (name_input, m-delta_m, m+delta_m)

      draw_out(binslog, binsper, \
      data[mask_len]["pho"],           min_max_pho, \
      data[mask_len]["vmean_par"],  min_max_velpar, \
      data[mask_len]["vmean_per"], min_max_velper, \
      data[mask_len]["vrms_par"],   min_max_rmspar, \
      data[mask_len]["vrms_per"],  min_max_rmsper, \
      lista_pho, lista_vel, lista_std, lin, LOG, \
      flag_vel, name, m, Nlevel=25)

      fig = plt.figure(figsize=(9, 7)) 
      ax = plt.subplot(1, 1, 1)
      fig.suptitle(name,x=0.47)
      draw_hist(ax,set_aux,"Delta en el medio",num_bins=20)
      lista_gr.append(fig)

      print '----------------------\n',m,len(Qq),'\n----------------------'
 
      for i in xrange(num_bins):

        q_mask = ((Qq>q[i])*(Qq<=q[i+1]))

        print '%d-quartil %f %f Num %d' % (i,q[i],q[i+1],len(data[mask_len]["pho"][q_mask]))
        print 'Media Min ',np.mean(MNod[mask_len][q_mask][:,0]),' Media Max ',np.mean(MNod[mask_len][q_mask][:,1])

        #if(data.shape[1]!=data.shape[2]):
        #  binslog = np.linspace(-0.5*m,1.5*m,data.shape[2]+1)
        #else:
        #  binslog = np.linspace(0.0,m,data.shape[2]+1)
        #binsper = np.linspace(0.0,5.0,data.shape[1]+1)
        #
        #binslog = 0.5*(binslog[1:]+binslog[:-1])
        #binsper = 0.5*(binsper[1:]+binsper[:-1])

        print len(binslog)
        print len(binsper)
        print data.shape
        if(i==0):
          name  = r'$%s\ Lenght\ \in\ [%.2f, %.2f]\ and\  %s\ q\ \in\ [%.2f, %.2f]$' % (name_input, m-delta_m, m+delta_m, name_input,q[i], q[i+1])
        else:
          name  = r'$%s\ Lenght\ \in\ [%.2f, %.2f]\ and\  %s\ q\ \in\ (%.2f, %.2f]$' % (name_input, m-delta_m, m+delta_m, name_input,q[i], q[i+1])

        draw_out(binslog, binsper, \
        data[mask_len]["pho"][q_mask],           min_max_pho, \
        data[mask_len]["vmean_par"][q_mask],  min_max_velpar, \
        data[mask_len]["vmean_per"][q_mask], min_max_velper,  \
        data[mask_len]["vrms_par"][q_mask],   min_max_rmspar, \
        data[mask_len]["vrms_per"][q_mask],  min_max_rmsper,  \
        lista_pho, lista_vel, lista_std, lin, LOG, flag_vel,  \
        name, m, Nlevel=25)

        #fig = plt.figure(figsize=(9, 7)) 
        #ax = plt.subplot(1, 1, 1)
        #fig.suptitle(name.replace('\n','  '),x=0.47)
        #draw_hist(ax,set_aux[q_mask],"Delta en el medio",num_bins=20)

        #test_statistic = lambda x: (np.mean(x), np.var(x))
        test_stat = np.median
        estimate_dM, bias, stderr_dM, conf_interval = jackknife_stats(set_aux[q_mask],test_stat)
        estimate_M0, bias, stderr_M0, conf_interval = jackknife_stats(set_aux_MNod[q_mask][:,0],test_stat)
        estimate_M1, bias, stderr_M1, conf_interval = jackknife_stats(set_aux_MNod[q_mask][:,1],test_stat)

        array = [0.5*(q[i]+q[i+1]), estimate_dM, stderr_dM, estimate_M0, stderr_M0, estimate_M1, stderr_M1] 

        lista_gr.append(array)
        #lista_gr.append(fig)

        del q_mask

      #leng     = leng[~mask_len]
      #rper     = rper[~mask_len]
      data     = data[~mask_len] 
      MNod     = MNod[~mask_len]
      longitud = longitud[~mask_len]
      set_fil  = set_fil[~mask_len]
      
      del mask_len,Qq

    #del leng,rper,data,MNod,longitud
    del flag,data,MNod,longitud
    ##############################################################

    return lista_pho,lista_vel,lista_std,lista_gr
