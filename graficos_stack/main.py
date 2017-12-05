import numpy as np
import lector as lee
import draw as graf
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys

if __name__ == '__main__':

  Path_prop, Name_prop, \
  Path_cil, Name_cil, \
  Nfile_cil, Ncil \
  = lee.init(sys.argv[1])

  POSFACTOR = 1000.  
  print 'Numero de cilindros ',Ncil

  Prop = lee.read_prop(Path_prop,Name_prop,POSFACTOR)
  
  mask_prop = (Prop["flag"]>=2)
  print "Mask Segmentos",sum(mask_prop)
  #mask_prop = (Prop["size"]>2)*mask_prop
  #print "Segmentos mask N>2",len(Prop["size"][mask_prop])

  #plt.ion()
  #ax = graf.draw_scatter(Prop["len"][mask_prop],Prop["elong"][mask_prop],x_log=True)
  #ax[0].set_title("Hexagon binning")
  #ax[0].set_ylabel('elong')
  #ax[0].set_xlabel('len')
  #ax[1].set_title("Log color scale")
  #ax[1].set_xlabel('len')

  #ax = graf.draw_scatter(Prop["len"][mask_prop],Prop["razon"][mask_prop],x_log=True)
  #ax[0].set_title("Hexagon binning")
  #ax[0].set_ylabel('razon')
  #ax[0].set_xlabel('len')
  #ax[1].set_title("Log color scale")
  #ax[1].set_xlabel('len')
 
  #ax = graf.draw_hist(Prop["razon"][mask_prop],num_bins=100)
  #ax.set_xlabel('razon')

  ##############################
  ########## LONGITUD ##########
  print 'Ingresar Media'
  mean = [float(x) for x in raw_input().split(',')]
  mean.sort()

  ##############################################################
  mask  = mask_prop*(Prop["len"]>mean[0]-1.0)*(Prop["len"]<mean[-1]+1.0)
  leng, flag, rper, data, MNod, longitud = \
  lee.read_cilindros(Path_cil, Name_cil, Nfile_cil, Ncil, \
  mask, norm_x=True, norm_y=True)
  longitud /= POSFACTOR
  ##############################################################

  num_bins = 5
  percent, lista_pho, lista_vel = [], [], []

  for m in mean:
 
    mask_len = (longitud>m-1.0)*(longitud<m+1.0)
    Qq = MNod[mask_len][:,0]/MNod[mask_len][:,1]

    wall = np.linspace(0.0,100.0,num_bins+1)
    q = np.percentile(Qq,wall)
    q[0] -= 1e-6; q[-1] += 1e-6

    print '----------------------'
    print m,len(Qq)
    print '----------------------'
 
    for i in xrange(num_bins):  

      q_mask = ((Qq>q[i])*(Qq<=q[i+1]))

      print '%d-quartil %f %f Num %d' % (i,q[i],q[i+1],len(data[mask_len]["pho"][q_mask]))
      print 'Media Min ',np.mean(MNod[mask_len][q_mask][:,0]),' Media Max ',np.mean(MNod[mask_len][q_mask][:,1])

      #if(i!=num_bins-1 and i!=0): continue      

      binslog = np.mean(leng[mask_len][q_mask],axis=0).ravel()
      binsper = np.mean(rper[mask_len][q_mask],axis=0).ravel()  

      name  = 'long > %.2f && long < %.2f\n' % (m-1.0, m+1.0)
      name += 'q > %.2f && q <= %.2f' % (q[i], q[i+1])

      graf.draw_out(binslog, binsper, data[mask_len]["pho"][q_mask], \
      data[mask_len]["vmean_par"][q_mask], data[mask_len]["vmean_perp"][q_mask], \
      data[mask_len]["vrms_par"][q_mask],  data[mask_len]["vrms_perp"][q_mask], \
      lista_pho, lista_vel, title=name,LOG=True, Nlevel=25)

      percent.append(i+1)

    leng     = leng[~mask_len]
    flag     = flag[~mask_len]
    rper     = rper[~mask_len]
    data     = data[~mask_len] 
    MNod     = MNod[~mask_len]
    longitud = longitud[~mask_len]

  percent   = np.unique(np.asarray(percent))
  bdim      = len(percent)
  m         = len(mean)

  lista_pho = np.asarray(lista_pho).reshape(m,bdim, order='C')
  lista_vel = np.asarray(lista_vel).reshape(m,bdim, order='C')

  Transpose=True

  if(Transpose==True):

    for i in range(bdim):

      pdf_pho = PdfPages('%.2d_%.2d_%.2d_pho_seg_q.pdf' % (m,percent[i],num_bins))
      pdf_vel = PdfPages('%.2d_%.2d_%.2d_vel_seg_q.pdf' % (m,percent[i],num_bins))

      for j in range(m):

        pdf_pho.savefig(lista_pho[j,i])
        pdf_vel.savefig(lista_vel[j,i])
        plt.close(lista_pho[j,i])
        plt.close(lista_vel[j,i])

      pdf_pho.close()
      pdf_vel.close()

  else:

    for i in range(m):

      pdf_pho = PdfPages('%.2f_%.2d_%.2d_pho_seg_long.pdf' % (mean[i],m,num_bins))
      pdf_vel = PdfPages('%.2f_%.2d_%.2d_vel_seg_long.pdf' % (mean[i],m,num_bins))

      for j in range(bdim):

        pdf_pho.savefig(lista_pho[i,j])
        pdf_vel.savefig(lista_vel[i,j])
        plt.close(lista_pho[i,j])
        plt.close(lista_vel[i,j])

      pdf_pho.close()
      pdf_vel.close()
