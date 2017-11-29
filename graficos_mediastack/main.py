import numpy as np
import lector as lee
import draw as graf
import matplotlib.pyplot as plt
import sys

if __name__ == '__main__':

  Path_prop, Name_prop, \
  Path_cil, Name_cil, \
  Nfile_cil, Ncil \
  = lee.init(sys.argv[1])

  POSFACTOR = 1000.  
  print 'Numero de cilindros ',Ncil

  Prop = lee.read_prop(Path_prop,Name_prop,POSFACTOR)
  
  mask = (Prop["flag"]>=2)
  print "Mask Segmentos",sum(mask)
  #mask = (Prop["size"]>2)*mask
  #print "Segmentos mask N>2",len(Prop["size"][mask])

  plt.ion()

  #ax = graf.draw_scatter(Prop["len"][mask],Prop["elong"][mask],x_log=True)
  #ax[0].set_title("Hexagon binning")
  #ax[0].set_ylabel('elong')
  #ax[0].set_xlabel('len')
  #ax[1].set_title("Log color scale")
  #ax[1].set_xlabel('len')

  #ax = graf.draw_scatter(Prop["len"][mask],Prop["razon"][mask],x_log=True)
  #ax[0].set_title("Hexagon binning")
  #ax[0].set_ylabel('razon')
  #ax[0].set_xlabel('len')
  #ax[1].set_title("Log color scale")
  #ax[1].set_xlabel('len')
 
  #ax = graf.draw_hist(Prop["razon"][mask],num_bins=100)
  #ax.set_xlabel('razon')

  ##############################
  ########## LONGITUD ##########
  print 'Ingresar Media'
  mean = [float(x) for x in raw_input().split(',')]

  for m in mean:

    mask *= (Prop["len"]>m-1.0)
    mask *= (Prop["len"]<m+1.0)

    ##############################################################
    leng, flag, rper, data, MNod = \
    lee.read_cilindros(Path_cil, Name_cil, Nfile_cil, Ncil, \
    mask, norm_x=True, norm_y=True)

    num_bins = 5
    wall = np.linspace(0.0,100.0,num_bins+1)
    q = np.percentile(Prop[mask]["razon"],wall)
    q[0] -= 1e-6; q[-1] += 1e-6

    Total = len(Prop[mask]); contador = 0
    for i in xrange(num_bins):  

      mask = ((MNod>q[i])*(MNod<=q[i+1]))
      contador += np.sum(mask)

      if(i!=num_bins-1 and i!=0): continue      

      binslog = np.mean(leng[mask],axis=0).ravel()
      binsper = np.mean(rper[mask],axis=0).ravel()  

      name  = 'long > %.2f && long < %.2f\n' % (m-1.0, m+1.0)
      name += 'q > %.2f && q <= %.2f' % (q[i], q[i+1])

      graf.draw_out(binslog, binsper, data["pho"][mask], \
      data["vmean_par"][mask], data["vmean_perp"][mask], \
       data["vrms_par"][mask],  data["vrms_perp"][mask], \
      title=name,LOG=True, Nlevel=25)

      print '%d %s' % (len(data["pho"][mask]),name.replace('\n',' '))

    assert(contador==Total)
