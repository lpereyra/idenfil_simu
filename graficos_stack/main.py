import numpy as np
import lector as lee
import draw as graf
import matplotlib.pyplot as plt
import sys

if __name__ == '__main__':

  Path_prop, Name_prop, \
  Path_cil, Name_cil, Nfile_cil, Ncil,  \
  Path_par, Name_par, Nfile_par \
  = lee.init(sys.argv[1])

  POSFACTOR = 1000.  
  print 'Numero de cilindros ',Ncil

  Prop = lee.read_prop(Path_prop,Name_prop,POSFACTOR)

  mask_prop = (Prop["flag"]>=2)
  print "Mask Segmentos",sum(mask_prop)
  #mask_prop = (Prop["size"]>2)*mask_prop
  #print "Segmentos mask N>2",len(Prop["size"][mask_prop])

  Scatter = False
  if(Scatter==True):

    graf.draw_scatter(Prop["len"][mask_prop],Prop["elong"][mask_prop],'len','elong',x_log=True)
    graf.draw_scatter(Prop["len"][mask_prop],Prop["razon"][mask_prop],'len','razon',x_log=True)

  Hist    = False
  if(Hist==True):
    graf.draw_hist(Prop["razon"][mask_prop],'razon',num_bins=100)

  del Prop

  ####### BINES LONGITUD #######
  print 'Ingresar Media'
  mean = [float(x) for x in raw_input().split(',')]
  mean.sort()

  num_bins = 5
  ##############################

  ##############################################################
  print 'Guardar Segmentos?'
  select = raw_input()
  
  if('s' in select or 'S' in select or 'y' in select or 'Y' in select):

    leng, flag, rper, data, MNod, longitud = \
    lee.read_cilindros(Path_cil,Name_cil,Nfile_cil,Ncil,extend=True)

    lista_pho, lista_vel = graf.realizo_graficos(POSFACTOR, num_bins, mean, leng, flag, rper, \
    data, MNod, longitud, 'SEG', LOG = False)

    percent = np.arange(num_bins)+1
    graf.save_pdf(lista_pho,lista_vel,'seg_q',num_bins,mean,percent)
    graf.save_pdf(lista_pho,lista_vel,'seg_long',num_bins,mean,percent)

    longitud /= POSFACTOR
    mask_len = (longitud>mean[0]-1.0)*(longitud<mean[0]+1.0)
    fil_nod = np.log10(np.concatenate((MNod[mask_len][:,0],MNod[mask_len][:,1])))

    del lista_pho,lista_vel
    ##############################################################

  print 'Guardar Pares?'
  select = raw_input()
  
  if('s' in select or 'S' in select or 'y' in select or 'Y' in select):

    ##############################################################
    leng, flag, rper, data, MNod, longitud = \
    lee.read_cilindros(Path_par,Name_par,Nfile_par,Ncil,extend=True)

    lista_pho, lista_vel = graf.realizo_graficos(POSFACTOR, num_bins, mean, leng, flag, rper, \
    data, MNod, longitud, 'PAR', LOG = True)

    percent = np.arange(num_bins)+1
    graf.save_pdf(lista_pho,lista_vel,'pares_q',num_bins,mean,percent)
    graf.save_pdf(lista_pho,lista_vel,'pares_long',num_bins,mean,percent)

    longitud /= POSFACTOR
    mask_len = (longitud>mean[0]-1.0)*(longitud<mean[0]+1.0)
    par_nod = np.log10(np.concatenate((MNod[mask_len][:,0],MNod[mask_len][:,1])))

    del lista_pho,lista_vel
    ##############################################################

  #bins = np.logspace(np.log10(par_nod.min()),np.log10(par_nod.max()),num=20)
  #f, ax = plt.subplots(1)
  #ax.hist(par_nod,bins,color="r",histtype='step',normed=True)
  #ax.hist(fil_nod,bins,color="b",histtype='step',normed=True)
  #plt.show()


