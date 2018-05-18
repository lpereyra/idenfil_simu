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

  if 'fixed' in Name_cil:
    if 'fixed' in Name_par:
      fixed = True
    else:
      print "SE ROMPE MAL LOS NOMBRES"
      assert(0)
  else: 
    if 'fixed' not in Name_par:
      fixed = False
    else:
      print "SE ROMPE MAL LOS NOMBRES"
      assert(0)

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

  par_nod, fil_nod = [], []
  Save_Seg ,Save_Pares = False, False
  ##############################################################

  print 'Guardar Segmentos?'
  select = raw_input() 

  if('s' in select or 'S' in select or 'y' in select or 'Y' in select):

    Save_Seg = True

    leng, flag, rper, data, MNod, longitud = \
    lee.read_cilindros(POSFACTOR, Path_cil,Name_cil,Nfile_cil,Ncil,mean,extend=True)

    title = "SEG"
    if "sin_vel" in Name_cil: title += "_sin_vel"
    if "media" in Name_cil:   title += "_media"
    lista_pho, lista_vel = graf.realizo_graficos(num_bins, mean, leng, flag, rper, \
    data, MNod, longitud, title, LOG = True)

    percent = np.arange(num_bins)+1

    title = "seg_q"
    if "sin_vel" in Name_cil: title += "_sin_vel"
    if "media" in Name_cil:   title += "_media"
    graf.save_pdf(lista_pho,lista_vel,title,num_bins,mean,percent,fixed)

    title = "seg_long"
    if "sin_vel" in Name_cil: title += "_sin_vel"
    if "media" in Name_cil:   title += "_media"
    graf.save_pdf(lista_pho,lista_vel,title,num_bins,mean,percent,fixed)

    for m in mean:
      mask_len = (longitud>m-1.0)*(longitud<m+1.0)*(flag==2) # TENGO QUE INCLUIR EL FLAG
      fil_nod.append(np.concatenate((MNod[mask_len][:,0],MNod[mask_len][:,1])))

    del lista_pho,lista_vel
    ##############################################################

  print 'Guardar Pares?'
  select = raw_input()
 
  if('s' in select or 'S' in select or 'y' in select or 'Y' in select):

    Save_Pares = True
    ##############################################################
    leng, flag, rper, data, MNod, longitud = \
    lee.read_cilindros(POSFACTOR, Path_par, Name_par, Nfile_par, Ncil, mean,extend=True)
 
    title = "PAR"
    if "sin_vel" in Name_par: title += "_sin_vel"
    if "media" in Name_par:   title += "_media"
    lista_pho, lista_vel = graf.realizo_graficos(num_bins, mean, leng, flag, rper, \
    data, MNod, longitud, title, LOG = True)

    percent = np.arange(num_bins)+1
    title = "pares_q"
    if "sin_vel" in Name_par: title += "_sin_vel"
    if "media" in Name_par:   title += "_media"
    graf.save_pdf(lista_pho,lista_vel,title,num_bins,mean,percent,fixed)

    title = "pares_long"
    if "sin_vel" in Name_par: title += "_sin_vel"
    if "media" in Name_par:   title += "_media"
    graf.save_pdf(lista_pho,lista_vel,title,num_bins,mean,percent,fixed)

    for m in mean:
      mask_len = (longitud>m-1.0)*(longitud<m+1.0)
      par_nod.append(np.concatenate((MNod[mask_len][:,0],MNod[mask_len][:,1])))

    del lista_pho,lista_vel
    ##############################################################


  if Save_Pares*Save_Seg == True:

    LOG = True

    for i in range(len(mean)):

      mmax = np.max([np.max(fil_nod[i]), np.max(par_nod[i])])
      mmin = np.min([np.min(fil_nod[i]), np.min(par_nod[i])])

      if LOG == True:
        fil_nod[i] = np.log10(fil_nod[i])
        par_nod[i] = np.log10(par_nod[i])
        mmax = np.max([np.max(fil_nod[i]), np.max(par_nod[i])])
        mmin = np.min([np.min(fil_nod[i]), np.min(par_nod[i])])
        bins = np.logspace(np.log10(mmin),np.log10(mmax),num=20)
      else:
        mmax = np.max([np.max(fil_nod[i]), np.max(par_nod[i])])
        mmin = np.min([np.min(fil_nod[i]), np.min(par_nod[i])])
        bins = np.linspace(mmin,mmax,num=20)

      f = plt.figure(figsize=(11,8))
      f.suptitle('len %f\n' % mean[i],x=0.47)
      f.subplots_adjust(wspace=0.5,hspace=0.3)
      
      Nmedia = len(fil_nod[i])/2

      ax = plt.subplot(3, 1, 1)
      ax.hist(par_nod[i][:Nmedia],bins,color="r",histtype='step',normed=True)
      ax.hist(fil_nod[i][:Nmedia],bins,color="b",histtype='step',normed=True)
      ax.set_title('NODO menos masivo')

      ax = plt.subplot(3, 1, 2)
      ax.hist(par_nod[i][Nmedia:],bins,color="r",histtype='step',normed=True)
      ax.hist(fil_nod[i][Nmedia:],bins,color="b",histtype='step',normed=True)
      ax.set_title('NODO mas masivo')

      ax = plt.subplot(3, 1, 3)
      ax.hist(par_nod[i],bins,color="r",histtype='step',normed=True)
      ax.hist(fil_nod[i],bins,color="b",histtype='step',normed=True)
      ax.set_title('DOS NODOS')


      plt.show()


