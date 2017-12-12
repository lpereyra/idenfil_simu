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
    lee.read_cilindros(Path_cil,Name_cil,Nfile_cil,Ncil,norm_x=True,norm_y=True)

    longitud /= POSFACTOR
    mask     = (longitud>mean[0]-1.0)*(longitud<mean[-1]+1.0)*(flag==2)
    leng     = leng[mask]
    flag     = flag[mask]
    rper     = rper[mask]
    data     = data[mask] 
    MNod     = MNod[mask]
    longitud = longitud[mask]

    data["pho"] += 1.0 ### PARA USAR DELTA+1.0

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

      data_mask = np.mean(data["vmean_par"][mask_len],axis=0)
      min_max_velpar[0] = min_max_velpar[0] if np.min(data_mask)>min_max_velpar[0] else np.min(data_mask)
      min_max_velpar[1] = min_max_velpar[1] if np.max(data_mask)<min_max_velpar[1] else np.max(data_mask)

      data_mask = np.mean(data["vmean_per"][mask_len],axis=0)
      min_max_velper[0] = min_max_velper[0] if np.min(data_mask)>min_max_velper[0] else np.min(data_mask)
      min_max_velper[1] = min_max_velper[1] if np.max(data_mask)<min_max_velper[1] else np.max(data_mask)

      data_mask = np.mean(data["vrms_par"][mask_len],axis=0)
      min_max_rmspar[0] = min_max_rmspar[0] if np.min(data_mask)>min_max_rmspar[0] else np.min(data_mask)
      min_max_rmspar[1] = min_max_rmspar[1] if np.max(data_mask)<min_max_rmspar[1] else np.max(data_mask)

      data_mask = np.mean(data["vrms_per"][mask_len],axis=0)
      min_max_rmsper[0] = min_max_rmsper[0] if np.min(data_mask)>min_max_rmsper[0] else np.min(data_mask)
      min_max_rmsper[1] = min_max_rmsper[1] if np.max(data_mask)<min_max_rmsper[1] else np.max(data_mask)

    #min_max_pho    = [0.0,   5000.]
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

      q = np.linspace(0.0,100.0,num_bins+1)
      q = np.percentile(Qq,q); q[0] -= 1e-6; q[-1] += 1e-6

      print '----------------------\n',m,len(Qq),'\n----------------------'
 
      for i in xrange(num_bins):

        q_mask = ((Qq>q[i])*(Qq<=q[i+1]))

        print '%d-quartil %f %f Num %d' % (i,q[i],q[i+1],len(data[mask_len]["pho"][q_mask]))
        print 'Media Min ',np.mean(MNod[mask_len][q_mask][:,0]),' Media Max ',np.mean(MNod[mask_len][q_mask][:,1])

        binslog = np.mean(leng[mask_len][q_mask],axis=0).ravel()
        binsper = np.mean(rper[mask_len][q_mask],axis=0).ravel()  

        name  = 'long > %.2f && long < %.2f\n' % (m-1.0, m+1.0)
        name += 'SEG\n'
        name += 'q > %.2f && q <= %.2f' % (q[i], q[i+1])

        graf.draw_out(binslog, binsper, \
        data[mask_len]["pho"][q_mask],           min_max_pho, \
        data[mask_len]["vmean_par"][q_mask],  min_max_velpar, \
        data[mask_len]["vmean_per"][q_mask], min_max_velper, \
        data[mask_len]["vrms_par"][q_mask],   min_max_rmspar, \
        data[mask_len]["vrms_per"][q_mask],  min_max_rmsper, \
        lista_pho, lista_vel, title=name,LOG=True, Nlevel=25)

        del q_mask

      leng     = leng[~mask_len]
      rper     = rper[~mask_len]
      data     = data[~mask_len] 
      MNod     = MNod[~mask_len]
      longitud = longitud[~mask_len]
      
      del mask_len,Qq

    del leng,rper,data,MNod,longitud

    percent = np.arange(num_bins)+1
    graf.save_pdf(lista_pho,lista_vel,'seg_q',num_bins,mean,percent)
    graf.save_pdf(lista_pho,lista_vel,'seg_long',num_bins,mean,percent)

    del lista_pho,lista_vel
    ##############################################################

  print 'Guardar Pares?'
  select = raw_input()
  
  if('s' in select or 'S' in select or 'y' in select or 'Y' in select):

    ###############################################################################################################################
    #
    # GUARDO LOS PARES
    #
    ###############################################################################################################################
    
    ##############################################################
    leng, flag, rper, data, MNod, longitud = \
    lee.read_cilindros(Path_par,Name_par,Nfile_par,Ncil,norm_x=True,norm_y=True)

    longitud /= POSFACTOR
    mask     = (longitud>mean[0]-1.0)*(longitud<mean[-1]+1.0)*(flag==2)
    leng     = leng[mask]
    flag     = flag[mask]
    rper     = rper[mask]
    data     = data[mask] 
    MNod     = MNod[mask]
    longitud = longitud[mask]

    data["pho"] += 1.0 ### PARA USAR DELTA+1.0

    del flag
    del mask
    ##############################################################

    ##############################################################
    lista_pho, lista_vel = [], []

    for m in mean:
 
      mask_len = (longitud>m-1.0)*(longitud<m+1.0)
      Qq = MNod[mask_len][:,0]/MNod[mask_len][:,1]

      q = np.linspace(0.0,100.0,num_bins+1)
      q = np.percentile(Qq,q); q[0] -= 1e-6; q[-1] += 1e-6

      print '----------------------\n',m,len(Qq),'\n----------------------'
 
      for i in xrange(num_bins):

        q_mask = ((Qq>q[i])*(Qq<=q[i+1]))

        print '%d-quartil %f %f Num %d' % (i,q[i],q[i+1],len(data[mask_len]["pho"][q_mask]))
        print 'Media Min ',np.mean(MNod[mask_len][q_mask][:,0]),' Media Max ',np.mean(MNod[mask_len][q_mask][:,1])

        binslog = np.mean(leng[mask_len][q_mask],axis=0).ravel()
        binsper = np.mean(rper[mask_len][q_mask],axis=0).ravel()  

        name  = 'long > %.2f && long < %.2f\n' % (m-1.0, m+1.0)
        name += 'PAR\n'
        name += 'q > %.2f && q <= %.2f' % (q[i], q[i+1])

        graf.draw_out(binslog, binsper, \
        data[mask_len]["pho"][q_mask],           min_max_pho, \
        data[mask_len]["vmean_par"][q_mask],  min_max_velpar, \
        data[mask_len]["vmean_per"][q_mask], min_max_velper, \
        data[mask_len]["vrms_par"][q_mask],   min_max_rmspar, \
        data[mask_len]["vrms_per"][q_mask],  min_max_rmsper, \
        lista_pho, lista_vel, title=name,LOG=True, Nlevel=25)

        del q_mask

      leng     = leng[~mask_len]
      rper     = rper[~mask_len]
      data     = data[~mask_len] 
      MNod     = MNod[~mask_len]
      longitud = longitud[~mask_len]
      
      del mask_len,Qq

    del leng,rper,data,MNod,longitud

    percent = np.arange(num_bins)+1
    graf.save_pdf(lista_pho,lista_vel,'pares_q',num_bins,mean,percent)
    graf.save_pdf(lista_pho,lista_vel,'pares_long',num_bins,mean,percent)

    del lista_pho,lista_vel
    ##############################################################

