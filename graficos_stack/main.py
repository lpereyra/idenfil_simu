import numpy as np
import matplotlib.pyplot as plt
import lector as lee
import draw as graf
import sys

if __name__ == '__main__':

  Path_cil  = "../"
  Path_snap = "/mnt/mirta3/marioagustin/1600_400Mpc/snapshots/"
  Name_snap = "snapshot_104"
  Name_cil_01  =  "104_84588_type_01_smooth_matrix_vcm_extend_30_10_10.00_Mpc_LOG_0.79_0.17"
  Name_cil_02  =  "104_84588_type_02_smooth_matrix_vcm_extend_30_10_10.00_Mpc_LOG_0.79_0.17"

  assert(np.sum([ord(a) != ord(b) for a,b in \
  zip(Name_cil_01.replace('type_01',''),Name_cil_02.replace('type_02',''))]) == 0)

  fixed, extend, smooth, vc, lin, flag_vel, rperp = \
  lee.set_variables(Name_cil_02) 

  lbox, Mpart, npart = lee.read_datos(128,Path_snap,Name_snap,1.0)

  POSFACTOR = 1000. #estan en Kpc

  ####### BINES LONGITUD #######
  print 'Ingresar Media'
  mean = [float(x) for x in raw_input().split(',')]
  mean.sort()

  num_bins = 5
  delta_m = 10.0
  mask_read = [[mean[i]-delta_m, mean[i]+delta_m] for i in xrange(len(mean))]

  ##############################################################
  ##############################################################
  ##############################################################

  #flag, longitud, MNod, data, m_middle = \
  #lee.read_new(POSFACTOR, Path_cil, Name_cil_02, mask_read)

  ##############################################################
  ##############################################################
  ##############################################################
  ##############################################################
  
  flag_02, longitud_02, MNod_02, data_02, m_middle_02 = \
  lee.read_new(POSFACTOR, Path_cil, Name_cil_02, mask_read)

  flag_01, longitud_01, MNod_01, data_01, m_middle_01 = \
  lee.read_new(POSFACTOR, Path_cil, Name_cil_01, mask_read)

  mask = (MNod_01[:,0]/MNod_01[:,1])>0.2

  flag_01     = flag_01[mask]
  longitud_01 = longitud_01[mask]
  MNod_01     = MNod_01[mask]
  data_01     = data_01[mask] 
  m_middle_01 = m_middle_01[mask]
  
  del mask

  print "rper %d" % (rperp)

  flag     = np.concatenate((flag_01,flag_02))
  longitud = np.concatenate((longitud_01, longitud_02)) 
  MNod     = np.concatenate((MNod_01, MNod_02)) 
  data     = np.concatenate((data_01,data_02))
  m_middle = np.concatenate((m_middle_01, m_middle_02)) 

  del flag_01, flag_02
  del longitud_01, longitud_02
  del MNod_01, MNod_02
  del data_01, data_02
  del m_middle_01, m_middle_02

  ##############################################################
  ##############################################################
  ##############################################################
  ##############################################################

  title = "FIL"
  lista_pho, lista_vel, lista_std, lista_gr = graf.realizo_graficos(num_bins,mean,delta_m,rperp,data,MNod,longitud,flag,m_middle,extend,lin,flag_vel,LOG=True,name_input=title)

  percent = np.arange(num_bins+1)+1

  title = "fil_long"
  graf.save_pdf(lista_pho,lista_vel,lista_std,lista_gr,title,num_bins+1,mean,percent,smooth,extend,fixed,vc,lin,flag_vel)

  del lista_pho, lista_vel, lista_std, lista_gr
