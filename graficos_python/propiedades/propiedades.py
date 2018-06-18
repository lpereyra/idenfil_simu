import numpy as np
import matplotlib.pyplot as plt
import sys

################################

propiedadestype = np.dtype([
    ("flag"       , np.int32), 
    ("npart"      , np.int32), 
    ("q"          , np.float32), 
    ("long"       , np.float32), 
    ("elong"      , np.float32), 
    ("rms"        , np.float32), 
    ])

################################

################################

def init(name):

  contenido = []

  archivo = open(name,"r")
  for linea in archivo.readlines():
    contenido.append(linea.strip('\n'))
  archivo.close()

  Ngad      = np.int32(contenido[0])
  Path_gad  = contenido[1]
  Name_gad  = contenido[2]
  Path_gr   = contenido[3]
  Name_gr   = contenido[4]
  Path_seg  = contenido[5]
  Name_seg  = contenido[6]

  return Ngad,Path_gad,Name_gad,\
         Path_gr,Name_gr,\
         Path_seg,Name_seg,\

################################

################################

def read_prop(Path_seg,Name_seg,POSFACTOR):

  filename = "%s%s" % (Path_seg,Name_seg)
  filename = filename.replace('segmentos','propiedades')

  print "Reading file: %s" % (filename)

  binario = open(filename,"rb")
  N = np.fromfile(binario, dtype=np.int32, count=1)[0]
  prop = np.fromfile(binario, dtype=propiedadestype)

  print "Lee propiedades de %d segmentos" % (N)

  #prop["rms"]  /= prop["long"]
  prop["long"] /= POSFACTOR

  return prop

################################

if __name__ == '__main__':

  Nsnap, Path_snap, Name_snap, \
  Path_gr, Name_gr, \
  Path_seg, Name_seg \
  = init(sys.argv[1])

  POSFACTOR = 1000. #estan en Kpc
  data = read_prop(Path_seg,Name_seg,POSFACTOR)

  num_bins = 100
  mask = (data["flag"]>=2)
  print "Mask Segmentos",sum(mask)
  mask = (data["npart"]>2)*mask
  print "Segmentos mask N>2",len(data["npart"][mask])
  
  #HISTOGRAMA LONGITUD
  log = 0
  len_min = (data["long"][mask].min())-1e-10
  len_max = (data["long"][mask].max())+1e-10
  print "len min", len_min
  print "len max", len_max
  if log==0:
    len_bins = np.linspace(len_min,len_max,num=num_bins)
  else:
    len_bins = np.logspace(np.log10(len_min),np.log10(len_max),num=num_bins)
  
  #HISTOGRAMA ELONGACION
  log = 0;
  elong_min = (data["elong"][mask].min())-1e-10
  elong_max = (data["elong"][mask].max())+1e-10
  print "elong min", elong_min
  print "elong max", elong_max
  if log==0:
    elong_bins = np.linspace(elong_min,elong_max,num=num_bins)
  else:
    elong_bins = np.logspace(np.log10(elong_min),np.log10(elong_max),num=num_bins)
  
  #HISTOGRAMA RMS
  log = 0;
  rms_min = (data["rms"][mask].min())-1e-10
  rms_max = (data["rms"][mask].max())+1e-10
  print "rms min", rms_min
  print "rms max", rms_max
  if log==0:
    rms_bins = np.linspace(rms_min,rms_max,num=num_bins)
  else:
    rms_bins = np.logspace(np.log10(rms_min),np.log10(rms_max),num=num_bins)
  
  #f, axarr = plt.subplots(3)
  f, axarr = plt.subplots(2)
  
  #HISTOGRAMA LONGITUD
  axarr[0].hist(data["long"][mask],len_bins,color="r",normed=True)
  axarr[0].set_xlabel("Longitud")
  axarr[0].set_xlim([len_min,len_max])
  
  #HISTOGRAMA ELONGACION
  axarr[1].hist(data["elong"][mask],elong_bins,color="g")
  axarr[1].set_xlabel("Elongacion")
  axarr[1].set_xlim([elong_min,elong_max])
  
  #HISTOGRAMA RMS
  #################axarr[2].hist(data["rms"][mask],rms_bins,color="b")
  #################axarr[2].set_xlabel("RMS")
  #################axarr[2].set_xlim([rms_min,rms_max])
  #################
  f.subplots_adjust(hspace=0.5)
  
  plt.show()
