import numpy as np
import matplotlib.pyplot as plt

ciltype = np.dtype([
        ("npart", np.int32),
        ("pho",   np.float32),
        ("vmean",  np.float32),
        ("vrms",   np.float32),
        ("vmean_par",  np.float32),
        ("vrms_par",   np.float32),
        ("vmean_per",  np.float32),
        ("vrms_per",   np.float32),
                  ])

def init(name):

  contenido = []

  archivo = open(name,"r")
  for linea in archivo.readlines():
    contenido.append(linea.strip('\n'))
  archivo.close()

  Path_prop = contenido[0]
  Name_prop = contenido[1]
  Path_cil  = contenido[2]
  Name_cil  = contenido[3]
  Nfile_cil = np.int32(contenido[4])
  Ncil      = np.int32(contenido[5])
  Path_par  = contenido[6]
  Name_par  = contenido[7]
  Nfile_par = np.int32(contenido[8])

  return Path_prop,Name_prop, \
  Path_cil,Name_cil,          \
  Nfile_cil,Ncil,             \
  Path_par,Name_par,Nfile_par

def some_read(pack):

  POSFACTOR = pack[0]
  filename  = pack[1]
  Ncil      = pack[2]
  extend    = pack[3]
  mask      = pack[4]

  aux_k, aux_flag, aux_data, aux_rper, aux_leng, aux_MNod, aux_lll \
  = [], [], [], [], [], [], []

  binario = open(filename,"rb")
  N = np.fromfile(binario,dtype=np.int32, count=1)[0]

  if extend == True:
  
    ext_filename = filename.replace("densidad","extend")
    ext_binario = open(ext_filename,"rb")
    ext_n = np.fromfile(ext_binario,dtype=np.int32, count=1)[0]
    assert(N==ext_n)

  print "Num",N

  for _ in range(N):

    idfil = np.fromfile(binario,dtype=np.int32,count=1)[0]

    fflag = np.fromfile(binario,dtype=np.int32,count=1)[0]
    longitud = np.fromfile(binario,dtype=np.float32,count=1)[0]/POSFACTOR
    M0 = np.fromfile(binario,dtype=np.float32,count=1)[0]
    M1 = np.fromfile(binario,dtype=np.float32,count=1)[0]
    lenbins = np.fromfile(binario,dtype=np.int32,count=1)[0]

    dist = np.fromfile(binario,dtype=np.float32,count=lenbins)/POSFACTOR

    rr, sub = [None]*Ncil, [None]*Ncil
    for j in range(Ncil):
      rr[j] = np.fromfile(binario,dtype=np.float32,count=1)[0]/POSFACTOR
      sub[j] = np.fromfile(binario,dtype=ciltype,count=lenbins)

    ###########################################################################
    if extend == True:
      ext_idfil = np.fromfile(ext_binario,dtype=np.int32,count=1)[0]
      assert(idfil==ext_idfil)
      ext_longitud = np.fromfile(ext_binario,dtype=np.float32,count=1)[0]/POSFACTOR
      assert(longitud==ext_longitud)
      ext_lenbins = np.fromfile(ext_binario,dtype=np.int32,count=1)[0]
      assert(ext_lenbins%2==0)
      ext_lenbins = ext_lenbins/2

      ########## INICIO #########
      ext_dist = np.fromfile(ext_binario,dtype=np.float32,count=ext_lenbins)/POSFACTOR
      
      ext_rr, ext_sub = [None]*Ncil, [None]*Ncil
      for j in range(Ncil):
        ext_rr[j]  = np.fromfile(ext_binario,dtype=np.float32,count=1)[0]/POSFACTOR
        ext_sub[j] = np.fromfile(ext_binario,dtype=ciltype,count=ext_lenbins)

      dist = np.concatenate((ext_dist[:-1],dist))
      rr   = np.mean((ext_rr,rr),axis=0)
      for j in range(Ncil):
        sub[j] = np.concatenate((ext_sub[j][:-1],sub[j]),axis=0)
      ###########################

      ########### FINAL #########
      ext_dist = np.fromfile(ext_binario,dtype=np.float32,count=ext_lenbins)/POSFACTOR
      
      ext_rr, ext_sub = [None]*Ncil, [None]*Ncil
      for j in range(Ncil):
        ext_rr[j]  = np.fromfile(ext_binario,dtype=np.float32,count=1)[0]/POSFACTOR
        ext_sub[j] = np.fromfile(ext_binario,dtype=ciltype,count=ext_lenbins)

      dist = np.concatenate((dist,ext_dist[1:]))
      rr   = np.mean((rr,ext_rr),axis=0)
      for j in range(Ncil):
        sub[j] = np.concatenate((sub[j],ext_sub[j][1:]),axis=0)
      ###########################

    if sum([longitud>=item[0] and longitud<=item[1] \
    for item in mask]) == 0: continue

    aux_flag.append(fflag)
    aux_leng.append(dist)
    aux_rper.append(rr)
    aux_data.append(sub)    
    aux_MNod.append([M0, M1])
    aux_lll.append(longitud)

    for j in range(Ncil):
      if(np.any(sub[j]["npart"]==0)):
        aux_k.append(idfil)
        break

  return [N, np.asarray(aux_leng), np.asarray(aux_flag), np.asarray(aux_rper), \
  np.asarray(aux_data), np.asarray(aux_MNod), np.asarray(aux_lll), np.asarray(aux_k)]


def read_cilindros(POSFACTOR, Path, Name, Nfile, Ncil, mask, extend=False):

  from multiprocessing import Pool

  Total = 0

  mask = [[item-1.0, item+1.0] for item in mask]

  args = []
  for ifile in range(Nfile):
    filename = "%s%s.%.2d.bin" % (Path,Name,ifile)
    args.append([POSFACTOR, filename, Ncil, extend, mask])

  pool = Pool()
  result = pool.map(some_read, args)
  pool.close() 

  flag, data, rper, leng, MNod, lll \
  = [], [], [], [], [], []

  while result:
    Total += result[-1][0]
    leng.append(result[-1][1])
    flag.append(result[-1][2])
    rper.append(result[-1][3])
    data.append(result[-1][4])
    MNod.append(result[-1][5])
    lll.append( result[-1][6])
    #k.append(result[-1][7])
    result.pop()

  leng = np.concatenate(leng)
  flag = np.concatenate(flag)
  rper = np.concatenate(rper)
  data = np.concatenate(data)
  MNod = np.concatenate(MNod)
  lll  = np.concatenate(lll)

  print "Total  ",Total
  #print "Quedan ",len(data)
  #print "Caidos ",len(k)

  return leng,flag,rper,data,MNod,lll


#################################################################

proptype = np.dtype([
        ("flag",  np.int32),
        ("size",  np.int32),
        ("razon", np.float32),
        ("len",   np.float32),
        ("elong", np.float32),
        ("rms",   np.float32),
                  ])

def read_prop(Path,Name,POSFACTOR=1.0):

  filename = "%s%s" % (Path,Name)
  binario = open(filename,"rb")
  N = np.fromfile(binario, dtype=np.int32, count=1)[0]
  print "Num Segmentos",N
  
  data = np.fromfile(binario, dtype=proptype)
  data["len"] = data["len"]/POSFACTOR

  return data
  
