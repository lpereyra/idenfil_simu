import numpy as np

ciltype = np.dtype([
        ("npart", np.int32),
        ("pho",   np.float32),
        ("vmean",  np.float32),
        ("vrms",   np.float32),
        ("vmean_par",  np.float32),
        ("vrms_par",   np.float32),
        ("vmean_perp",  np.float32),
        ("vrms_perp",   np.float32),
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
 
  return Path_prop,Name_prop,\
  Path_cil,Name_cil,\
  Nfile_cil,Ncil

def read_cilindros(Path, Name, Nfile, Ncil, mask, norm_x=False, norm_y=False):

  Total = 0
  k, flag, data, rper, leng, MNod = [], [], [], [], [], []

  if(norm_x==True): print "x normed"
  if(norm_y==True): print "y normed" 

  for ifile in range(Nfile):

    filename = "%s%s.%.2d.bin" % (Path,Name,ifile)
    print filename
    binario = open(filename,"rb")
    N = np.fromfile(binario,dtype=np.int32, count=1)[0]
    Total+=N
    print "Num",N

    for _ in range(N):

      idfil = np.fromfile(binario,dtype=np.int32,count=1)[0]
      fflag = np.fromfile(binario,dtype=np.int32,count=1)[0]
      longitud = np.fromfile(binario,dtype=np.float32,count=1)[0]
      M0 = np.fromfile(binario,dtype=np.float32,count=1)[0]
      M1 = np.fromfile(binario,dtype=np.float32,count=1)[0]
      lenbins = np.fromfile(binario,dtype=np.int32,count=1)[0]

      dist = []
      for j in range(lenbins):
        dist.append(np.fromfile(binario,dtype=np.float32,count=1))

      rr   = [None]*Ncil
      sub  = [None]*Ncil
      for j in range(Ncil):
        rr[j]  = np.fromfile(binario,dtype=np.float32,count=1)[0]        
        sub[j] = np.fromfile(binario,dtype=ciltype,count=lenbins)

      if(mask[idfil]==False): continue

      if(norm_x==True): dist /= longitud
      if(norm_y==True): rr   /= longitud*0.5

      flag.append(fflag)
      leng.append(dist)
      rper.append(rr)
      data.append(sub)    
      MNod.append(np.array([M0,M1]))

      for j in range(Ncil):
        if(np.any(sub[j]["npart"]==0)):
          k.append(idfil)
          break

  leng = np.asarray(leng)
  flag = np.asarray(flag)
  rper = np.asarray(rper)
  data = np.asarray(data)
  MNod = np.asarray(MNod)

  print "Total  ",Total
  print "Quedan ",len(data)
  print "Caidos ",len(k)

  return leng,flag,rper,data,MNod


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
  
