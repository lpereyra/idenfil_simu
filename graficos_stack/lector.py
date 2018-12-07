import numpy as np
import matplotlib.pyplot as plt

ciltype = np.dtype([
        ("npart", np.int32),
        ("pho",   np.float32),
        #("vmean",  np.float32),
        #("vrms",   np.float32),
        ("vmean_par",  np.float32),
        ("vrms_par",   np.float32),
        ("vmean_per",  np.float32),
        ("vrms_per",   np.float32),
                  ])

propiedadestype = np.dtype([
    ("flag"       , np.int32), 
    ("npart"      , np.int32), 
    ("Mass",      (np.float32,2)),
    ("Vnodos",    (np.float32,6)),
    ("q"          , np.float32), 
    ("long"       , np.float32), 
    ("elong"      , np.float32), 
    ("rms"        , np.float32), 
    ])

stdgrup = np.dtype([
  ("Save",      np.uint32),
  ("Id",        np.uint32),
  ("Pos",      (np.float32,3)),
  ("NumPart",   np.uint32)
  ])

io_header = np.dtype([
            ("npart",      (np.uint32,6)),
            ("mass",     (np.float64,6)),
            ("time",         np.float64),
            ("redshift",     np.float64),
            ("flag_sfr",       np.uint32),
            ("flag_feedback",  np.uint32),
            ("npartTotal", (np.uint32,6)),
            ("flag_cooling",   np.uint32),
            ("num_files",      np.uint32),
            ("BoxSize",      np.float64),
            ("Omega0",       np.float64),              
            ("OmegaLambda",  np.float64),
            ("HubbleParam",  np.float64),
            ("fill", (np.character,256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8))
            ])

stdset = np.dtype([
  ("r",     np.float32),
  ("mass",  np.float32)
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
  Path_par  = contenido[5]
  Name_par  = contenido[6]
  Nfile_par = np.int32(contenido[7])

  return Path_prop,Name_prop, \
  Path_cil,Name_cil,          \
  Nfile_cil,             \
  Path_par,Name_par,Nfile_par

def set_variables(Name):

  if 'fixed' in Name:
    fixed = True
  else: 
    fixed = False

  if 'extend' in Name:
    extend = True
  else: 
    extend = False

  if 'smooth' in Name:
    smooth = True
  else: 
    smooth = False

  if 'vc' in Name:
    vc = True
  else: 
    vc = False

  if 'LIN' in Name:
    lin = True
  else: 
    lin = False

  if 'vcm' in Name:
    flag_vel = 'vcm'
  else: 
    flag_vel = 'crudo'

  rperp = np.asarray(Name.split('_'))
  rperp = rperp[np.where('Mpc' == rperp)[0]-1]
  rperp = np.float32(rperp)[0]

  return fixed, extend, smooth, vc, lin, flag_vel, rperp

def read_new(POSFACTOR, Path, Name, mask):

  filename = "%s%s.bin" % (Path,Name)

  binario = open(filename,"rb")
  N  = np.fromfile(binario,dtype=np.int32, count=1)[0]
  nx = np.fromfile(binario,dtype=np.int32, count=1)[0]
  ny = np.fromfile(binario,dtype=np.int32, count=1)[0]
  lenbins = nx*ny

  aux_flag, aux_data, aux_leng, aux_MNod, aux_lll, aux_middle \
  = [], [], [], [], [], []

  for _ in range(N):

    idfil = np.fromfile(binario,dtype=np.int32,count=1)[0]
    fflag = np.fromfile(binario,dtype=np.int32,count=1)[0]
    longitud = np.fromfile(binario,dtype=np.float32,count=1)[0]/POSFACTOR
    M0 = np.fromfile(binario,dtype=np.float32,count=1)[0]
    M1 = np.fromfile(binario,dtype=np.float32,count=1)[0]
    m_middle = np.fromfile(binario,dtype=np.float32,count=1)[0]
    sub = np.fromfile(binario,dtype=ciltype,count=lenbins)

    if sum([longitud>=item[0] and longitud<=item[1] \
    for item in mask]) == 0: continue

    aux_flag.append(fflag)
    aux_lll.append(longitud)
    aux_MNod.append([M0, M1])
    aux_data.append(sub.reshape(ny,nx,order='C'))
    aux_middle.append(m_middle)

  return np.asarray(aux_flag), np.asarray(aux_lll), np.asarray(aux_MNod), np.asarray(aux_data), np.asarray(aux_middle)

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

  if(Nfile == 1):
    filename = "%s%s.bin" % (Path,Name)
    args.append([POSFACTOR, filename, Ncil, extend, mask])
  else:
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

def leeheader(name,POSFACTOR):

  binario = open(name,"rb")

  d1 = np.fromfile(binario,dtype=np.int32, count=1)[0]
  header = np.fromfile(binario,dtype=io_header, count=1)[0]
  d2 = np.fromfile(binario,dtype=np.int32, count=1)[0]
  assert(d1==d2)

  binario.close()

  lbox     = header["BoxSize"]/POSFACTOR
  npart    = np.sum(header["npartTotal"])
  redshift = header["redshift"]
  omegam   = header["Omega0"]
  omegal   = header["OmegaLambda"]
  hparam   = header["HubbleParam"]
  Mpart    = header["mass"][1]

  return lbox, npart, redshift, omegam, omegal, hparam, Mpart

def read_datos(Nfile, Path_snap, Name_snap, POSFACTOR):

  if Nfile > 1:
    filename = "%s%s.0" % (Path_snap,Name_snap)
  else:
    filename = "%s%s" % (Path_snap,Name_snap)

  lbox, npart, redshift, omegam, omegal, \
  hparam, Mpart = leeheader(filename,POSFACTOR)

  print "***********************************"
  print "*   Parametros de la simulacion   *"
  print "***********************************"
  print "  Numero de particulas = %u" % (npart)
  print "  Lado del box = %g"         % (lbox)
  print "  Redshift = %g"             % (redshift)
  print "  Omega Materia = %g"        % (omegam)
  print "  Omega Lambda = %g"         % (omegal)
  print "  Parametro de Hubble = %g"  % (hparam)
  print "  Masa por particula = %g"   % (Mpart)
  #print "  Softening = %g",           %  (soft)
  print "***********************************"
  print "***********************************"

  lbox /= POSFACTOR

  return lbox, Mpart, npart 

def read_seg(Path_seg,Name_seg,POSFACTOR):

  filename = "%s%s" % (Path_seg,Name_seg)
  print "Reading file: %s" % (filename)

  binario = open(filename,"rb")
  N = np.fromfile(binario, dtype=np.int32, count=1)[0]

  segment = []
  for i in range(N):
    chain = np.fromfile(binario, dtype=np.int32, count=1)[0]
    idp = np.fromfile(binario, dtype=np.int32, count=chain)
    segment.append(idp)

  print "%d segmentos" % (N)

  # REEMPLAZA LA PALABRA PARA BUSCAR EL ARCHIVO #
  filename = filename.replace('segmentos','propiedades')

  print "Reading file: %s" % (filename)

  binario = open(filename,"rb")
  N = np.fromfile(binario, dtype=np.int32, count=1)[0]
  prop = np.fromfile(binario, dtype=propiedadestype)

  print "Lee propiedades de %d segmentos" % (N)

  segment = np.asarray(segment) # transforma
  #prop["rms"]  /= prop["long"]
  prop["long"] /= POSFACTOR

  return segment,prop

def read_seg_smooth(Path_seg,Name_seg,POSFACTOR):

  filename = "%s%s" % (Path_seg,Name_seg)
  print "Reading file: %s" % (filename)

  binario = open(filename,"rb")
  N = np.fromfile(binario, dtype=np.int32, count=1)[0]

  segment = []
  for i in range(N):
    chain = np.fromfile(binario, dtype=np.int32, count=1)[0]
    Pos = np.fromfile(binario, dtype=np.float32, count=3*chain).ravel()
    Pos = np.reshape(Pos/POSFACTOR,(chain,3)) 
    segment.append(Pos)

  print "%d Segmentos" % (N)

  # REEMPLAZA LA PALABRA PARA BUSCAR EL ARCHIVO #
  filename = filename.replace('segmentos','propiedades')

  print "Reading file: %s" % (filename)

  binario = open(filename,"rb")
  N = np.fromfile(binario, dtype=np.int32, count=1)[0]
  prop = np.fromfile(binario, dtype=propiedadestype)

  print "Lee propiedades de %d Segmentos" % (N)

  segment = np.asarray(segment) # transforma
  #prop["rms"]  /= prop["long"]
  prop["long"] /= POSFACTOR

  return segment,prop

def read_grup(Path_gr,Name_gr,POSFACTOR):

  filename = "%s%s" % (Path_gr,Name_gr)
  binario = open(filename,"rb")

  print "Reading file: %s" % (filename)

  N = np.fromfile(binario,dtype=np.int32, count=1)[0]
  data = np.fromfile(binario,dtype=stdgrup, count=N).ravel()
  binario.close()

  data["Pos"] /= POSFACTOR

  return data

def crea_set_inter(Seg,PSeg,Gr,lbox,Mpart,mean):

  mask = (PSeg["flag"]==2)
  Seg, PSeg = Seg[mask], PSeg[mask]

  set_inter = []

  for m in mean:
 
    mask = ((PSeg["long"]>(m-1.0)) * (PSeg["long"]<(m+1.0)))
 
    for i in xrange(len(Seg[mask])):

      fil = np.copy(Seg[mask][i][:-1])

      dx = Gr[fil]["Pos"][:,2]   
      dy = Gr[fil]["Pos"][:,1]
      dz = Gr[fil]["Pos"][:,0]

      fil = fil[1:]

      aux = np.empty(len(fil),dtype=stdset)

      dx = np.diff(dx)
      dy = np.diff(dy)
      dz = np.diff(dz)

      dx[dx>   0.5*lbox] -= lbox
      dx[dx<= -0.5*lbox] += lbox

      dy[dy>   0.5*lbox] -= lbox
      dy[dy<= -0.5*lbox] += lbox

      dz[dz>   0.5*lbox] -= lbox
      dz[dz<= -0.5*lbox] += lbox

      aux["r"] = np.cumsum(np.sqrt(dx*dx+dy*dy+dz*dz))
      aux["r"] /= PSeg[mask][i]["long"]
      aux["mass"] = Mpart*Gr[fil]["NumPart"]
      
      set_inter.append(aux)

    Seg, PSeg = Seg[~mask], PSeg[~mask]    

  return np.asarray(set_inter)
