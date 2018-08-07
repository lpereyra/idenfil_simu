import numpy as np

################################

stdpart = np.dtype([
          ("Pos",    (np.float32,3)),
          ("Id",      np.uint32)
          ])

################################

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

################################

stdgrup = np.dtype([
  ("Save",      np.uint32),
  ("Id",        np.uint32),
  ("Pos",      (np.float32,3)),
  ("NumPart",   np.uint32)
  ])

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

def read_gadget(Nfile, Path_snap, Name_snap, POSFACTOR, sli):

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

  ind = 0

  #P = np.empty(npart,dtype=stdpart)
  P = []

  for ifile in range(Nfile):

    if(Nfile>1):
      filename = "%s%s.%d" % (Path_snap,Name_snap,ifile)
    else:
      filename = "%s%s" % (Path_snap,Name_snap)

    #P, ind = lee(filename,P,ind,POSFACTOR,sli)
    P, ind = lee(filename,P,ind,POSFACTOR,sli)

  #assert(ind==npart)

  P = np.concatenate(P)

  return P, lbox, Mpart 

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

################################

def lee(filename,Q,ind,POSFACTOR,sli):

  binario = open(filename,"rb")

  print "Reading file: %s" % (filename)

  d1 = np.fromfile(binario,dtype=np.int32, count=1)[0]
  header = np.fromfile(binario,dtype=io_header, count=1)[0]
  d2 = np.fromfile(binario,dtype=np.int32, count=1)[0]
  assert(d1==d2)

  for k in xrange(6):
    if(k == 1): #ONLY KEEP DARK MATTER PARTICLES
      n = header["npart"][k]

  P = np.empty(n,dtype=stdpart)

  d1 = np.fromfile(binario,dtype=np.int32, count=1)[0]
  for k in xrange(6):
   n = header["npart"][k]
   r = np.fromfile(binario,dtype=np.float32,count=3*n).ravel()
   if(k == 1): #ONLY KEEP DARK MATTER PARTICLES
     r = np.reshape(r*POSFACTOR,(n,3)) 
     mask = (r[:,2]>sli[0]) * (r[:,2]<sli[1])
     r = r[mask]
     n = sum(mask)
     #Q[ind:ind+n]["Pos"] = r
     P[0:n]["Pos"] = r
     pc = n
  d2 = np.fromfile(binario,dtype=np.int32, count=1)[0]
  assert(d1==d2)

  # SALTA LAS VELOCIDADES #
  d1 = np.fromfile(binario,dtype=np.int32, count=1)[0]
  binario.seek(d1,1)
  d2 = np.fromfile(binario,dtype=np.int32, count=1)[0]
  assert(d1==d2)

  d1 = np.fromfile(binario,dtype=np.int32, count=1)[0]
  for k in xrange(6):
   n = header["npart"][k]
   r = np.fromfile(binario,dtype=np.int32, count=n).ravel()
   if(k == 1): #ONLY KEEP DARK MATTER PARTICLES
     r = r[mask]
     n = sum(mask)
     #Q[ind:ind+n]["Id"] = r
     P[0:n]["Id"] = r
     pc = n
  d2 = np.fromfile(binario,dtype=np.int32, count=1)[0]
  assert(d1==d2)
  
  binario.close()

  Q.append(P)
  
  del P
  #return Q, ind+pc
  return Q, sum(mask)

################################

def read_grup(Path_gr,Name_gr,POSFACTOR):

  filename = "%s%s" % (Path_gr,Name_gr)
  binario = open(filename,"rb")

  print "Reading file: %s" % (filename)

  N = np.fromfile(binario,dtype=np.int32, count=1)[0]
  data = np.fromfile(binario,dtype=stdgrup, count=N).ravel()
  binario.close()

  data["Pos"] /= POSFACTOR

  return data

################################

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

  return lbox, Mpart 

