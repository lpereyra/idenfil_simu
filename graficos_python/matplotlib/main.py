import numpy as np
import carga_datos as dat
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import sys

if __name__ == '__main__':

  colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')  

  Nsnap, Path_snap, Name_snap, \
  Path_gr, Name_gr, \
  Path_seg, Name_seg \
  = dat.init(sys.argv[1])

  print Nsnap, Path_snap, Name_snap

  sli = [0,30.0] #in Mpc

  #P, lbox, Mpart = dat.read_gadget(Nsnap,Path_snap,Name_snap,POSFACTOR,sli)
  #plt.scatter(P["Pos"][:,0], P["Pos"][:,1], alpha=0.5, s=0.5, c='b', edgecolors=None)

  POSFACTOR = 1. #estan en Mpc
  lbox, Mpart = dat.read_datos(Nsnap,Path_snap,Name_snap,POSFACTOR)

  mass_cut  = 2000.  # Mcut
  POSFACTOR = 1000. # estan en Kpc
  Gr = dat.read_grup(Path_gr,Name_gr,POSFACTOR)
  mask = (Gr["Pos"][:,2]>sli[0]) * (Gr["Pos"][:,2]<sli[1])
  mask_num = mask * (Mpart*Gr["NumPart"]>mass_cut)
  
  m = Gr["Pos"][:,0]>lbox
  Gr["Pos"][:,0][m] = Gr["Pos"][:,0][m]-lbox
  m = Gr["Pos"][:,1]>lbox
  Gr["Pos"][:,1][m] = Gr["Pos"][:,1][m]-lbox
  m = Gr["Pos"][:,2]>lbox
  Gr["Pos"][:,2][m] = Gr["Pos"][:,2][m]-lbox

  m = Gr["Pos"][:,0]<0
  Gr["Pos"][:,2][m] = Gr["Pos"][:,2][m]+lbox
  m = Gr["Pos"][:,1]<0
  Gr["Pos"][:,1][m] = Gr["Pos"][:,1][m]+lbox
  m = Gr["Pos"][:,2]<0
  Gr["Pos"][:,0][m] = Gr["Pos"][:,0][m]+lbox

  #Tree = cKDTree(Gr["Pos"][mask],boxsize=lbox) 
  #counts = Tree.query_ball_point(Gr["Pos"][mask], 5.0)	# vecinos dentro de r Mpc
  #counts = 1.0*np.asarray(map(len,counts))/np.max(counts)

  # slice
  plt.scatter(Gr["Pos"][mask_num][:,0],Gr["Pos"][mask_num][:,1],s=5.0,c='g', edgecolors=None)
  print "Quedan %d Grupos Mayor a %.2g" % (np.sum(mask_num),mass_cut*1.e10)
  mass_cut  = 100.  # Mcut
  mask_num = mask * (Mpart*Gr["NumPart"]>mass_cut)
  plt.scatter(Gr["Pos"][mask_num][:,0],Gr["Pos"][mask_num][:,1],s=3.0,c='b', edgecolors=None)


  del mask_num

  Seg, PSeg, R = dat.read_seg(Path_seg,Name_seg,Gr["NumPart"],POSFACTOR)

  mask_fil = (PSeg["flag"]==2)      #Tipo 2
  Seg      = Seg[mask_fil] 
  PSeg     = PSeg[mask_fil] 
  R        = R[mask_fil] 
 
  print "Quedan %d Segmentos Tipo 2" % (np.sum(mask_fil))
 
  del mask_fil

  contador = 0
  for i in xrange(len(Seg)):

    fil = Seg[i]

    if(np.all(mask[fil])==False): continue
    
    #print i, mask[fil], fil
    contador = contador + 1 
    
    zg = Gr[fil]["Pos"][:,2]   
    yg = Gr[fil]["Pos"][:,1]
    xg = Gr[fil]["Pos"][:,0]

    dx = np.diff(xg)
    dy = np.diff(yg)
    dz = np.diff(zg)
    N = len(dx)

    for r in xrange(N):

      #condicion en x
      if(dx[r]>   0.5*lbox): xg[r+1:] -= lbox
      if(dx[r]<= -0.5*lbox): xg[r+1:] += lbox

      #condicion en y
      if(dy[r]>   0.5*lbox): yg[r+1:] -= lbox
      if(dy[r]<= -0.5*lbox): yg[r+1:] += lbox

      #condicion en z
      if(dz[r]>   0.5*lbox): zg[r+1:] -= lbox
      if(dz[r]<= -0.5*lbox): zg[r+1:] += lbox

    #if(np.all(zg>5.0)==True): continue

    color = colors[i % len(colors)]
    plt.plot(xg,yg,c=color,linewidth=1)
    #plt.scatter(xg,yg,s=Gr[fil]["NumPart"],c='b')

  print "Sobreviven ",contador," filamentos"

  plt.xlim(0,lbox)
  plt.ylim(0,lbox)
  plt.show()
