import numpy as np
import lector as lee
import draw as graf
import sys

if __name__ == '__main__':

  Path_prop, Name_prop, \
  Path_cil, Name_cil, \
  Nfile_cil, Ncil \
  = lee.init(sys.argv[1])

  POSFACTOR = 1000.  
  
  Prop = lee.read_prop(Path_prop,Name_prop,POSFACTOR)
  
  mask = (Prop["flag"]>=2)
  print "Mask Segmentos",sum(mask)
  mask = (Prop["size"]>3)*mask
  print "Segmentos mask N>2",len(Prop["size"][mask])


  #graf.draw_scatter(Prop["len"][mask],Prop["elong"][mask],  \
  #         x_label='len',y_label='elong',x_log=True)

  #graf.draw_scatter(Prop["len"][mask],Prop["razon"][mask],  \
  #         x_label='len',y_label='razon',x_log=True)
  
  #graf.draw_hist(Prop["razon"][mask], \
  #            x_label='razon',num_bins=100)


  leng, flag, rper, data, MNod = \
  lee.read_cilindros(Path_cil, Name_cil, Nfile_cil, Ncil, mask)

  ##############################################################

  q = np.percentile(Prop["razon"][mask],50)
  mask = (Prop["razon"][mask]>q)

  ##############################################################
 
  binslog = np.mean(leng,axis=0).ravel()
  binsper = np.mean(rper,axis=0).ravel()  

  print len(data["pho"][mask])
  graf.draw_out(binslog, binsper, data["pho"][mask], \
  data["vmean_par"][mask],  data["vmean_perp"][mask], \
  data["vrms_par"][mask],  data["vrms_perp"][mask], \
  LOG=True, Nlevel=20)

  print len(data["pho"][~mask])
  graf.draw_out(binslog, binsper, data["pho"][~mask], \
  data["vmean_par"][~mask],  data["vmean_perp"][~mask], \
  data["vrms_par"][~mask],  data["vrms_perp"][~mask], \
  LOG=True, Nlevel=20)
