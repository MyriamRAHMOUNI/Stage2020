import sys 
import pandas as pd
import os

#inputs
mon_fichier = os.path.basename(sys.argv[1])
milgenome_fichier = os.path.basename(sys.argv[2])
print('Je cherche '+mon_fichier+' dans '+ milgenome_fichier)

#merging
mf=pd.read_csv(sys.argv[1], sep=" ")
gf=pd.read_csv(sys.argv[2], sep=" ")
merged=gf.merge(mf,how='right')
del mf
del gf
res = merged.iloc[:,[2,3]]
del merged
res = res.dropna()
list_rs = res.iloc[:,[0]]
#output
nom_resultat= mon_fichier[:-15]+"rs_pval"
list_rs_resultat= mon_fichier[:-15]+"list_rs"
list_rs.to_csv(list_rs_resultat, sep = " ",index=False, header=False)
res.to_csv(nom_resultat, sep = " ",index=False)
