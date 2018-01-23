fichier = open("fichier.txt" , "r" )
fichier2= open("list_ECK.txt" , "r")
fichier3= open("list_gene.txt" , "w")
lignes= fichier.readlines()
lignes2=fichier2.readlines()
mydict={}

for uneligne in lignes:
    champs=uneligne.split('\t')
    mydict[champs[0]]=champs[1]
#print mydict

for l in lignes2:
	eck = l[:-1]
	if eck == "":
		continue
	if eck in mydict:
		print mydict[eck]
		fichier3.write(mydict[eck]+'\n')
	else:
		
		print 'eck pas dans dict', eck
	
"""
for valeur in  mydict.values():
	for ligne in lignes2:
		champs2=ligne.split()
			if mydict[champs[0]] == champs2:
				print mydict.values()
"""
