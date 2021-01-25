from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import gzip
import sys
import os
import re
import numpy as np
from operator import attrgetter
import urllib.request
import gzip
import shutil

#Lecture nom molecule et id chebi
def nom_id(fichier, numeroligne):
    ligne_prec = "> <ChEBI Name>" #la ligne qui précède le nom
    ligne_prec2 = "> <NAME>" #Peut etre ecrit de cette maniere
    ligne = fichier.readline(numeroligne)
    while ligne_prec not in ligne and ligne_prec2 not in ligne:
        if "> <ChEBI ID>" in ligne or "> <ID>" in ligne:
            idmolecule = fichier.readline(numeroligne + 1)
            idmolecule = idmolecule.strip('%s\n')
        ligne = fichier.readline(numeroligne)
        numeroligne = numeroligne + 1
    nommolecule = fichier.readline(numeroligne + 1) #On récupère le nom
    numeroligne += 1
    result = [idmolecule, nommolecule, str(numeroligne)]
    return result

#Calcul formule de la molécule
def formule_molecule(molecule):
    formule = rdMolDescriptors.CalcMolFormula(molecule)
    formule2 = re.findall(r'([A-Z][a-z]?)(\d*)', formule)
            
    formuledecomp = ""
    for q in formule2:
        if q[0] != "H": #suppression des hydrogènes
            formuledecomp += q[0] + "\n"
            if q[1] == "":
                formuledecomp += "1\n"
            else:
                formuledecomp += q[1] + "\n"
    return formuledecomp

#Calcul degré de chaque atome de la molécule (tableau d)
def degre_atomes(molecule):
    degres = ""
    atomes = molecule.GetAtoms()
    for a in atomes:
        degres += str(a.GetDegree()) + " "
    return degres

#Calcul indice de début de chaque atome dans e (construction de v)
def indice_e(molecule):
    atomes = molecule.GetAtoms()
    indices = "\n0 "
    cpt = 0
    ret = 0
    if molecule.GetNumAtoms() > 1: #on ne traite pas les voisins si l'atome est seul
        for a in atomes:
            if cpt == 0:
                indices += str(a.GetDegree()) + " "
                ret = a.GetDegree()
            elif cpt > 0 and cpt < molecule.GetNumAtoms()-1:
                indices += str(ret + a.GetDegree()) + " "
                ret += a.GetDegree()
            cpt += 1
                
    indices += "\n"
    return indices

#Recherche des voisins de chaque atome de la molécule (construction de e)
def voisins_atomes(molecule):
    atomes = molecule.GetAtoms()
    cptE = 0 #compteur taille de e
    e = ""
    taille_e = ""
    
    for a in atomes:
        voisins = a.GetNeighbors()
        for v in voisins:
            e += str(v.GetIdx()) + " "
            cptE += 1
    taille_e = str(cptE)
    result = [taille_e, e]
    return result
    
#Nombre de types d'atomes différents
def types_atomes(molecule):
    typesatomes = [] #Tableau des types d'atome deja trouvés
    typesatomesnonperiodique = ["R", "R1", "R2", "R3", "R4","R5","R6","R7","R8","R9","R10","R11","R12","R13","R14","R15","R16","R17","R18","R19","R20","R21","R22","R23","R24","R25","R26","R27","R28","R29","R30","R31","R32","*","R#","hv","X","A","D","ACP","Ps","Mu","Mu-","0","e","?"]
    atomes = molecule.GetAtoms()
    cpttype = 0 #Compteur de types d'atomes
    types = "" #Nombre de types differents
    for a in atomes:
        if a.GetSymbol() not in typesatomes and a.GetSymbol() not in typesatomesnonperiodique:
            typesatomes.append(a.GetSymbol())
            cpttype += 1
    types = "\n" + str(cpttype)
    return types

#Coloration des atomes (lab et ptn)
def coloration_atomes(molecule):          
    coloration = "" #lab
    colorationatome = "" #ptn      
    atomes = molecule.GetAtoms()
    listeatomes = []
    for g in atomes:
        listeatomes.append(g.GetSymbol())
    listecouleurs = np.unique(listeatomes)
            
    for h in listecouleurs:
        colorationliste = []
        for v in atomes:
            if h == v.GetSymbol():
                colorationatome += str(v.GetIdx()) + " "
                colorationliste.append("1")
        if len(atomes) > 1:
            del colorationliste[len(colorationliste)-1]
            colorationliste.append("0")
        for j in colorationliste:
            coloration += j + " "
    result = [colorationatome, coloration]
    return result

#Téléchargement d'une molécule spécifique depuis chebi (format sdf)
def telechargement_molecule(idmol):
    loc = "MoleculesSDF/"+idmol + ".sdf" 
    url = "https://www.ebi.ac.uk/chebi/saveStructure.do?sdf=true&chebiId=" + idmol + "&imageId=0" #le lien de téléchargement A COMPLETER
    req = urllib.request.urlretrieve(url, loc)

#Téléchargement, extraction du fichier des 50000 molécules 3 stars de Chebi (format sdf)
def telechargement_50000():
    loc = "chebi.sdf.gz" #le chemin absolu vers ton fichier en local
    url = "ftp://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete_3star.sdf.gz" #le lien de téléchargement
    req = urllib.request.urlretrieve(url, loc)

    #on décompresse le fichier téléchargé avec ça
    new_file = "ChEBI_complete_3star.sdf" #le chemin absolu vers ton fichier décompressé en local
    with gzip.open(loc, 'rb') as file_in:
        with open(new_file, 'wb') as file_out:
            shutil.copyfileobj(file_in, file_out)

#Créer un dossier pour contenir les molécules au format personnalisé s'il n'existe pas
def creation_dossier_mol():
    try:
       os.mkdir('Molecules')
    except OSError:
        if not os.path.isdir('Molecules'):
            Raise
#Créer un dossier pour contenir les molécules au format SDF s'il n'existe pas
def creation_dossier_molSDF():
    try:
        os.mkdir('MoleculesSDF')
    except OSError:
        if not os.path.isdir('MoleculesSDF'):
            Raise

#Fonction principale de parsing qui recueille les différentes informations importantes des molécules
#pour notre fichier personnalisé
def parsing(mode):
    i = j = 1 #Compteur de molécules
    infoerror = "Molécules non lues: \n"
    cpterror = 0
    numeroligne = 0 #compteur de lignes pour retrouver les noms de molécules et id chebi
    
    creation_dossier_mol()
    creation_dossier_molSDF()
    
    if mode == -1: #Parsing sur toutes les molécules 3 stars
        telechargement_50000()
        suppl = Chem.SDMolSupplier('ChEBI_complete_3star.sdf', sanitize = False) #Lecture fichier sdf et transformation en tableau de mol
        fichier = open('ChEBI_complete_3star.sdf', 'r') #Ouverture fichier sdf pour récupérer ids chebi et noms molécules
    else: #Parsing sur une molécule 3 stars
        telechargement_molecule(mode) 
        suppl = Chem.SDMolSupplier("MoleculesSDF/"+str(mode) + ".sdf", sanitize = False) #Lecture fichier sdf et transformation en tableau de mol
        fichier = open("MoleculesSDF/"+str(mode) + '.sdf', 'r') #Ouverture fichier sdf pour récupérer ids chebi et noms molécules
    
    for mol in suppl:

        #Lecture nom molecule et id chebi
        res = []
        res = nom_id(fichier, numeroligne)
        idmolecule = res[0]
        nommolecule = res[1]
        numeroligne = int(res[2])
        
        try:
            mol.UpdatePropertyCache(strict=False)
            mol = Chem.RemoveAllHs(mol, sanitize=False)
        except:
            pass
        
        if mol is not None:
            #Passage en SMILES puis en mol
            #smiles = Chem.MolToSmiles(mol)
            #m = Chem.MolFromSmiles(smiles)
                        
            #Calcul formule
            formuledecomp = formule_molecule(mol)

            name = "Molecule_" + str(i) + ".txt"
            repertoire = "Molecules/" + idmolecule[6:] + ".txt"
            #Récupération nombre atomes et liaisons
            nbatom = mol.GetNumAtoms()
            nbbond = mol.GetNumBonds()
            infos = ""

            infosDebut = "1\n" + str(nbatom) + "\n" + str(nbbond*2) + "\n" + str(nbatom) + "\n" + str(nbatom) + "\n"

            #Récupération degrés atomes (construction de d)
            infos += "\n" + degre_atomes(mol)

            #Coloration   
            res = coloration_atomes(mol)
            colorationatome = res[0] #lab
            coloration = res[1] #ptn

            #Calcul indice de début de chaque atome dans e (construction de v)
            infos += indice_e(mol)

            #Recherche des voisins de chaque atome (construction de e)
            res = voisins_atomes(mol)
            infosDebut += res[0] #Taille de e
            infos += res[1] #e

            #Nombre de types d'atomes différents
            infos += types_atomes(mol)

            #Ajout de la formule de molécule, du nom et de la coloration
            infos += "\n" + formuledecomp + nommolecule + colorationatome + "\n" + coloration

            #Ce qu'on  ecrit dans le fichier
            #infosDebut contient un bit si molécule est lue ou non, nombre d'atomes de la molécule, de liaisons et taille de d, v et e
            #infos contient ce qui suit: d,v et e, le nombre de types d'atomes différents, formule chimique, nom molécule et coloration
            infosFinal = infosDebut + infos

            #Ecriture dans le fichier
            outf = open(repertoire, 'w+')
            outf.write(infosFinal)
            outf.close()
        else:
            #Molécule pas lue
            name = "Molecule_" + str(i) + ".txt"
            repertoire = "Molecules/" + idmolecule[6:] + ".txt"
            outf = open(repertoire, 'w+')
            st = "0 \nErreur de lecture molecule " + str(i) + "\n" + idmolecule + "\n" + nommolecule
            outf.write(st)
            outf.close()
            cpterror = cpterror + 1
            infoerror += "Molécule " + str(i) + " " + idmolecule + " ->" + nommolecule + "\n"
                
        i = i + 1
    #Ecriture fichier molécules non lues
    outerror = open("Molecules/ErreurMolecule", 'w+')
    outerror.write(infoerror)
    outerror.close()
    print("Il y a {} molécules dans le fichier".format(i-1) + " dont " + str(cpterror) + " molécules non lues")
    fichier.close()
