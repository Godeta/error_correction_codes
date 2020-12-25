# Installation :
# python -m pip install bchlib

import bchlib
import hashlib
import os
import random

# on créer un objet BCH, le polynome est un polynome primitif du champ de Galois de 2 puissance 13,14 ou 15
# https://link.springer.com/content/pdf/bbm%3A978-1-4615-1509-8%2F1.pdf
BCH_POLYNOME = 16659 
# nombre de bits que l'on va changer pour mettre des erreurs
BCH_BITS = 50
bch = bchlib.BCH(BCH_POLYNOME, BCH_BITS)

# données aléatoires
data = bytearray(os.urandom(512))

# encodage dans un paquet
ecc = bch.encode(data)
paquet = data + ecc

# affichage du hash du paquet
sha1_initial = hashlib.sha1(paquet)
print('sha1: %s' % (sha1_initial.hexdigest(),))

"""Fonction pour ajouter aléatoirement des erreurs dans notre paquet"""
def bitflip(paquet):
    byte_num = random.randint(0, len(paquet) - 1)
    bit_num = random.randint(0, 7)
    paquet[byte_num] ^= (1 << bit_num)

# on rajoute des erreurs 
for _ in range(BCH_BITS):
    bitflip(paquet)

# affichage du hash du paquet
sha1_corrupt = hashlib.sha1(paquet)
print('sha1: %s' % (sha1_corrupt.hexdigest(),))

# de-paquetize
data, ecc = paquet[:-bch.ecc_bytes], paquet[-bch.ecc_bytes:]

# correction
bitflips = bch.decode_inplace(data, ecc)
print('bitflips: %d' % (bitflips))

# paquetize
paquet = data + ecc

# affichage du hash du paquet
sha1_corrected = hashlib.sha1(paquet)
print('sha1: %s' % (sha1_corrected.hexdigest(),))

if sha1_initial.digest() == sha1_corrected.digest():
    print('Correction effectuee !')
else:
    print('Echec')
