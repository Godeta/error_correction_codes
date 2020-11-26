# Entrer les données à transmettre
# Calculer le nombre de bits redondants à ajouter
# Déterminer la parité des bits
# Créer des erreurs dans les données pour tester
# Vérifier les erreurs

def hammingCodeTest():
    "Test un exemple d'utilisation de Hamming code"
    # Entrer les données à transmettre 
    data = '1011001'
    
    # Calculer le nombre de bits redondants à ajouter 
    m = len(data) 
    r = calcRedundantBits(m) 
    
    # Déterminer les positions des bits redondants 
    arr = posRedundantBits(data, r) 
    
    # Déterminer la parité des bits 
    arr = calcParityBits(arr, r) 
    
    # Les données à transférer
    print("Les donnees transferees sont " + arr)   
    
    # Simuler une erreur dans la transmission en changeant une valeur binaire
    # 10101001110 -> 11101001110, erreur en dixième position. 
    
    arr = '10111001110'
    print("Les donnees erronees sont " + arr) 
    # detectError renvoie la position en partant du dernier bit, de façon plus intuitive pour renvoyer la position du bit en prenant le premier pour 1 
    # j'ai mis la longueur totale des données je sosutrait la position du bit en partant de la fin et j'ajoute 1 (car sinon première position =0)
    correction = len(arr)-detectError(arr, r) +1 
    print("La position de l'erreur est " + str(correction))
    print ("Le bit incorrect est donc "+arr[correction])
    
def calcRedundantBits(m): 
    "Calcule le nombre de bits redondants, autrement dit les bits qui ne portent pas de données mais servent à vérifier la parité"
    # Utilisation de la formule 2 ^ r >= m + r + 1 
    # Iteration de 0 à m et retourne la valeur 
  
    for i in range(m): 
        if(2**i >= m + i + 1): 
            return i 
  
  
def posRedundantBits(data, r): 
    "Renvoie un tableau contenant la position des bits redondants"
    # Redundancy bits are placed at the positions 
    # which correspond to the power of 2. 
    j = 0
    k = 1
    m = len(data) 
    res = '' 
  
    # If position is power of 2 then insert '0' 
    # Else append the data 
    for i in range(1, m + r+1): 
        if(i == 2**j): 
            res = res + '0'
            j += 1
        else: 
            res = res + data[-1 * k] 
            k += 1
  
    # The result is reversed since positions are 
    # counted backwards. (m + r+1 ... 1) 
    return res[::-1] 
  
  
def calcParityBits(arr, r): 
    "Renvoie les données binaires avec les bits de parités ajoutés"
    n = len(arr) 
  
    # For finding rth parity bit, iterate over 
    # 0 to r - 1 
    for i in range(r): 
        val = 0
        for j in range(1, n + 1): 
  
            # If position has 1 in ith significant 
            # position then Bitwise OR the array value 
            # to find parity bit value. 
            if(j & (2**i) == (2**i)): 
                val = val ^ int(arr[-1 * j]) 
                # -1 * j is given since array is reversed 
  
        # String Concatenation 
        # (0 to n - 2^r) + parity bit + (n - 2^r + 1 to n) 
        arr = arr[:n-(2**i)] + str(val) + arr[n-(2**i)+1:] 
    return arr 
  
  
def detectError(arr, nr): 
    "Détecte si il y a une erreur dans des données binaires"
    n = len(arr) 
    res = 0
  
    # Calculate parity bits again 
    for i in range(nr): 
        val = 0
        for j in range(1, n + 1): 
            if(j & (2**i) == (2**i)): 
                val = val ^ int(arr[-1 * j]) 
  
        # Create a binary no by appending 
        # parity bits together. 
  
        res = res + val*(10**i) 
  
    # Convert binary to decimal 
    return int(str(res), 2) 
  
# code executé
hammingCodeTest()