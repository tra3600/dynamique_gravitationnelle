def smul(scalar, lst):
    """
    Multiplie chaque élément de la liste `lst` par `scalar` et renvoie une nouvelle liste.
    
    Paramètres:
    scalar (int ou float): Le nombre par lequel multiplier chaque élément de la liste.
    lst (list of int ou float): La liste de nombres à multiplier.
    
    Retourne:
    list of int ou float: Une nouvelle liste avec chaque élément multiplié par `scalar`.
    """
    return [scalar * x for x in lst]

# Exemple de test
print(smul(2, [1, 2, 3]))  # [2, 4, 6]

def vsom(lst1, lst2):
    """
    Calcule la somme terme à terme de deux listes de nombres de même longueur.
    
    Paramètres:
    lst1 (list of int ou float): La première liste de nombres.
    lst2 (list of int ou float): La deuxième liste de nombres.
    
    Retourne:
    list of int ou float: Une nouvelle liste avec la somme terme à terme des deux listes.
    """
    return [x + y for x, y in zip(lst1, lst2)]

# Exemple de test
print(vsom([1, 2, 3], [4, 5, 6]))  # [5, 7, 9]

def vdif(lst1, lst2):
    """
    Calcule la différence terme à terme de deux listes de nombres de même longueur (lst1 - lst2).
    
    Paramètres:
    lst1 (list of int ou float): La première liste de nombres.
    lst2 (list of int ou float): La deuxième liste de nombres.
    
    Retourne:
    list of int ou float: Une nouvelle liste avec la différence terme à terme des deux listes.
    """
    return [x - y for x, y in zip(lst1, lst2)]

# Exemple de test
print(vdif([1, 2, 3], [4, 5, 6]))  # [-3, -3, -3]