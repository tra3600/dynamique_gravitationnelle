
def forceN(j, m, pos):
    """
    Calcule la force totale exercée sur le corps j par tous les autres corps du système.
    
    Paramètres:
    j (int): Indice du corps considéré.
    m (list of float): Liste des masses des N corps.
    pos (list of list of float): Liste des positions des N corps, chaque position étant une liste [x, y, z].
    
    Retourne:
    list of float: Force totale exercée sur le corps j sous forme [Fx, Fy, Fz] en newtons.
    """
    N = len(m)
    force_total = [0.0, 0.0, 0.0]
    
    for k in range(N):
        if k != j:
            force = force2(m[j], pos[j], m[k], pos[k])
            force_total = [force_total[i] + force[i] for i in range(3)]
    
    return force_total

# Exemple de test
masses = [5.972e24, 7.348e22]  # masses de la Terre et de la Lune en kg
positions = [[0, 0, 0], [384400000, 0, 0]]  # positions de la Terre et de la Lune en m
print(forceN(0, masses, positions))
print(forceN(1, masses, positions))