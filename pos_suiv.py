def force2(m1, p1, m2, p2):
    G = 6.67e-11  # constante de gravitation universelle en N·m²/kg²
    r = [p2[i] - p1[i] for i in range(3)]
    r_magnitude = math.sqrt(sum([r[i]**2 for i in range(3)]))
    force_magnitude = G * m1 * m2 / r_magnitude**3
    force = [force_magnitude * r[i] for i in range(3)]
    return force

def forceN(j, m, pos):
    N = len(m)
    force_total = [0.0, 0.0, 0.0]
    for k in range(N):
        if k != j:
            force = force2(m[j], pos[j], m[k], pos[k])
            force_total = [force_total[i] + force[i] for i in range(3)]
    return force_total

def pos_suiv(m, pos, vit, h):
    """
    Calcule les positions des N corps à l'instant t_{i+1} en utilisant le schéma de Verlet.
    
    Paramètres:
    m (list of float): Liste des masses des N corps (en kilogrammes).
    pos (list of list of float): Liste des positions des N corps à l'instant t_i (en mètres).
    vit (list of list of float): Liste des vitesses des N corps à l'instant t_i (en mètres par seconde).
    h (float): Pas d'intégration (en secondes).
    
    Retourne:
    list of list of float: Liste des positions des N corps à l'instant t_{i+1}.
    """
    N = len(m)
    new_pos = []

    for j in range(N):
        F_j = forceN(j, m, pos)
        new_pos_j = [
            pos[j][k] + h * vit[j][k] + 0.5 * h**2 * F_j[k] / m[j]
            for k in range(3)
        ]
        new_pos.append(new_pos_j)
    
    return new_pos

# Exemple de test
masses = [5.972e24, 7.348e22]  # masses de la Terre et de la Lune en kg
positions = [[0, 0, 0], [384400000, 0, 0]]  # positions initiales en m
vitesses = [[0, 0, 0], [0, 1022, 0]]  # vitesses initiales en m/s (orbite circulaire simplifiée)
h = 60  # pas d'intégration en secondes
print(pos_suiv(masses, positions, vitesses, h))