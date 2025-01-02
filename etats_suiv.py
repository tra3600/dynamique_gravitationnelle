def etat_suiv(m, pos, vit, h):
    """
    Calcule les positions et vitesses des N corps à l'instant t_{i+1} en utilisant le schéma de Verlet.
    
    Paramètres:
    m (list of float): Liste des masses des N corps (en kilogrammes).
    pos (list of list of float): Liste des positions des N corps à l'instant t_i (en mètres).
    vit (list of list of float): Liste des vitesses des N corps à l'instant t_i (en mètres par seconde).
    h (float): Pas d'intégration (en secondes).
    
    Retourne:
    tuple: Deux listes - (positions des N corps à l'instant t_{i+1}, vitesses des N corps à l'instant t_{i+1}).
    """
    N = len(m)
    new_pos = pos_suiv(m, pos, vit, h)
    new_vit = []

    for j in range(N):
        F_j = forceN(j, m, pos)
        F_j_next = forceN(j, m, new_pos)
        new_vit_j = [
            vit[j][k] + 0.5 * h * (F_j[k] + F_j_next[k]) / m[j]
            for k in range(3)
        ]
        new_vit.append(new_vit_j)
    
    return new_pos, new_vit

# Exemple de test
print(etat_suiv(masses, positions, vitesses, h))