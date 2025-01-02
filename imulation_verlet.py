import math

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

def verlet_step(m, pos, vit, h):
    N = len(m)
    new_pos = []
    new_vit = []

    for j in range(N):
        F_j = forceN(j, m, pos)
        new_pos_j = [
            pos[j][k] + h * vit[j][k] + 0.5 * h**2 * F_j[k] / m[j]
            for k in range(3)
        ]
        new_pos.append(new_pos_j)
    
    for j in range(N):
        F_j = forceN(j, m, pos)
        F_j_next = forceN(j, m, new_pos)
        new_vit_j = [
            vit[j][k] + 0.5 * h * (F_j[k] + F_j_next[k]) / m[j]
            for k in range(3)
        ]
        new_vit.append(new_vit_j)
    
    return new_pos, new_vit

def simulation_verlet(deltat, n):
    """
    Simule les positions des corps en utilisant le schéma d'intégration de Verlet.
    
    Paramètres:
    deltat (float): Incrément de temps en secondes.
    n (int): Nombre d'itérations.
    
    Retourne:
    list of list of list of float: Liste des positions des corps pour chaque instant t0, t0 + deltat, ..., t0 + n*deltat.
    """
    positions = [p0]
    current_pos = p0
    current_vit = v0

    for _ in range(n):
        new_pos, new_vit = verlet_step(masse, current_pos, current_vit, deltat)
        positions.append(new_pos)
        current_pos, current_vit = new_pos, new_vit
    
    return positions

# Exemple de test
t0 = 0  # date des conditions initiales
p0 = [[0, 0, 0], [1, 0, 0]]  # positions initiales en unité astronomique
v0 = [[0, 0, 0], [0, 1, 0]]  # vitesses initiales en km/s
masse = [5.972e24, 7.348e22]  # masses en kg
deltat = 60  # incrément de temps en secondes
n = 100  # nombre d'itérations

positions = simulation_verlet(deltat, n)
print(positions)