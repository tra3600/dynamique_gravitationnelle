import math

G = 6.67e-11  # constante de gravitation universelle en N·m²/kg²

def force2(m1, p1, m2, p2):
    """
    Calcule la force gravitationnelle exercée par le corps 2 sur le corps 1.
    
    Paramètres:
    m1 (float): Masse du corps 1 en kilogrammes.
    p1 (list of float): Position du corps 1 sous forme [x1, y1, z1] en mètres.
    m2 (float): Masse du corps 2 en kilogrammes.
    p2 (list of float): Position du corps 2 sous forme [x2, y2, z2] en mètres.
    
    Retourne:
    list of float: Force exercée par le corps 2 sur le corps 1 sous forme [Fx, Fy, Fz] en newtons.
    """
    # Calculer le vecteur position de p2 par rapport à p1
    r = [p2[i] - p1[i] for i in range(3)]
    # Calculer la distance entre les deux corps
    r_magnitude = math.sqrt(sum([r[i]**2 for i in range(3)]))
    # Calculer la force gravitationnelle
    force_magnitude = G * m1 * m2 / r_magnitude**3
    # Calculer les composantes de la force
    force = [force_magnitude * r[i] for i in range(3)]
    return force

# Exemple de test
m1 = 5.972e24  # masse de la Terre en kg
p1 = [0, 0, 0]  # position de la Terre en m
m2 = 7.348e22  # masse de la Lune en kg
p2 = [384400000, 0, 0]  # position de la Lune en m (384,400 km de la Terre)
print(force2(m1, p1, m2, p2))