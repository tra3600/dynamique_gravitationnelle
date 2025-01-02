def euler(f, y0, z0, tmin, tmax, n):
    """
    Applique le schéma d'Euler explicite pour résoudre l'équation différentielle y''(t) = f(y(t)).
    
    Paramètres:
    f (function): La fonction f(y(t)) définissant l'équation différentielle.
    y0 (float): La valeur initiale de y(t) à tmin.
    z0 (float): La valeur initiale de y'(t) à tmin.
    tmin (float): La borne inférieure de l'intervalle de temps.
    tmax (float): La borne supérieure de l'intervalle de temps.
    n (int): Le nombre de points de discrétisation.
    
    Retourne:
    tuple: Deux listes correspondant aux valeurs approchées de y(t) et z(t) aux points de discrétisation.
    """
    h = (tmax - tmin) / (n - 1)
    t_values = [tmin + i * h for i in range(n)]
    y_values = [y0]
    z_values = [z0]
    
    for i in range(1, n):
        y_prev = y_values[-1]
        z_prev = z_values[-1]
        
        y_next = y_prev + h * z_prev
        z_next = z_prev + h * f(y_prev)
        
        y_values.append(y_next)
        z_values.append(z_next)
    
    return y_values, z_values

# Exemple de test
import math

def f(y):
    omega = 2 * math.pi
    return -omega**2 * y

y0 = 3
z0 = 0
tmin = 0
tmax = 3
n = 100
y_values, z_values = euler(f, y0, z0, tmin, tmax, n)
print(y_values)
print(z_values)