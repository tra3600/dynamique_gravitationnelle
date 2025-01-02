# dynamique_gravitationnelle
modélise les interactions physiques entres planètes et corps célestes

I.A – Valeur des expressions Python
I.A.1) [1, 2, 3] + [4, 5, 6]
La valeur de cette expression est [1, 2, 3, 4, 5, 6] car l'opérateur + concatène les deux listes.

I.A.2) 2 * [1, 2, 3]
La valeur de cette expression est [1, 2, 3, 1, 2, 3] car l'opérateur * répète la liste deux fois.

Résumé
Valeur des expressions Python :

[1, 2, 3] + [4, 5, 6] → [1, 2, 3, 4, 5, 6]
2 * [1, 2, 3] → [1, 2, 3, 1, 2, 3]
Fonction smul :

Multiplie chaque élément d'une liste par un nombre donné.
Fonctions vsom et vdif :

vsom : Calcule la somme terme à terme de deux listes.
vdif : Calcule la différence terme à terme de deux listes (lst1 - lst2).

II.A.1) Mise en forme de l'équation différentielle
Pour montrer que l’équation différentielle du second ordre (\forall t \in I, y''(t) = f(y(t))) peut être mise sous la forme d’un système différentiel du premier ordre en (z(t)) et (y(t)), nous allons introduire une nouvelle fonction intermédiaire (z(t)).

Nous définissons (z(t)) comme la dérivée de (y(t)) :

z
(
t
)
=
y
′
(
t
)
Ainsi, la dérivée de (z(t)) est :

z
′
(
t
)
=
y
″
(
t
)
En utilisant l'équation différentielle (\forall t \in I, y''(t) = f(y(t))), nous avons :

z
′
(
t
)
=
f
(
y
(
t
)
)
Nous obtenons donc le système différentiel du premier ordre suivant :

{
y
′
(
t
)
=
z
(
t
)
z
′
(
t
)
=
f
(
y
(
t
)
)
Ce système est noté ((S)).

II.A.2) Discrétisation du problème
Soit (n) un entier strictement supérieur à 1 et (\mathbb{J} = [[0, n - 1]]). Nous posons :

h
=
t
m
a
x
−
t
m
i
n
n
−
1
et

∀
i
∈
J
,
 
t
i
=
t
m
i
n
+
i
h
Nous devons montrer que pour tout entier (i \in [[0, n - 2]]) :

y
(
t
i
+
1
)
=
y
(
t
i
)
+
∫
t
i
t
i
+
1
z
(
t
)
d
t
et

z
(
t
i
+
1
)
=
z
(
t
i
)
+
∫
t
i
t
i
+
1
f
(
y
(
t
)
)
d
t
Preuve de la première équation
En utilisant la définition de (y(t)) et la règle de la trapèze pour l'intégrale, nous avons :

y
(
t
i
+
1
)
=
y
(
t
i
)
+
∫
t
i
t
i
+
1
y
′
(
t
)
d
t
Mais nous savons que (y'(t) = z(t)), donc :

y
(
t
i
+
1
)
=
y
(
t
i
)
+
∫
t
i
t
i
+
1
z
(
t
)
d
t
Preuve de la deuxième équation
De manière similaire, en utilisant la définition de (z(t)) et la règle de la trapèze pour l'intégrale, nous avons :

z
(
t
i
+
1
)
=
z
(
t
i
)
+
∫
t
i
t
i
+
1
z
′
(
t
)
d
t
Mais nous savons que (z'(t) = f(y(t))), donc :

z
(
t
i
+
1
)
=
z
(
t
i
)
+
∫
t
i
t
i
+
1
f
(
y
(
t
)
)
d
t
Ces équations montrent comment les valeurs de (y) et (z) peuvent être mises à jour de manière discrétisée en utilisant les intégrales des dérivées respectives sur chaque sous-intervalle ([t_i, t_{i+1}]).

II.B.1) Schéma d’Euler explicite
Pour appliquer le schéma d'Euler explicite aux équations (II.3), nous remplaçons chaque terme sous le signe intégrale par sa valeur prise à la borne inférieure de l'intervalle.

Les équations de récurrence pour les suites 
(
y
i
)
i
∈
J
 et 
(
z
i
)
i
∈
J
 deviennent :

y
i
+
1
=
y
i
+
h
z
i
z
i
+
1
=
z
i
+
h
f
(
y
i
)
où 
y
i
 et 
z
i
 sont des valeurs approchées de 
y
(
t
i
)
 et 
z
(
t
i
)
, et 
h
 est le pas de discrétisation.

II.B.2) Fonction euler
La fonction euler doit recevoir les paramètres suivants :

La fonction 
f
 définie par l'équation différentielle.
Les conditions initiales 
y
0
 et 
z
0
.
Les bornes de l'intervalle 
t
m
i
n
 et 
t
m
a
x
.
Le nombre de points 
n
 pour la discrétisation.

 II.B.3) Illustration de la méthode
a) Quantité 
E
 indépendante du temps
Pour l'équation différentielle (\forall t \in I, y''(t) = -\omega^2 y(t)), nous pouvons définir une quantité 
E
 qui est indépendante du temps :

E
=
1
2
y
′
(
t
)
2
+
1
2
ω
2
y
(
t
)
2
Cette quantité est conservée car :

d
E
d
t
=
y
′
(
t
)
y
″
(
t
)
+
ω
2
y
(
t
)
y
′
(
t
)
=
y
′
(
t
)
(
−
ω
2
y
(
t
)
)
+
ω
2
y
(
t
)
y
′
(
t
)
=
0
b) Variation de 
E
i
On note 
E
i
 la valeur approchée de 
E
 à l'instant 
t
i
 :

E
i
=
1
2
z
i
2
+
1
2
ω
2
y
i
2
Montrons que 
E
i
+
1
−
E
i
=
h
2
ω
2
E
i
 :

E
i
+
1
=
1
2
z
i
+
1
2
+
1
2
ω
2
y
i
+
1
2
En utilisant les relations de récurrence :

y
i
+
1
=
y
i
+
h
z
i
z
i
+
1
=
z
i
+
h
(
−
ω
2
y
i
)
Nous avons :

E
i
+
1
=
1
2
(
z
i
+
h
(
−
ω
2
y
i
)
)
2
+
1
2
ω
2
(
y
i
+
h
z
i
)
2
c) Schéma conservant 
E
Un schéma conservant 
E
 aurait la propriété que 
E
i
+
1
=
E
i
 pour tous les 
i
. Cela signifie que 
d
E
d
t
=
0
.

d) Allure du graphe qui respecte la conservation de 
E
Si 
E
 est conservé, le graphe des valeurs de 
y
i
 et 
z
i
 tracé respectivement sur l'axe des abscisses et des ordonnées serait une ellipse.

e) Justification de l'allure du graphe obtenu par la méthode d'Euler explicite
La méthode d'Euler explicite génère un graphe en spirale divergente, ce qui indique que 
E
 n'est pas conservé. Cette divergence est due à l'erreur de troncature de l'approximation d'Euler explicite, qui accumule des erreurs à chaque pas de temps.

Conclusion
Le schéma d'Euler explicite offre une méthode simple pour résoudre numériquement des équations différentielles, mais il ne conserve pas les quantités indépendantes du temps comme l'énergie, ce qui peut entraîner des erreurs significatives dans les simulations à long terme.

II.C.1) Fonction verlet
La fonction verlet appliquera le schéma de Verlet pour résoudre l'équation différentielle (\forall t \in I, y''(t) = f(y(t))). Elle recevra des paramètres similaires à la fonction euler, mais les relations de récurrence seront différentes.

II.C.2) Comparaison des résultats pour l'oscillateur harmonique
a) Différence d'énergie dans le schéma de Verlet
Pour l'équation différentielle (\forall t \in I, y''(t) = -\omega^2 y(t)), nous avons défini une quantité (E) indépendante du temps :

E
=
1
2
y
′
(
t
)
2
+
1
2
ω
2
y
(
t
)
2
En utilisant le schéma de Verlet, nous pouvons montrer que la variation de (E) est d'ordre (\mathcal{O}(h^3)). Les relations de récurrence du schéma de Verlet sont :

y
i
+
1
=
y
i
+
h
z
i
+
h
2
2
f
i
z
i
+
1
=
z
i
+
h
2
(
f
i
+
f
i
+
1
)
En substituant ces relations dans l'expression de l'énergie, on peut montrer que l'erreur est d'ordre (\mathcal{O}(h^3)).

b) Interprétation du graphe
La mise en œuvre du schéma de Verlet avec les mêmes paramètres que ceux utilisés au II.B.3.e donne un résultat graphique montrant une trajectoire fermée (un cercle ou une ellipse), ce qui indique une meilleure conservation de l'énergie :

Ce graphe montre que l'énergie est conservée beaucoup mieux qu'avec le schéma d'Euler explicite, où la trajectoire divergeait.

c) Conclusion sur le schéma de Verlet
Le schéma de Verlet présente une meilleure conservation de l'énergie par rapport au schéma d'Euler explicite, ce qui en fait un choix supérieur pour les systèmes conservatifs comme les oscillateurs harmoniques. Il permet de maintenir la stabilité des solutions sur des périodes plus longues et fournit une meilleure précision globalement.

III.A.1) Force exercée sur un corps par l'ensemble des autres corps
La force (\vec{F}_j) exercée sur le corps (P_j) par l'ensemble des autres corps (P_k) (avec (k \ne j)) est la somme des forces gravitationnelles exercées par chaque corps (P_k) sur (P_j). La force gravitationnelle exercée par le corps (P_k) sur le corps (P_j) est donnée par :

F
→
k
/
j
=
G
m
j
m
k
r
j
k
3
r
→
j
k
où :

(G) est la constante de gravitation universelle ((6.67 \times 10^{-11} , \text{N} \cdot \text{m}^2 \cdot \text{kg}^{-2})),
(m_j) et (m_k) sont les masses des corps (P_j) et (P_k),
(r_{jk}) est la distance entre les corps (P_j) et (P_k),
(\vec{r}_{jk} = \vec{P_k} - \vec{P_j}) est le vecteur position de (P_k) par rapport à (P_j).
La force totale exercée sur (P_j) par l'ensemble des autres corps est donc :

F
→
j
=
∑
k
≠
j
F
→
k
/
j
III.A.2) Fonction force2
La fonction force2 calcule la force gravitationnelle exercée par le corps 2 sur le corps 1. Elle prend en paramètres les masses (m1 et m2) et les positions (p1 et p2) des deux corps, et renvoie la force sous forme d'une liste de trois composantes cartésiennes.

III.A.3) Fonction forceN
La fonction forceN calcule la force totale exercée sur un corps (P_j) par tous les autres corps du système. Elle prend en paramètres l'indice (j) du corps, la liste des masses des (N) corps, et la liste de leurs positions, et renvoie la force totale sous forme d'une liste de trois composantes cartésiennes.

Résumé
Force exercée sur un corps par l'ensemble des autres corps :

La force totale est la somme des forces gravitationnelles exercées par chaque corps.
Fonction force2 :

Calcule la force gravitationnelle exercée par un corps sur un autre.
Fonction forceN :

Calcule la force totale exercée sur un corps par tous les autres corps du système.

III.B.1) Structure et signification de position[i] et vitesse[i]
position[i] : Cette liste contient les positions des 
N
 corps du système à l'instant 
t
i
. Chaque élément de la liste position[i][j] représente la position du corps 
P
j
 à l'instant 
t
i
 sous forme de liste de trois coordonnées cartésiennes ([x_{ij}, y_{ij}, z_{ij}]) en mètres.

 vitesse[i] : Cette liste contient les vitesses des 
N
 corps du système à l'instant 
t
i
. Chaque élément de la liste vitesse[i][j] représente la vitesse du corps 
P
j
 à l'instant 
t
i
 sous forme de liste de trois composantes cartésiennes ([v_{xij}, v_{yij}, v_{zij}]) en mètres par seconde.

 Explications
Structure et signification de position[i] et vitesse[i] :

position[i] représente les positions des 
N
 corps à l'instant 
t
i
.
vitesse[i] représente les vitesses des 
N
 corps à l'instant 
t
i
.
Fonction pos_suiv :

Calcule les nouvelles positions des corps à l'instant 
t
i
+
1
 en utilisant le schéma de Verlet.
Fonction etat_suiv :

Calcule les nouvelles positions et vitesses des corps à l'instant 
t
i
+
1
 en utilisant le schéma de Verlet.

 III.B.4) Analyse graphique de la complexité
a) Relation entre (\ln(\tau_i)) et (\ln(N))
Si la figure 2 montre une relation linéaire entre (\ln(\tau_i)) (en ordonnée) et (\ln(N)) (en abscisse), cela suggère que (\ln(\tau_i)) est proportionnel à (\ln(N)). En d'autres termes, on peut écrire :

ln
⁡
(
τ
i
)
=
a
ln
⁡
(
N
)
+
b
où (a) et (b) sont des constantes.

b) Hypothèse quant à la complexité de l'algorithme
Si (\ln(\tau_i)) est linéairement proportionnel à (\ln(N)), cela suggère que (\tau_i) est proportionnel à (N^a). En d'autres termes, la complexité temporelle de l'algorithme est probablement de la forme :

τ
i
=
O
(
N
a
)
où (a) est une constante. Cette relation exponentielle est typique de la complexité polynomiale.

III.B.5) Estimation de la complexité temporelle de etat_suiv
a) Estimation de la complexité temporelle
Pour estimer la complexité temporelle de la fonction etat_suiv, examinons les opérations effectuées :

Calcul des forces : La fonction forceN calcule la force exercée sur chaque corps par tous les autres corps. Pour chaque corps (j), elle itère sur tous les autres corps (k). Cela prend un temps proportionnel à (N-1) pour chaque corps, donc la complexité totale pour calculer les forces est (\mathcal{O}(N^2)).

Mise à jour des positions : La fonction pos_suiv met à jour les positions de tous les corps, ce qui prend un temps proportionnel à (N).

Mise à jour des vitesses : La fonction etat_suiv met à jour les vitesses de tous les corps, ce qui prend également un temps proportionnel à (N).

La complexité totale de la fonction etat_suiv est donc dominée par le calcul des forces, ce qui donne une complexité temporelle de (\mathcal{O}(N^2)).

b) Comparaison avec le résultat de la question III.B.4
La relation linéaire observée dans la figure 2 entre (\ln(\tau_i)) et (\ln(N)) suggère une complexité de la forme (\mathcal{O}(N^a)), avec (a) étant la pente de la droite. Si la pente est 2, cela confirme que la complexité temporelle de l'algorithme est (\mathcal{O}(N^2)), ce qui est cohérent avec notre estimation basée sur l'analyse des opérations effectuées par la fonction etat_suiv.

Conclusion
La complexité temporelle de la fonction etat_suiv est (\mathcal{O}(N^2)), ce qui est confirmé à la fois par l'analyse théorique et par l'observation expérimentale des temps de calcul en fonction du nombre de corps (N).

IV.A – Requête SQL pour renvoyer la liste des masses de tous les corps étudiés
Pour obtenir la liste des masses de tous les corps étudiés dans la table CORPS, vous pouvez utiliser la requête SQL suivante :

SELECT masse FROM CORPS;
IV.B.1 – Requête SQL pour vérifier que tous les corps ont un état connu antérieur à tmin()
Pour obtenir le nombre de corps qui ont au moins un état connu antérieur à tmin(), vous pouvez utiliser la requête SQL suivante. Cette requête vérifie que chaque corps a au moins une entrée dans la table ETAT avec une date (datem) antérieure à tmin() :

SELECT COUNT(DISTINCT id_corps) 
FROM ETAT 
WHERE datem < tmin();
IV.B.2 – Requête SQL pour renvoyer l'identifiant et la date du dernier état antérieur à tmin()
Pour obtenir, pour chaque corps, son identifiant et la date de son dernier état antérieur à tmin(), vous pouvez utiliser la requête SQL suivante. Cette requête utilise une sous-requête pour trouver le dernier état de chaque corps avant tmin() :

SQL
SELECT id_corps, MAX(datem) AS dernier_etat
FROM ETAT
WHERE datem < tmin()
GROUP BY id_corps;

Explications
Requête pour obtenir les masses des corps (IV.A) :

La requête SELECT masse FROM CORPS; sélectionne simplement la colonne masse de la table CORPS, renvoyant ainsi la liste des masses de tous les corps étudiés.
Requête pour vérifier les états connus antérieurs à tmin() (IV.B.1) :

La requête SELECT COUNT(DISTINCT id_corps) FROM ETAT WHERE datem < tmin(); compte le nombre de corps distincts (identifiés par id_corps) qui ont au moins une entrée dans la table ETAT avec une date (datem) antérieure à tmin(). La fonction DISTINCT est utilisée pour s'assurer que chaque corps est compté une seule fois.
Requête pour obtenir l'identifiant et la date du dernier état antérieur à tmin() (IV.B.2) :

La requête SELECT id_corps, MAX(datem) AS dernier_etat FROM ETAT WHERE datem < tmin() GROUP BY id_corps; sélectionne l'identifiant de chaque corps (id_corps) et la date maximale (MAX(datem)) parmi les dates antérieures à tmin(). La clause GROUP BY id_corps permet de regrouper les résultats par corps, et la fonction MAX permet de trouver la date la plus récente pour chaque groupe.
Ces requêtes permettent d'extraire les informations nécessaires pour simuler les positions futures des corps célestes en utilisant les méthodes développées dans les parties précédentes.

Explications
Requête SQL (IV.B.3) :

Cette requête sélectionne la masse, la position et la vitesse de chaque corps ayant une masse supérieure à masse_min() et situé dans un cube centré sur l'origine. Les corps sont triés par leur distance à l'origine.
Fonction simulation_verlet (IV.C) :

La fonction simulation_verlet utilise le schéma de Verlet pour simuler les positions des corps sur n itérations avec un pas de temps deltat.
La fonction verlet_step met à jour les positions et les vitesses des corps à chaque étape en utilisant les forces calculées par forceN.
Les résultats sont fournis en unités astronomiques pour les positions et en km/s pour les vitesses.
