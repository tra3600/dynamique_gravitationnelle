SELECT 
    c.masse, 
    e.x, 
    e.y, 
    e.z, 
    e.vx, 
    e.vy, 
    e.vz
FROM 
    CORPS c
JOIN 
    ETAT e ON c.id_corps = e.id_corps
JOIN 
    date_mesure d ON c.id_corps = d.id_corps AND e.datem = d.date_der
WHERE 
    c.masse >= masse_min()
    AND ABS(e.x) <= arete() / 2 
    AND ABS(e.y) <= arete() / 2
    AND ABS(e.z) <= arete() / 2
ORDER BY 
    SQRT(e.x * e.x + e.y * e.y + e.z * e.z);