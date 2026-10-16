def generar_punto_diagrama_interaccion(c, b, h, fck, fy, armadura_capas):
    """
    c: profundidad del eje neutro
    armadura_capas: lista con [(Area_i, distancia_d_i), ...]
    """
    epsilon_cu = 0.003
    Es = 200000  # MPa (módulo del acero)
    
    # 1. Resultante del hormigón (Bloque rectangular de Whitney)
    beta1 = 0.85 if fck <= 30 else max(0.65, 0.85 - 0.05 * (fck - 30) / 7)
    a = beta1 * c
    Cc = 0.85 * fck * a * b  # Fuerza en el hormigón
    
    # 2. Sumatoria de fuerzas en el acero
    F_acero_total = 0
    M_acero_total = 0
    h_baricentro = h / 2
    
    for As_i, d_i in armadura_capas:
        # Deformación en la capa i por compatibilidad
        epsilon_si = epsilon_cu * (c - d_i) / c
        
        # Tensión en el acero (Ley constitutiva elasto-plástica)
        sigma_si = Es * epsilon_si
        if sigma_si > fy:
            sigma_si = fy
        elif sigma_si < -fy:
            sigma_si = -fy
            
        # Fuerza de la capa (positiva compresión, negativa tracción)
        F_i = As_i * sigma_si
        F_acero_total += F_i
        
        # Momento de la capa respecto al baricentro
        brazo = h_baricentro - d_i
        M_acero_total += F_i * brazo
        
    # 3. Equilibrio de la sección (Carga nominal Pn y Momento nominal Mn)
    Pn = Cc + F_acero_total
    
    # Momento del bloque de hormigón respecto al baricentro
    M_cc = Cc * (h_baricentro - a / 2)
    Mn = M_cc + M_acero_total
    
    return Pn, Mn