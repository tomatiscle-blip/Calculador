from pynite import FEModel3D

def analizar_tablero_losa(lx, ly, q_carga, espesor, E_c, nu_c):
    """
    Crea un modelo espacial en PyNite para un tablero rectangular de losa,
    aplica una carga uniforme y devuelve las reacciones perimetrales.
    """
    # 1. Crear el modelo
    modelo = FEModel3D()
    
    # 2. Definir materiales y secciones (o espesor de placa)
    # PyNite permite agregar placas con 'add_quad' o mallados automáticos.
    # (Para simplificar, defines los nodos de la malla y los elementos Quad).
    
    # 3. Aplicar condiciones de borde en el perímetro
    # Restringir los desplazamientos verticales (y rotaciones según corresponda)
    # en los nodos que coinciden con la viga de borde.
    
    # 4. Aplicar cargas
    # Carga uniformemente repartida 'q' sobre los elementos Quad.
    
    # 5. Resolver el modelo
    modelo.analyze()
    
    # 6. Extraer reacciones nodales de apoyo
    # Esto te dará: Fz (corte vertical), Mx y My (momentos torsores y flectores en borde).
    reacciones = {}
    # Recorrer nodos de borde y extraer modelo.nodes[nombre].Rx, Ry, Rz, etc.
    
    return reacciones