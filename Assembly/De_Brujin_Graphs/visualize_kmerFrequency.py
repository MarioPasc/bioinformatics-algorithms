import pandas as pd
import matplotlib.pyplot as plt

# Simulación de entrada: string multilínea como si fuera la salida de la consola
data = """
K-mero: CG, Frecuencia: 1833101
K-mero: GC, Frecuencia: 2271922
K-mero: TG, Frecuencia: 4276311
K-mero: CT, Frecuencia: 3511102
K-mero: TA, Frecuencia: 6362665
K-mero: GG, Frecuencia: 1917435
K-mero: TT, Frecuencia: 8178103
K-mero: AT, Frecuencia: 7462032
K-mero: GT, Frecuencia: 3481031
K-mero: AA, Frecuencia: 8201271
K-mero: TC, Frecuencia: 3814819
K-mero: CA, Frecuencia: 4127012
K-mero: AC, Frecuencia: 3380148
K-mero: AG, Frecuencia: 3199112
K-mero: CC, Frecuencia: 1625630
K-mero: GA, Frecuencia: 3546511
"""

# Procesar los datos
data = data.strip().split("\n")
kmer_data = [line.split(", ") for line in data]
kmer_info = [(item[0].split(": ")[1], int(item[1].split(": ")[1])) for item in kmer_data]

# Convertir a DataFrame
df = pd.DataFrame(kmer_info, columns=["K-mero", "Frecuencia"])

# Ordenar los datos por frecuencia para mejor visualización
df = df.sort_values(by="K-mero", ascending=True)
# Crear un gráfico de barras
plt.figure(figsize=(10, 8))
plt.barh(df["K-mero"], df["Frecuencia"], color='skyblue')
plt.xlabel('Frecuencia')
plt.ylabel('K-mero')
plt.title('Frecuencia de K-mers')
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()
