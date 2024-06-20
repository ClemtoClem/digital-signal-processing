import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Chemin vers le fichier CSV
csv_file = './data/noiseRMS.csv'  # Remplacez par le chemin correct de votre fichier

# Lecture des données depuis le fichier CSV
df = pd.read_csv(csv_file)

# Extraction des colonnes Nb points par période et Noise RMS
nb_points = df['Nb points par période']
noise_rms = df['Noise RMS']

# Création du plot
plt.figure(figsize=(10, 6))
plt.semilogx(nb_points, noise_rms, marker='o', linestyle='-', color='b', label='Noise RMS')

# Ajout des titres et labels
plt.title('Noise RMS en fonction du nombre de points par période')
plt.xlabel('Nb points par période')
plt.ylabel('Noise RMS')
plt.grid(True)
plt.legend()

# Affichage du plot
plt.tight_layout()
plt.show()