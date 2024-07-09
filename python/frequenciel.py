import numpy as np
import matplotlib.pyplot as plt
from parseCSV import *

if __name__ == "__main__":
    # Spécifiez le chemin vers votre fichier CSV
    fichier_csv = "./data/test_FFT_signal.csv"
    
    # Lecture des données à partir du fichier CSV
    signals: dict = read_csv(fichier_csv)

    # Appel de la fonction pour tracer les données
    ax = plt.subplot(111)
    plot_data(ax, signals, "Frequency", axis_form = "loglog", title = "Transformée de fourirer discrète du signal")
    plt.show()