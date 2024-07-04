import pandas as pd
import matplotlib.pyplot as plt
from parseCSV import *

if __name__ == "__main__":
    # Spécifiez le chemin vers votre fichier CSV
    fichier_csv = "./data/demodulation.csv"
    
    # Lecture des données à partir du fichier CSV
    signals: dict = read_csv(fichier_csv)

    # Appel de la fonction pour tracer les données
    plot_data(signals, "Time", axis_form = "linear", title = "Affichage des signaux", unit="m")
