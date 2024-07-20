import numpy as np
import matplotlib.pyplot as plt
from parseCSV import *

if __name__ == "__main__":
	folder = input("Nom dossier contenant les fichiers CSV : ")
	fichier_csv = "./data/" + folder + "/spectrums.csv"
	
	# Lecture des données à partir du fichier CSV
	fft_signals, _ = read_csv(fichier_csv)

	ignore_signals = []

	# Appel de la fonction pour tracer les données
	fig, ax = plt.subplots(figsize=(16/1.2, 9/1.8))
	plot_data(ax, fft_signals, "frequency", axis_form = "loglog", unit="k", usedColorMap=False, ignore_signals = ignore_signals)
	plt.show()
