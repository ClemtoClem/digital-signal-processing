import numpy as np
import matplotlib.pyplot as plt
from parseCSV import *

if __name__ == "__main__":
    # Spécifiez le chemin vers votre fichier CSV
    fichier_csv = "./data/test_temporel.csv"
    
    # Lecture des données à partir du fichier CSV
    signals: dict = read_csv(fichier_csv)

    time = signals["time"]
    signal = signals["output(t)"]

    fft = np.fft.fft(signal, n=len(signal))
    fft_freq = np.fft.fftfreq(len(signal), d=time[1]-time[0])

    print(len(fft_freq))
    print(max(fft_freq))

    plt.loglog(fft_freq, np.abs(fft) * (np.sqrt(2.0)/len(signal)))
    plt.xlabel("Frequence [Hz]")
    plt.ylabel("Amplitude [V]")
    plt.title("Transformée de Fourier discrète du signal")
    plt.show()