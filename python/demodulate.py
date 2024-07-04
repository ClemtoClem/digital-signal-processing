import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Paramètres généraux
decimation = 16
num_samples = 16384  # Nombre d'échantillons
Fs = 125e6 / decimation  # Fréquence d'échantillonnage en Hz après décimation
T = num_samples / Fs     # Durée du signal en secondes (nombre d'échantillons / Fréquence d'échantillonnage)
f_coupure = 3e3          # Fréquence de coupure des filtres en Hz

# Générer les différentes fréquences de signal source à tester (de 100 Hz à 1 MHz)
f_signal_values = np.logspace(2, 6, num=20)  # De 100 Hz à 1 MHz, échelle logarithmique

# Initialiser un tableau pour stocker les FFT des amplitudes
fft_amplitudes = np.zeros((num_samples, len(f_signal_values)))

# Boucle sur les différentes fréquences de signal source
for idx, f_signal in enumerate(f_signal_values):
    # Générer le signal sinusoïdal
    t = np.arange(0, T, 1/Fs)
    signal_source = np.sin(2 * np.pi * f_signal * t)

    # Générer la porteuse à une fréquence fixe
    f_porteuse = 3 * f_signal  # Choisir une fréquence de porteuse par exemple
    porteuse = np.sin(2 * np.pi * f_porteuse * t)

    # Démodulation par multiplication avec la porteuse
    signal_demod_I = signal_source * porteuse
    signal_demod_Q = signal_source * np.cos(2 * np.pi * f_porteuse * t)

    # Filtrage avec un filtre passe-bas Butterworth d'ordre 4
    b, a = signal.butter(4, f_coupure / (Fs / 2), btype='low')
    x_filtered = signal.filtfilt(b, a, signal_demod_I)
    y_filtered = signal.filtfilt(b, a, signal_demod_Q)

    # Calcul de l'amplitude
    amplitude = np.sqrt(2 * x_filtered**2 + 2 * y_filtered**2)

    # Calcul de la FFT de l'amplitude et stockage dans fft_amplitudes
    fft_amplitudes[:, idx] = np.abs(np.fft.fft(amplitude))

# Création du plot en log-log
plt.figure(figsize=(12, 6))

# Définir l'axe des fréquences en Hz
freq_axis = np.linspace(0, Fs/8, num=num_samples)

# Boucle sur les différentes fréquences de signal source pour tracer les FFT
for idx in range(len(f_signal_values)):
    plt.loglog(freq_axis, fft_amplitudes[:, idx])

plt.title('FFT de l\'amplitude du signal démodulé')
plt.xlabel('Fréquence du signal source (Hz)')
plt.ylabel('Amplitude de la FFT (log scale)')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()
