import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button

# Fonction pour lire les signaux à partir du fichier CSV
def read_csv(filename):
    signals = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        if len(lines) < 2:  # Vérifier si la ligne contient suffisamment d'éléments
            return 0, []
        sample_frequency = int(lines[0])
        for line in lines[1:]:
            row = line.strip().split(" ")
            if len(row) < 2 or len(row) % 2 != 0:  # Vérifier si la ligne contient suffisamment d'éléments et si elle a un nombre pair d'éléments
                continue
            name = row[0]
            size = int(row[1])
            values = [complex(float(row[i]), float(row[i+1])) for i in range(2, len(row), 2)]
            signals.append([name, size, values])
    return sample_frequency, signals

# Fonction pour afficher les signaux
def plot_signals(ax, sample_frequency, signals):
    if len(signals) == 0:
        return

    time = np.arange(0, signals[0][1]) / sample_frequency
    for signal in signals:
        ax.plot(time, np.real(signal[2]), label=signal[0] + ' (Real)')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude')
    ax.set_title('Signals')
    ax.legend()
    ax.grid(True)

# Fonction pour afficher les FFT avec des axes logarithmiques
def plot_fft_signals(ax_mag, ax_phase, sample_frequency, fft_signals):
    if len(fft_signals) == 0:
        return

    for fft_signal in fft_signals:
        freqs = np.fft.fftfreq(fft_signal[1]) * sample_frequency
        magnitude_spectrum = np.abs(fft_signal[2])
        phase_spectrum = np.angle(fft_signal[2])
        ax_mag.semilogx(freqs, magnitude_spectrum, label=fft_signal[0])
        ax_phase.semilogx(freqs, phase_spectrum, label=fft_signal[0])

    ax_mag.set_xlabel('Frequency (Hz)')
    ax_mag.set_ylabel('Magnitude')
    ax_mag.set_title('FFT Magnitude')
    ax_mag.legend()
    ax_mag.grid(True)

    ax_phase.set_xlabel('Frequency (Hz)')
    ax_phase.set_ylabel('Phase (radians)')
    ax_phase.set_title('FFT Phase')
    ax_phase.legend()
    ax_phase.grid(True)

# Fonction pour basculer entre l'affichage des signaux et des FFT
def toggle_plot(event):
    global current_plot_index
    current_plot_index = (current_plot_index + 1) % 2
    ax.clear()
    ax_mag.clear()
    ax_phase.clear()
    if current_plot_index == 0:
        plot_signals(ax, sample_frequency, signals)
        ax.set_visible(True)
        ax_mag.set_visible(False)
        ax_phase.set_visible(False)
        toggle_button.label.set_text("Voir FFT")
    else:
        plot_fft_signals(ax_mag, ax_phase, sample_frequency, fft_signals)
        ax.set_visible(False)
        ax_mag.set_visible(True)
        ax_phase.set_visible(True)
        toggle_button.label.set_text("Voir signaux")
    plt.draw()  # Rafraîchir l'écran

# Nom du fichier CSV contenant les signaux
filename = 'test.csv'

# Nom du fichier CSV contenant les FFT
fft_filename = 'test_fft.csv'

# Lecture des signaux à partir du fichier CSV
sample_frequency, signals = read_csv(filename)

# Lecture des FFT à partir du fichier CSV
_, fft_signals = read_csv(fft_filename)

# Création de la figure et des sous-graphiques
fig, ax = plt.subplots()
gs = fig.add_gridspec(2, hspace=0.3)
ax_mag = fig.add_subplot(gs[0, 0])
ax_phase = fig.add_subplot(gs[1, 0])

# Ajout du bouton pour basculer entre les affichages des signaux et des FFT
toggle_button_ax = plt.axes([0.88, 0.01, 0.1, 0.075])
toggle_button = Button(toggle_button_ax, 'Voir FFT')
toggle_button.on_clicked(toggle_plot)

# Affichage initial des signaux
current_plot_index = 0
plot_signals(ax, sample_frequency, signals)
ax_mag.set_visible(False)
ax_phase.set_visible(False)

plt.show()
