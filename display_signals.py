import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button

# Fonction pour lire les signaux à partir du fichier CSV
def read_csv(filename):
    signals = {}
    titles = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        for title in lines[0].strip('\n').split(","):
            titles.append(title)
            signals[title] = []

        for line in lines[1:]:  # Commencer à partir de la deuxième ligne pour les données
            values = line.strip('\n').split(",")
            for i, val in enumerate(values):
                try:
                    if titles[i] == "time" or val[-1] != "j":
                        signals[titles[i]].append(float(val))
                    else:
                        signals[titles[i]].append(complex(val))
                except:
                    signals[titles[i]].append(complex(0))
    return signals

# Fonction pour afficher les signaux
def plot_signals(ax, signals: dict):
    if len(signals) == 0 or "time" not in signals.keys():
        return

    for signal_name, signal_data in signals.items():
        print(signal_name)
        if signal_name != "time":
            ax.plot(signals["time"], np.real(signal_data), label=signal_name)  # Afficher l'amplitude absolue
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude')
    ax.set_title('Signals')
    ax.legend()
    ax.grid(True)

# Fonction pour afficher les FFT avec des axes logarithmiques
def plot_fft_signals(ax_mag, ax_phase, fft_signals):
    if len(fft_signals) == 0 or "time" not in fft_signals.keys():
        return

    sample_frequency = 1/(fft_signals["time"][1] - fft_signals["time"][0])
    for signal_name, signal_data in fft_signals.items():
        if signal_name is ["signal", "tempA_fter", "tempA"]:
            freqs = np.fft.fftfreq(len(signal_data), d=1/sample_frequency)
            magnitude_spectrum = np.abs(np.imag(signal_data))
            phase_spectrum = np.angle(np.imag(signal_data))
            ax_mag.semilogx(freqs, magnitude_spectrum, label=signal_name)
            ax_phase.semilogx(freqs, phase_spectrum, label=signal_name)

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
        plot_signals(ax, signals)
        ax.set_visible(True)
        ax_mag.set_visible(False)
        ax_phase.set_visible(False)
        toggle_button.label.set_text("Voir FFT")
    else:
        plot_fft_signals(ax_mag, ax_phase, fft_signals)
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
signals = read_csv(filename)

# Lecture des FFT à partir du fichier CSV
fft_signals = read_csv(fft_filename)

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
plot_signals(ax, signals)
ax_mag.set_visible(False)
ax_phase.set_visible(False)

plt.show()
