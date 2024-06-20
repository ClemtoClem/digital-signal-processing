import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button, TextBox
from scipy.signal import find_peaks

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
def plot_signals(ax: plt.Axes, signals: dict):
    if len(signals) == 0 or "time" not in signals.keys():
        return

    for signal_name, signal_data in signals.items():
        print(signal_name)
        if signal_name.lower() != "time":
            ax.plot(signals["time"], np.real(signal_data), label=signal_name)  # Afficher l'amplitude absolue
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude')
    ax.set_title('Signals')
    ax.legend()
    ax.grid(True)

# Fonction pour afficher les FFT avec des axes logarithmiques
def plot_fft_signals(ax: plt.Axes, fft_signals: dict, min_freq=None, max_freq=None):
    if len(fft_signals) == 0 or "freq" not in fft_signals.keys():
        return

    sample_frequency = 1/(fft_signals["freq"][1] - fft_signals["freq"][0])
    for signal_name, signal_data in fft_signals.items():
        if signal_name.lower() != "freq":
            magnitude_spectrum = np.abs(np.imag(signal_data))
            ax.loglog(fft_signals["freq"], magnitude_spectrum, label=signal_name)

            if min_freq is not None and max_freq is not None:
                indices = np.where((fft_signals["freq"] >= min_freq) & (fft_signals["freq"] <= max_freq))
                selected_freqs = fft_signals["freq"][indices]
                selected_magnitudes = magnitude_spectrum[indices]
                
                peaks, _ = find_peaks(selected_magnitudes)
                peak_freqs = selected_freqs[peaks]
                peak_magnitudes = selected_magnitudes[peaks]

                ax.plot(peak_freqs, peak_magnitudes, 'ro', label='Peaks')
                for i, freq in enumerate(peak_freqs):
                    ax.text(freq, peak_magnitudes[i], f'{freq:.1f} Hz', color='red')

                # Ajouter des lignes pointillées aux positions des fréquences min et max
                ax.axvline(min_freq, color='r', linestyle='--', label='Min Freq')
                ax.axvline(max_freq, color='g', linestyle='--', label='Max Freq')

    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('dB')
    ax.set_title('FFT dB')
    ax.legend()
    ax.grid(True)

def calculate_rise_time(time, signal):
    min_val = np.min(signal)
    max_val = np.max(signal)
    low_threshold = min_val + 0.1 * (max_val - min_val)
    high_threshold = min_val + 0.9 * (max_val - min_val)

    low_index = np.where(signal >= low_threshold)[0][0]
    high_index = np.where(signal >= high_threshold)[0][0]

    rise_time = time[high_index] - time[low_index]
    
    return rise_time, low_index, high_index

def plot_signal_with_rise_time(ax: plt.Axes, signals: dict, signal_name: str):
    time = signals['time']
    if not(signal_name in signals.keys()):
        return
    signal = signals[signal_name]

    rise_time, low_index, high_index = calculate_rise_time(time, signal)
    print("rise time :", rise_time, "s, low_index :", time[low_index], "s, high_index :", time[high_index], "s")

    ax.plot(time, np.real(signal), label=signal_name)
    ax.axvline(time[low_index], color='r', linestyle='--', label='10% threshold')
    ax.axvline(time[high_index], color='g', linestyle='--', label='90% threshold')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude')
    ax.set_title(f'Signal with Rise Time: {rise_time} s')
    ax.legend()
    ax.grid(True)

# Fonction pour basculer entre l'affichage des signaux et des FFT
def toggle_plot(event, next=True):
    global current_plot_index
    if next:
        current_plot_index = (current_plot_index + 1) % 3
    ax.clear()
    ax.set_visible(True)
    set_widgets_visibility()
    if current_plot_index == 0:
        toggle_button.label.set_text("Voir temps de montée")
        plot_signals(ax, signals)
    elif current_plot_index == 1:
        toggle_button.label.set_text("Voir FFT")
        plot_signal_with_rise_time(ax, signals, "amplitude(t)")
    else:
        toggle_button.label.set_text("Voir signaux")
        plot_fft_signals(ax, fft_signals)
    plt.draw()  # Rafraîchir l'écran

def set_widgets_visibility():
    visible: bool = (current_plot_index == 2)
    min_freq_text_box_ax.set_visible(visible)
    max_freq_text_box_ax.set_visible(visible)
    apply_button_ax.set_visible(visible)

def update_freq_range(event):
    global min_freq, max_freq
    try:
        min_freq = float(min_freq_text_box.text)
        max_freq = float(max_freq_text_box.text)
        toggle_plot(None, False)  # Mettre à jour l'affichage des FFT avec la nouvelle plage de fréquences
    except ValueError:
        print("Entrée non valide pour les fréquences")

# Nom du fichier CSV contenant les signaux
filename = './data/test_signals.csv'

# Nom du fichier CSV contenant les FFT
fft_filename = './data/test_spectrums.csv'

# Lecture des signaux à partir du fichier CSV
signals = read_csv(filename)

# Lecture des FFT à partir du fichier CSV
fft_signals = read_csv(fft_filename)

# Création de la figure et des sous-graphiques
fig, ax = plt.subplots()
gs = fig.add_gridspec(2, hspace=0.3)

# Ajout du bouton pour basculer entre les affichages des signaux et des FFT
toggle_button_ax = plt.axes([0.88, 0.01, 0.1, 0.075])
toggle_button = Button(toggle_button_ax, 'Voir FFT')
toggle_button.on_clicked(toggle_plot)

# Ajout des zones de saisie pour la plage de fréquences
min_freq_text_box_ax = plt.axes([0.1, 0.01, 0.1, 0.05])
min_freq_text_box = TextBox(min_freq_text_box_ax, 'Min Freq (Hz)', initial="0.0")
max_freq_text_box_ax = plt.axes([0.3, 0.01, 0.1, 0.05])
max_freq_text_box = TextBox(max_freq_text_box_ax, 'Max Freq (Hz)', initial="1000.0")
apply_button_ax = plt.axes([0.5, 0.01, 0.1, 0.05])
apply_button = Button(apply_button_ax, 'Apply')
apply_button.on_clicked(update_freq_range)

# Variables globales pour la plage de fréquences
min_freq = 0.0
max_freq = 1000.0

# Affichage initial des signaux
current_plot_index = -1
toggle_plot(None)

plt.show()
