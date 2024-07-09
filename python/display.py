import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button, TextBox
from scipy.signal import find_peaks
from parseCSV import *

IGNORE_SIGNALS = ["amplitude(t)", "phase(t)"]

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
    for i in range(len(titles)):
        signals[titles[i]] = np.array(signals[titles[i]])
    return signals

# Fonction pour afficher les signaux
def plot_signals(ax: plt.Axes, signals: dict, title = "Signals", xlabel = "Time", ylabel = "Amplitude"):
    if len(signals) == 0 or "time" not in signals.keys():
        return

    plot_data(ax, signals, 'Time', 'linear', 'time', title, unit='m', ignore_signals=IGNORE_SIGNALS)

# Fonction pour afficher les FFT avec des axes logarithmiques
def plot_FFT(ax: plt.Axes, fft_signals: dict, title = "Tranformée de Fourier", ylabel = "Amplitude (V)", xlabel = "Frequency (Hz)"):
    if len(fft_signals) == 0 or "freq" not in fft_signals.keys():
        return

    plot_data(ax, fft_signals, 'Frequency', 'loglog', 'freq', title)

def calculate_rise_time(time, signal):
    # Définir les seuils bas et haut
    min_signal = np.min(signal)
    max_signal = np.max(signal)
    low_threshold = min_signal + 0.1 * (max_signal - min_signal)
    high_threshold = min_signal + 0.9 * (max_signal - min_signal)

    # Trouver les index où le signal dépasse les seuils
    try:
        low_index = np.where(signal > low_threshold)[0][0]
        high_index = np.where(signal > high_threshold)[0][0]

        # Assurer que high_index est après low_index
        if high_index <= low_index:
            raise ValueError("L'index de seuil haut est avant ou égal à l'index de seuil bas.")

        # prendre à 5 tau
        high_index += 5 * (high_index - low_index)

        # Calculer le temps de montée
        rise_time = time[high_index] - time[low_index]
    except IndexError:
        raise ValueError("Le signal ne dépasse pas les seuils spécifiés.")
    
    return rise_time, low_index, high_index
    

def plot_signal_with_rise_time(ax: plt.Axes, signals: dict, signal_name: str):
    time = signals['time']
    if not signal_name in signals:
        raise ValueError(f"Le signal '{signal_name}' n'existe pas dans les données fournies.")
    signal = signals[signal_name]

    rise_time, low_index, high_index = calculate_rise_time(time, signal)
    print("rise time (5 tau) :", rise_time*1000, "ms, low_index :", low_index, ", high_index :", high_index)

    ax.plot(time * 1000, np.real(signal), label=signal_name)
    #ax.axvline(time[low_index] * 1000, color='r', linestyle='--', label='10% threshold')
    ax.axvline(time[high_index] * 1000, color='g', linestyle='--', label='5 tau')
    ax.set_xlabel('Time (ms)', fontsize=14)
    ax.set_ylabel('Amplitude', fontsize=14)
    ax.set_title(f'Signal with Rise Time: {(rise_time*1000):.2e} ms', fontsize=16)
    ax.legend()
    ax.grid(True)

    # Ajuster la taille de la police des numérotations des axes
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=14)

# Fonction pour basculer entre l'affichage des signaux et des FFT
def toggle_plot(event, next=True):
    global current_plot_index
    if next:
        current_plot_index = (current_plot_index + 1) % 4
    ax.clear()
    ax.set_visible(True)
    if current_plot_index == 0:
        toggle_button.label.set_text("Voir temps de montée")
        plot_signals(ax, signals, title = "Signaux")
    elif current_plot_index == 1:
        toggle_button.label.set_text("Voir FFT signal")
        plot_signal_with_rise_time(ax, signals, "amplitude(t)")
    elif current_plot_index == 2:
        toggle_button.label.set_text("Voir FFT amplitude")
        plot_FFT(ax, FFT_signal, title = "Tranformée de Fourier du signal")
    elif current_plot_index == 3:
        toggle_button.label.set_text("Voir signaux")
        plot_FFT(ax, FFT_amplitude, title = "Tranformée de Fourier de l'amplitude")
    plt.draw()  # Rafraîchir l'écran

# Nom du fichier CSV contenant les signaux et les FFT
filename1 = './data/test_temporel.csv'
filename2 = './data/test_FFT_signal.csv'
filename3 = './data/test_FFT_amplitude.csv'

# Lecture des fichiers CSV
signals         = read_csv(filename1)
FFT_signal      = read_csv(filename2)
FFT_amplitude   = read_csv(filename3)

# Création de la figure et des sous-graphiques
fig, ax = plt.subplots()
gs = fig.add_gridspec(2, hspace=0.3)

# Ajout du bouton pour basculer entre les affichages des signaux et des FFT
toggle_button_ax = plt.axes([0.895, 0.96, 0.1, 0.035])
toggle_button = Button(toggle_button_ax, 'Voir FFT')
toggle_button.on_clicked(toggle_plot)

# Affichage initial des signaux
current_plot_index = -1
toggle_plot(None)

plt.show()
