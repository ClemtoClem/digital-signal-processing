import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from scipy import signal

N = 1024

# Fréquence d'échantillonnage
fe = 126000  # Hz

# Période d'échantillonnage
Te = 1/fe  # s

# Fréquence de nyquist
f_nyq = fe / 2.0  # Hz

n_vect = np.arange(N+1) # 0 1 2 3 ... N

# Génération d'un signal à filtrer
x = np.linspace(0, Te*N, 1000)
f1 = 1000
f2 = 10000
sig = np.sin(2 * np.pi * f1 * x) + 0.2 * np.sin(2 * np.pi * f2 * x)

# Fréquence de coupure
fc = 5000  # Hz
tau = 1/(2*np.pi*fc)

# Fenêtre Idéale
# !! dans numpy, sinc(x)=sin(pi.x)/(pi.x)
h_d = 2*(fc/fe)*np.sinc((2*(fc/fe))*(n_vect-int((N+1)/2)))

# ======== Passe bas ========
# Préparation de la liste de sortie
s_pb = []
s_pb.append(0)

# Application du filtre
for i in range(1, len(sig)):
    s_pb.append(s_pb[i-1]+Te/tau*(sig[i-1]-s_pb[i-1]))

# Préparation du filtre de Butterworth en passe-low
b, a = signal.iirfilter(4, Wn=fc, fs=fe, btype="low", ftype="butter")

# Application du filtre
s_but_pb = signal.filtfilt(b, a, sig)

# ======== Passe haut ========
# Préparation de la liste de sortie
s_ph = []
s_ph.append(0)

# Application du filtre
for i in range(1, len(sig)):
    s_ph.append(s_ph[i-1]*(1-Te/tau)+sig[i]-sig[i-1])

# Préparation du filtre de Butterworth en passe-haut
b, a = signal.iirfilter(4, Wn=fc, fs=fe, btype="high", ftype="butter")

# Application du filtre
s_but_ph = signal.filtfilt(b, a, sig)

def plot_signals(axis):
    # Affichage du signal filtré
    axis[0].plot(x, sig, color='silver', label='Signal')
    axis[0].plot(x, s_pb, label='1er ordre')
    axis[0].plot(x, s_but_pb, color='#cc0000', label='Butterworth')
    axis[0].set_title("Filtre passe-bas numérique")

    axis[1].plot(x, sig, color='silver', label='Signal')
    axis[1].plot(x, s_ph, label='1er ordre')
    axis[1].plot(x, s_but_ph, color='#cc0000', label='Butterworth')
    axis[1].set_title("Filtre passe-haut numérique")

def plot_frequency(axis):
    # Affichage du signal filtré en fréquence

    w = np.fft.fftfreq(N, d=Te)

    H_pb = np.fft.fft(s_pb, N)
    H_but_b = np.fft.fft(s_but_pb, N)
    
    axis[0].plot(w[:N//2], 20*np.log10(np.abs(H_pb)[:N//2]), label='1er ordre')
    axis[0].plot(w[:N//2], 20*np.log10(np.abs(H_but_b)[:N//2]), color='#cc0000', label='Butterworth')
    axis[0].grid(True, which='both')
    axis[0].set_xlabel("Frequence [Hz]")
    axis[0].set_ylabel("Module $H(e^{j\omega})\ en\ dB$")
    axis[0].set_title("Filtre passe-bas numérique")

    axis[1].plot(w[:N//2], np.unwrap(np.angle(H_pb)[:N//2]), label='1er ordre')
    axis[1].plot(w[:N//2], np.unwrap(np.angle(H_but_b)[:N//2]), color='#cc0000', label='Butterworth')
    axis[1].grid(True, which='both')
    axis[1].set_xlabel("Frequence [Hz]")
    axis[1].set_ylabel("Phase $H(e^{j\omega})\ en\ rad$")
    axis[1].set_title("Filtre passe-bas numérique")

    H_ph = np.fft.fft(s_ph, N)
    H_but_h = np.fft.fft(s_but_ph, N)
    
    axis[2].plot(w[:N//2], 20*np.log10(np.abs(H_ph)[:N//2]), label='1er ordre')
    axis[2].plot(w[:N//2], 20*np.log10(np.abs(H_but_h)[:N//2]), color='#cc0000', label='Butterworth')
    axis[2].grid(True, which='both')
    axis[2].set_xlabel("Frequence [Hz]")
    axis[2].set_ylabel("Module $H(e^{j\omega})\ en\ dB$")
    axis[2].set_title("Filtre passe-haut numérique")

    axis[3].plot(w[:N//2], np.unwrap(np.angle(H_ph)[:N//2]), label='1er ordre')
    axis[3].plot(w[:N//2], np.unwrap(np.angle(H_but_h)[:N//2]), color='#cc0000', label='Butterworth')
    axis[3].grid(True, which='both')
    axis[3].set_xlabel("Frequence [Hz]")
    axis[3].set_ylabel("Phase $H(e^{j\omega})\ en\ rad$")
    axis[3].set_title("Filtre passe-haut numérique")

# Fonction pour basculer entre l'affichage des signaux et des FFT
def toggle_plot(event):
    global current_plot_index
    current_plot_index = 1-current_plot_index

    if current_plot_index == 0:
        for ax in axis_B:
            ax.clear()
        plot_signals(axis_A)
        for ax in axis_B:
            ax.set_visible(False)
        for ax in axis_A:
            ax.set_visible(True)
            ax.legend()
            ax.grid(True, which='both')
        toggle_button.label.set_text("Voir FFT")
    else:
        for ax in axis_A:
            ax.clear()
        plot_frequency(axis_B)
        for ax in axis_A:
            ax.set_visible(False)
        for ax in axis_B:
            ax.set_visible(True)
            ax.legend()
            ax.grid(True, which='both')
        toggle_button.label.set_text("Voir signaux")
    
    plt.draw()  # Rafraîchir l'écran

# Création de la figure et des sous-graphiques
fig, ax = plt.subplots(2, 2)
gs = fig.add_gridspec(2, 2, hspace=0.3)
axis_A = ax.flatten()[:2]
axis_B = ax.flatten()

# Ajout du bouton pour basculer entre les affichages des signaux et des FFT
toggle_button_ax = plt.axes([0.88, 0.01, 0.1, 0.075])
toggle_button = Button(toggle_button_ax, 'Voir FFT')
toggle_button.on_clicked(toggle_plot)

# Affichage initial des signaux
current_plot_index = 1
toggle_plot(None)

plt.show()
