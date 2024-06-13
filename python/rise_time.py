import numpy as np
import matplotlib.pyplot as plt

def calculate_rise_time(time, signal):
    min_val = np.min(signal)
    max_val = np.max(signal)
    low_threshold = min_val + 0.1 * (max_val - min_val)
    high_threshold = min_val + 0.9 * (max_val - min_val)

    low_index = np.where(signal >= low_threshold)[0][0]
    high_index = np.where(signal >= high_threshold)[0][0]

    rise_time = time[high_index] - time[low_index]
    
    return rise_time, low_index, high_index

def plot_signal_with_rise_time(time, signal):
    rise_time, low_index, high_index = calculate_rise_time(time, signal)

    plt.plot(time, signal, label='Signal')
    plt.axvline(time[low_index], color='r', linestyle='--', label='10% threshold')
    plt.axvline(time[high_index], color='g', linestyle='--', label='90% threshold')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title(f'Signal with Rise Time: {rise_time:.4f} s')
    plt.legend()
    plt.grid(True)
    plt.show()

    print(f'Rise Time: {rise_time:.4f} seconds')

# Exemple d'utilisation avec des données fictives
time = np.linspace(0, 1, 1000)
signal = np.tanh(5 * (time - 0.5))  # Signal fictif pour démonstration

plot_signal_with_rise_time(time, signal)
