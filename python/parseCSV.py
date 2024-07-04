import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_csv(filename, typePlot = ""):
    """
        Fonction pour lire un fichier CSV et le convertir en dictionnaire
        :param filename: Le chemin vers le fichier CSV
        :param typePlot: Le type de graphique à tracer (Time, Frequency, PointsPerPeriod)
        :return: Un dictionnaire avec les titres comme clés et les données comme valeurs
    """
    axis_name = ""
    if typePlot == "Time":
        axis_name = 'time'
    elif typePlot == "Frequency":
        axis_name = 'freq'
    elif typePlot == "PointsPerPeriod":
        axis_name = 'pointsPerPeriod'

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
                    if titles[i] == axis_name or val[-1] != "j":
                        signals[titles[i]].append(float(val))
                    else:
                        signals[titles[i]].append(complex(val))
                except:
                    signals[titles[i]].append(0)
    for key, value in signals.items():
        signals[key] = np.array(value)
    return signals

def plot_data(signals, typePlot, axis_form = "linear", axis_name = "", title = "", ignore_signals=[], unit = ""):
    """
        Fonction pour tracer les données
        :param signals: Le dictionnaire avec les données à tracer
        :param typePlot: Le type de graphique à tracer (Time, Frequency, PointsPerPeriod)
        :param axis_form: Le type d'axe à utiliser pour tracer les données (linear, semilogx, semilogy, loglog)
        :param axis_name: Le nom de l'axe à utiliser pour tracer les données
        :param title: Le titre du graphique
        :param ignore_signals: La liste des signaux à ignorer lors du tracé
        :param unit: L'unité à utiliser pour tracer les données (k, m, u, n)
    """
    
    if axis_name == "":
        if typePlot == "Time":
            axis_name = 'time'
        elif typePlot == "Frequency":
            axis_name = 'freq'
        elif typePlot == "PointsPerPeriod":
            axis_name = 'pointsPerPeriod'
    axis = signals[axis_name]

    unit_multiplicator = 1
    if (unit == "k"):
        unit_multiplicator = 1/1000
    elif (unit == "m"):
        unit_multiplicator = 1000
    elif (unit == "u"):
        unit_multiplicator = 1000000
    elif (unit == "n"):
        unit_multiplicator = 1000000000
    
    # Parcourir les colonnes restantes pour tracer chaque série de données
    for name, signal in signals.items():
        if name != axis_name:
            if len(axis) == len(signal) and (name not in ignore_signals):
                if axis_name == "freq":
                    signal = np.absolute(signal)*2/(len(signal))
                if axis_form == "linear":
                    plt.plot(axis*unit_multiplicator, signal, label=name)
                elif axis_form == "semilogx":
                    plt.semilogx(axis*unit_multiplicator, signal, label=name)
                elif axis_form == "semilogy":
                    plt.semilogy(axis*unit_multiplicator, signal, label=name)
                elif axis_form == "loglog":
                    plt.loglog(axis*unit_multiplicator, signal, label=name)
                else:
                    print("Error: axis_form not recognized")
                    print("The possible values are: regular, semilogx, semilogy, loglog")
                    exit(1)
    
    # Ajouter des étiquettes et une légende
    if typePlot == "Time":
        plt.xlabel('Temps [' + unit + 's]')
        plt.ylabel('Amplitude [V]')
    elif typePlot == "Frequency":
        plt.xlabel('Fréquence [' + unit + 'Hz] (échelle log)')
        plt.ylabel('Amplitude [V]')
    elif typePlot == "PointsPerPeriod":
        plt.xlabel('Nombre de périodes')
        plt.ylabel('Amplitude [V]')
    else:
        print("Error: typePlot not recognized")
        print("The possible values are: Time, Frequency, PointsPerPeriod")
        exit(1)
    if title != "":
        plt.title(title)
    else:
        plt.title(f"Graphique à partir du fichier CSV ({typePlot})")
    plt.legend()
    plt.grid(True)
    
    # Afficher le plot
    plt.show()