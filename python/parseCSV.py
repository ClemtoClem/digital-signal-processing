import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def read_csv(filename, axis_name=""):
    """
    Fonction pour lire un fichier CSV et le convertir en dictionnaire
    :param filename: Le chemin vers le fichier CSV
    :param typePlot: Le type de graphique à tracer (time, freq, ...)
    :return: Un dictionnaire avec les titres comme clés et les données comme valeurs
    """

    signals = {}
    titles = []

    with open(filename, 'r') as file:
        lines = file.readlines()

        # Les lignes suivantes contiennent les valeurs des signaux
        for line in lines:
            values = line.split(",")
            title = values[0].strip()  # Supprimer les espaces autour du titre
            signals[title] = []
            titles.append(title)

            # Conversion des données en float ou en complex
            for val in values[1:]:
                val = val.strip()  # Supprimer les espaces autour de la valeur
                try:
                    if title == axis_name or 'j' not in val:
                        signals[title].append(float(val))
                    else:
                        signals[title].append(complex(val))
                except ValueError:
                    signals[title].append(0)

    # Conversion des listes en arrays numpy
    for key, value in signals.items():
        signals[key] = np.array(value)

    return signals, titles

def plot_data(ax: plt.Axes, signals:dict, plot_type, axis_form="linear", title="", ignore_signals=[], unit="", 
			  title_fontsize=16, axis_label_fontsize=14, tick_label_fontsize=14, usedColorMap = False, axis_name=""):
	"""
	Fonction pour tracer les données
	:param ax: L'axe à utiliser pour tracer les données
	:param signals: Le dictionnaire avec les données à tracer
	:param plot_type: Le type de graphique à tracer (Time, Frequency, PointsPerPeriod)
	:param axis_form: Le type d'axe à utiliser pour tracer les données (linear, semilogx, semilogy, loglog)
	:param title: Le titre du graphique
	:param ignore_signals: La liste des signaux à ignorer lors du tracé
	:param unit: L'unité à utiliser pour tracer les données (k, m, u, n)
	:param title_fontsize: La taille de la police du titre du graphique
	:param axis_label_fontsize: La taille de la police des étiquettes des axes
	:param tick_label_fontsize: La taille de la police des numérotations des axes
	:param usedColorMap: La colormap à utiliser pour tracer les données
	:param axis_name: Le nom de l'axe à utiliser pour tracer les données
	"""

	axis_type = plot_type
	if plot_type == "freq_radian":
		axis_type = "frequency"
	elif plot_type == "time_radian":
		axis_type = "time"
	elif plot_type == "loop":
		axis_type = "loop_axis"
	axis = None
	if axis_name != "":
		axis = signals[axis_name]
	else:
		axis = signals[axis_type]

	unit_multiplicator = 1
	if unit == "k":
		unit_multiplicator = 1 / 1000
	elif unit == "m":
		unit_multiplicator = 1000
	elif unit == "u":
		unit_multiplicator = 1000000
	elif unit == "n":
		unit_multiplicator = 1000000000

	cmap = None
	norm = None
	color = None
	if usedColorMap:
		# Définir une colormap du rouge au bleu
		colors = [(0/255, 255/255, 0/255), (89/255, 0/255, 140/255)]
		cmap = mcolors.LinearSegmentedColormap.from_list('cmap', colors, N=len(signals)-1)
		norm = mcolors.Normalize(vmin=0, vmax=len(signals)-2)

	# Parcourir les colonnes restantes pour tracer chaque série de données
	for idx, (name, signal) in enumerate(signals.items()):
		if name != axis_type and name != axis_name:
			if len(axis) == len(signal) and (name not in ignore_signals):
				if axis_type == "frequency":
					signal = np.abs(signal) / len(signal)
				if usedColorMap:
					color = cmap(norm(idx))
				if axis_form == "linear":
					ax.plot(axis * unit_multiplicator, signal, label=name, color=color)
				elif axis_form == "semilogx":
					ax.semilogx(axis * unit_multiplicator, signal, label=name, color=color)
				elif axis_form == "semilogy":
					ax.semilogy(axis * unit_multiplicator, signal, label=name, color=color)
				elif axis_form == "loglog":
					ax.loglog(axis * unit_multiplicator, signal, label=name, color=color)
				else:
					print("Error: axis_form not recognized")
					print("The possible values are: regular, semilogx, semilogy, loglog")
					return
				
	# Ajouter des étiquettes et une légende
	if plot_type == "time":
		ax.set_xlabel(f'Temps [{unit}s]', fontsize=axis_label_fontsize)
		ax.set_ylabel('Amplitude [V]', fontsize=axis_label_fontsize)
	elif plot_type == "frequency":
		ax.set_xlabel(f'Fréquence [{unit}Hz]', fontsize=axis_label_fontsize)
		ax.set_ylabel(f'Amplitude [V]', fontsize=axis_label_fontsize)
	elif plot_type == "pointsPerPeriod":
		ax.set_xlabel('Nombre de périodes', fontsize=axis_label_fontsize)
		ax.set_ylabel('Amplitude [V]', fontsize=axis_label_fontsize)
	elif plot_type == "freq_radian":
		ax.set_xlabel(f'Fréquence [{unit}Hz]', fontsize=axis_label_fontsize)
		ax.set_ylabel('Phase [rad]', fontsize=axis_label_fontsize)
	elif plot_type == "time_radian":
		ax.set_xlabel(f'Temps [{unit}s]', fontsize=axis_label_fontsize)
		ax.set_ylabel('Phase [rad]', fontsize=axis_label_fontsize)
	elif plot_type == "loop":
		ax.set_xlabel(f'Nb de tests', fontsize=axis_label_fontsize)
		ax.set_ylabel('Temps [s]', fontsize=axis_label_fontsize)

	if title != "":
		ax.set_title(title, fontsize=title_fontsize)
	else:
		ax.set_title(f"Affichage du spectre", fontsize=title_fontsize)
	ax.legend()
	ax.grid(True)

	# Ajuster la taille de la police des numérotations des axes
	ax.tick_params(axis='both', which='major', labelsize=tick_label_fontsize)
	ax.tick_params(axis='both', which='minor', labelsize=tick_label_fontsize)

