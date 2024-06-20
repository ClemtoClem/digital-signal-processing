import pandas as pd
import matplotlib.pyplot as plt

# Chargement des données CSV
def plot_csv_data(file_path):
    # Charger le fichier CSV en utilisant pandas
    df = pd.read_csv(file_path)
    
    # Séparation de la première colonne comme l'axe des temps
    temps = df.iloc[:, 0]
    
    # Parcourir les colonnes restantes pour tracer chaque série de données
    for col in df.columns[1:]:
        plt.plot(temps, df[col], label=col)
    
    # Ajouter des étiquettes et une légende
    plt.xlabel('Temps')
    plt.ylabel('Amplitude')
    plt.title('Graphique à partir du fichier CSV')
    plt.legend()
    
    # Afficher le plot
    plt.show()

# Exemple d'utilisation
if __name__ == "__main__":
    # Spécifiez le chemin vers votre fichier CSV
    fichier_csv = "./data/test_PID.csv"
    
    # Appel de la fonction pour tracer les données
    plot_csv_data(fichier_csv)
