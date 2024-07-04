import numpy
import matplotlib.pyplot as plt
from parseCSV import *

def inversion_bits(i,p): # p = nombre de bits
    c = "{0:0%db}"%p
    binaire = c.format(i)
    inv = binaire[::-1]
    return int(inv,2)

def fft(u,p):
    N=len(u)
    if N!=2**p:
        print("Erreur de taille")
        return
    A = numpy.zeros(N,dtype=numpy.complex64)
    B = numpy.zeros(N,dtype=numpy.complex64)
    for k in range(N):
        j = inversion_bits(k,p)
        A[j] = u[k]
    for q in range(1,p+1): # FFT à 2**q éléments
        taille = 2**q
        taille_precedente = 2**(q-1)
        nombre_tfd = 2**(p-q) # nombre de TFD à calculer
        phi = -1j*2*numpy.pi/taille
        for m in range(nombre_tfd):
            position = m*taille
            for i in range(taille_precedente):
                W = numpy.exp(phi*i)
                B[position+i] = A[position+i] + W*A[position+taille_precedente+i]
            for i in range(taille_precedente,taille):
                W = numpy.exp(phi*i)
                B[position+i] = A[position+i-taille_precedente] + W*A[position+i]
        (A,B)=(B,A) # échange des références des tableaux
    # normalise
    for k in range(N):
        A[k] = 2.0*numpy.abs(A[k])/N
    return A

if __name__ == "__main__":
    fichier_csv = "./data/test.csv"

    # Lecture des données à partir du fichier CSV
    signals: dict = read_csv(fichier_csv)

    time = signals["time"]
    signal = signals["signal(t)"]

    N = len(signal)
    P = numpy.int32(numpy.log2(N))
    fs = 1/time[1]
    freq = []
    for k in range(N):
        freq.append(fs*k/N)

    spectre = fft(signal,P)
            

    plt.figure(figsize=(10,4))
    plt.plot(time, signal,'r')
    plt.figure(figsize=(10,4))
    plt.loglog(freq, spectre,'r')
    plt.show()