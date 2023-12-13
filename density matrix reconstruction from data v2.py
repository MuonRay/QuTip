# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 02:25:27 2022

@author: cosmi
"""

from matplotlib import pyplot
import qutip as q
import numpy as np
import matplotlib

measurement = [1.0, 0.5, 0.5]

p1, p2, p3 = measurement
density_matrx = np.array([[1-p1, p3 - .5 + 1.j * (p2-.5)],[p3 - .5 - 1.j * (p2-.5), p1]])

real = np.real(density_matrx)
imag = np.imag(density_matrx)

#make the figure
#this uses the qutip routine as the baseline but then makes modifications for bigger font sizes
fig = pyplot.figure()

labels = ['|D>','|S>']
zlim = [-1,1]
norm = matplotlib.colors.Normalize(zlim[0], zlim[1])
#first plot
ax = fig.add_subplot(121, projection = '3d')
ax.set_title('Real', fontsize = 30)
q.matrix_histogram(real,  limits = zlim, fig = fig, ax = ax, colorbar = False)
ax.set_xticklabels(labels, fontsize = 22)
ax.set_yticklabels(labels, fontsize = 22)
ax.tick_params(axis='z', labelsize=22)
cax, kw = matplotlib.colorbar.make_axes(ax, shrink=.75, pad=.0)
cb1 = matplotlib.colorbar.ColorbarBase(cax, norm = norm)
cl = pyplot.getp(cax, 'ymajorticklabels') 
pyplot.setp(cl, fontsize=22)
#next subplot
ax = fig.add_subplot(122, projection = '3d')
ax.set_title('Imaginary', fontsize = 30)
q.matrix_histogram(imag, limits = zlim, fig = fig, ax = ax, colorbar = False)
ax.set_xticklabels(labels, fontsize = 22)
ax.set_yticklabels(labels, fontsize = 22)
ax.tick_params(axis='z', labelsize=22)
cax, kw = matplotlib.colorbar.make_axes(ax, shrink=.75, pad=.0)



def errorBarSimple(trials, prob):
    #look at wiki http://en.wikipedia.org/wiki/Checking_whether_a_coin_is_fair
    '''returns 1 sigma error bar on each side i.e 1 sigma interval is val - err < val + err'''
    Z = 1.0
    s = np.sqrt(prob * (1.0 - prob) / float(trials))
    err = Z * s
    return err

err = errorBarSimple
trials = 100