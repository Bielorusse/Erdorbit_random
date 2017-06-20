'''
python script for Erdorbit_random
author: Thibaut Voirand
'''

# importing Tkinter GUI framework
from tkinter import *
import math
import random

# declaring global variables
global pos, alpha, beta, delta

# defining functions

# defining conversion from orbital parameters to cartesian coordinates
def orb2car(a, e, i, RAAN, om, t, mu_planet):

    # converting angles from degrees to radians
    i = i * math.pi / 180
    RAAN = RAAN * math.pi / 180
    om = om * math.pi / 180

    # computing mean anomaly
    n = math.sqrt(mu_planet / math.pow(a, 3))
    M = n * t

    # computing eccentric anomaly
    E = [M]
    for j in range(100):
        E.append(E[j] + (M - E[j] + e * math.sin(E[j])) / (1 - e * math.cos(E[j])))
        if(abs(E[j+1] - E[j]) < 1e-8):
            E = E[j+1]
            break

    # computing true anomaly
    nu = 2 * math.atan2(
            math.sqrt(1 + e) * math.sin(E / 2),
            math.sqrt(1 - e) * math.cos(E / 2)
        ) % (math.pi * 2)

    # computing radius
    r = a * (1 -math.pow(e, 2)) / (1 + e * math.cos(nu))

    # computing position vector
    pos = [
        r *
            (math.cos(om + nu) * math.cos(RAAN) -
                math.sin(om + nu) * math.sin(RAAN) * math.cos(i)),
        r *
            (math.cos(om + nu) * math.sin(RAAN) -
                math.sin(om + nu) * math.cos(RAAN) * math.cos(i)),
        r * (math.sin(om + nu) * math.sin(i))
    ]

    # returning result
    return pos

# defining projection of 3D vector on a plane
def proj(vec, alpha, beta, delta):

    vec = [
        math.cos(beta) * vec[0] - math.sin(beta) * vec[1],
        -math.sin(beta) * math.sin(alpha) * vec[0] -
            math.cos(beta) * math.sin(alpha) * vec[1] +
            math.cos(alpha) * vec[2]
    ]

    vec = [
        math.cos(delta) * vec[0] + math.sin(delta) * vec[1],
        math.cos(delta) * vec[1] - math.sin(delta) * vec[0]
    ]

    return vec

# rotation reference frame
def rot_frame(vec, t, rot_planet):

    # converting rotation velocity from deg to radians
    rot_planet = rot_planet * math.pi / 180

    vec_rot = [
        math.cos(rot_planet * t) * vec[0] + math.sin(rot_planet * t) * vec[1],
        -math.sin(rot_planet * t) * vec[0] + math.cos(rot_planet * t) * vec[1],
        vec[2]
    ]

    return vec_rot

# finding extremum in an array
def getextr(array):

    extr = array[0][0]

    for j in range(len(array)):
        for k in range(len(array[j])):
            if (extr < abs(array[j][k])):
                extr = abs(array[j][k])

    return extr

# multiplying each cell of an array by a scalar
def array_scal(array, scal):

    array_output = []

    for j in range(len(array)):
        array_output.append(array[j])
        for k in range(len(array[j])):
            array_output[j][k] = array[j][k] * scal

    return array_output

# callback function of the compute button
def computeCallBack():

    # declaring global variables
    global pos, alpha, beta, delta

    # obtaining parameters

    # semi-major axis
    a = random.uniform(float(aMin_input.get()), float(aMax_input.get()))
    # eccentricity
    e = random.uniform(float(eMin_input.get()), float(eMax_input.get()))
    # inclination
    i = random.uniform(float(iMin_input.get()), float(iMax_input.get()))
    # right ascension of the ascending node
    RAAN = 0.0
    # argument of periapsis
    om = random.uniform(float(omMin_input.get()), float(omMax_input.get()))
    # duration
    dur = 86400.0 * random.uniform(float(durMin_input.get()), float(durMax_input.get()))
    # number of steps
    step_nb = int(step_nb_input.get())
    # step time
    step_time = dur / step_nb

    # planet parameters
    mu_planet = 398600.4418 # gravitational parameter
    rot_planet = 360.0 / 86164.1004 # rotation velocity

    # plane angles
    alpha = random.uniform(float(alphaMin_input.get()), float(alphaMax_input.get())) * math.pi / 180
    beta = random.uniform(float(betaMin_input.get()), float(betaMax_input.get())) * math.pi / 180
    delta = random.uniform(float(deltaMin_input.get()), float(deltaMax_input.get())) * math.pi / 180

    # computing cartesian positions from orbital parameters for each step
    pos = []
    for j in range(step_nb):
        pos.append(orb2car(a, e, i, RAAN, om, j * step_time, mu_planet))
        pos[j] = rot_frame(pos[j], j * step_time, rot_planet)

    # getting size of drawing
    size = getextr(pos) * 2

    # resizing positions to fit the graph
    pos = array_scal(pos, (canvas_height / size * 9 / 10))

    # projection 3D positions on a 2D plane
    pos_2d = []
    for j in range(len(pos)):
        pos_2d.append(proj(pos[j], alpha, beta, delta))

    # converting coordinates to fit the canvas reference frame
    pos_canvas = []
    for j in range(len(pos_2d)):
        pos_2d[j][0] = pos_2d[j][0] + canvas_height / 2
        pos_2d[j][1] = -pos_2d[j][1] + canvas_height / 2
        for k in range(len(pos_2d[j])):
            pos_canvas.append(pos_2d[j][k])

    # clearing canvas
    canvas.delete('all')

    # drawing positions on canvas
    line = canvas.create_line(pos_canvas)

    # displaying parameters
    randomParametersFrame = Frame(resultsFrame)
    randomParametersFrame.grid(row = 2, column = 1)
    aRandom_display = Label(randomParametersFrame, text='a = ' + str(a))
    aRandom_display.grid(row = 0, column = 0)
    eRandom_display = Label(randomParametersFrame, text='e = ' + str(e))
    eRandom_display.grid(row = 1, column = 0)
    iRandom_display = Label(randomParametersFrame, text='i = ' + str(i))
    iRandom_display.grid(row = 2, column = 0)
    omRandom_display = Label(randomParametersFrame, text='om = ' + str(om))
    omRandom_display.grid(row = 3, column = 0)
    durRandom_display = Label(randomParametersFrame, text='dur = ' + str(dur / 86400))
    durRandom_display.grid(row = 4, column = 0)

# callback function of up button
def upCallBack():

    # declaring global variables
    global pos, alpha, beta, delta

    # increase projection plane angle
    alpha = alpha + 5 * math.pi / 180

    # projection 3D positions on a 2D plane
    pos_2d = []
    for j in range(len(pos)):
        pos_2d.append(proj(pos[j], alpha, beta, delta))

    # converting coordinates to fit the canvas reference frame
    pos_canvas = []
    for j in range(len(pos_2d)):
        pos_2d[j][0] = pos_2d[j][0] + canvas_height / 2
        pos_2d[j][1] = -pos_2d[j][1] + canvas_height / 2
        for k in range(len(pos_2d[j])):
            pos_canvas.append(pos_2d[j][k])

    # clearing canvas
    canvas.delete('all')

    # drawing positions on canvas
    line = canvas.create_line(pos_canvas)

# creating, dimensioning, and centering root window
root = Tk()
root.geometry('%dx%d+%d+%d' %
    (800,
    600,
    (root.winfo_screenwidth()/2)-(800/2),
    (root.winfo_screenheight()/2)-(600/2)))

# creating  and positioning header
header = Frame(root)
header.grid(row = 0, column = 0)
headerLabel = Label(header, text = 'Erdorbit', font = 'Helvetica 20 bold')
headerLabel.pack()

# creating and positioning simulation frame
simulationFrame = Frame(root)
simulationFrame.grid(row = 1, column = 0)

# creating and positioning footer
footer = Frame(root)
footer.grid(row = 2, column = 0)
footerLabel = Label(footer, text = 'Erdorbit_random v0.0.2')
footerLabel.pack()

# creating parameters input frame

# orbital parameters
orbitalParametersFrame = Frame(simulationFrame)
orbitalParametersFrame.grid(row = 0, column = 0)
orbitalParametersLabel = Label(orbitalParametersFrame,
        text = 'Orbital Parameters',
        font = 'Helvetica 16 bold'
    )
orbitalParametersLabel.grid(row = 0, column = 0)
orbitalMinLabel = Label(orbitalParametersFrame, text = 'min')
orbitalMinLabel.grid(row = 0, column = 1)
orbitalMaxLabel = Label(orbitalParametersFrame, text = 'max')
orbitalMaxLabel.grid(row = 0, column = 2)

# semi-major axis
a_label = Label(orbitalParametersFrame, text = 'Semi-major axis (km, must be >0)')
a_label.grid(row = 1, column = 0)
aMin_input = Entry(orbitalParametersFrame, width = 6)
aMin_input.insert(0, '6471')
aMin_input.grid(row = 1, column = 1)
aMax_input = Entry(orbitalParametersFrame, width = 6)
aMax_input.insert(0, '50000')
aMax_input.grid(row = 1, column = 2)

# eccentricity
e_label = Label(orbitalParametersFrame, text = 'Eccentricity (must be between 0 and 1)')
e_label.grid(row = 2, column = 0)
eMin_input = Entry(orbitalParametersFrame, width = 6)
eMin_input.insert(0, '0')
eMin_input.grid(row = 2, column = 1)
eMax_input = Entry(orbitalParametersFrame, width = 6)
eMax_input.insert(0, '0.99')
eMax_input.grid(row = 2, column = 2)

# inclination
i_label = Label(orbitalParametersFrame, text = 'Inclination (deg)')
i_label.grid(row = 3, column = 0)
iMin_input = Entry(orbitalParametersFrame, width = 6)
iMin_input.insert(0, '0')
iMin_input.grid(row = 3, column = 1)
iMax_input = Entry(orbitalParametersFrame, width = 6)
iMax_input.insert(0, '360')
iMax_input.grid(row = 3, column = 2)

# argument of periapsis
om_label = Label(orbitalParametersFrame, text = 'Argument of periapsis (deg)')
om_label.grid(row = 4, column = 0)
omMin_input = Entry(orbitalParametersFrame, width = 6)
omMin_input.insert(0, '0')
omMin_input.grid(row = 4, column = 1)
omMax_input = Entry(orbitalParametersFrame, width = 6)
omMax_input.insert(0, '360')
omMax_input.grid(row = 4, column = 2)

# Simulation parameters
simulationParametersFrame = Frame(simulationFrame)
simulationParametersFrame.grid(row = 0, column = 1)
simulationParametersLabel = Label(simulationParametersFrame,
        text = 'Simulation Parameters',
        font = 'Helvetica 16 bold'
    )
simulationParametersLabel.grid(row = 0, column = 0)
simulationMinLabel = Label(simulationParametersFrame, text = 'min')
simulationMinLabel.grid(row = 0, column = 1)
simulationMaxLabel = Label(simulationParametersFrame, text = 'max')
simulationMaxLabel.grid(row = 0, column = 2)

# duration
dur_label = Label(simulationParametersFrame, text = 'Duration (days)')
dur_label.grid(row = 1, column = 0)
durMin_input = Entry(simulationParametersFrame, width = 6)
durMin_input.insert(0, '0.1')
durMin_input.grid(row = 1, column = 1)
durMax_input = Entry(simulationParametersFrame, width = 6)
durMax_input.insert(0, '100')
durMax_input.grid(row = 1, column = 2)

# number of steps
step_nb_label = Label(simulationParametersFrame, text = 'Number of steps (-)')
step_nb_label.grid(row = 2, column = 0)
step_nb_input = Entry(simulationParametersFrame, width = 6)
step_nb_input.insert(0, '30000')
step_nb_input.grid(row = 2, column = 1)

#plane angles
alpha_label = Label(simulationParametersFrame, text='Projection angle alpha (deg)')
alpha_label.grid(row = 3, column = 0)
alphaMin_input = Entry(simulationParametersFrame, width = 6)
alphaMin_input.insert(0, '0')
alphaMin_input.grid(row = 3, column = 1)
alphaMax_input = Entry(simulationParametersFrame, width = 6)
alphaMax_input.insert(0, '90')
alphaMax_input.grid(row = 3, column = 2)
beta_label = Label(simulationParametersFrame, text='Projection angle beta (deg)')
beta_label.grid(row = 4, column = 0)
betaMin_input = Entry(simulationParametersFrame, width = 6)
betaMin_input.insert(0, '0')
betaMin_input.grid(row = 4, column = 1)
betaMax_input = Entry(simulationParametersFrame, width = 6)
betaMax_input.insert(0, '90')
betaMax_input.grid(row = 4, column = 2)
delta_label = Label(simulationParametersFrame, text='Projection angle delta (deg)')
delta_label.grid(row = 5, column = 0)
deltaMin_input = Entry(simulationParametersFrame, width = 6)
deltaMin_input.insert(0, '0')
deltaMin_input.grid(row = 5, column = 1)
deltaMax_input = Entry(simulationParametersFrame, width = 6)
deltaMax_input.insert(0, '0')
deltaMax_input.grid(row = 5, column = 2)

# creating results frame
resultsFrame = Frame(simulationFrame)
resultsFrame.grid(row = 1, column = 0)
resultsLabel = Label(resultsFrame, text = 'Results', font = 'Helvetica 16 bold')
resultsLabel.grid(row = 0, column = 0)

# creating canvas
canvas_width = 300
canvas_height = 300
canvas = Canvas(resultsFrame, width = str(canvas_width), height = str(canvas_height))
canvas.grid(row = 2, column = 0)

# creating 'compute' button
button = Button(resultsFrame, text='compute', command = computeCallBack)
button.grid(row = 1, column = 0)

# creating 'up' button
upButton = Button(resultsFrame, text='up', command = upCallBack)
upButton.grid(row = 1, column = 1)

root.mainloop()
