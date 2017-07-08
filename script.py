'''
python script for Erdorbit_random
author: Thibaut Voirand
'''

from tkinter import *
import math
import random

global pos, alpha, beta, a, e, i, RAAN, om, numberOfSteps, dt, dur

DRAWING_SIZE_FACTOR = 9 / 10
RAAN = 0.0 # right ascension of the ascending node
MU_PLANET = 398600.4418 # gravitational parameter
ROT_PLANET = - 360.0 / 86164.1004 # rotation velocity
NB_OF_SECONDS_PER_DAY = 86400.0

'''defining math, geometry and canvas related functions---------------------------------------------
'''

def fromOrbitalToCartesianCoordinates(a, e, i, RAAN, om, t, MU_PLANET):
    '''
	Converting from orbital parameters to cartesian coordinates
	- Inputs:
			a       	semi-major axis (km)
			e       	eccentricity (-)
			i       	inclination (deg)
			RAAN    	right ascension of the ascending node (deg)
			om      	argument of periapsis (deg)
			t       	time spent since passage at periapsis (s)
			MU_PLANET	gravitational parameter of the planet	(km3/s2)
	- Outputs:
			pos   	  position vector of the orbiting object (km)
    '''


    # converting angles from degrees to radians
    i = i * math.pi / 180
    RAAN = RAAN * math.pi / 180
    om = om * math.pi / 180

    # computing mean anomaly
    n = math.sqrt(MU_PLANET / math.pow(a, 3))
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

    return pos

def fromJ2000ToECEF(posJ2000, t, ROT_PLANET):
    '''
	Converting coordinates from standard J2000 to Earth-centered Earth-fixed reference frame
      - Inputs:
          posJ2000    position vector in rotating frame (km)
          t 		 	    time (s)
          ROT_PLANET	planet rotational velocity (deg/s)
      - Outputs:
          posECEF	    position vector in planet fixed ref. frame (km)
    '''

    ROT_PLANET = ROT_PLANET * math.pi / 180

    posECEF = [
        math.cos(ROT_PLANET * t) * posJ2000[0] + math.sin(ROT_PLANET * t) * posJ2000[1],
        -math.sin(ROT_PLANET * t) * posJ2000[0] + math.cos(ROT_PLANET * t) * posJ2000[1],
        posJ2000[2]
    ]

    return posECEF

def getExtremum(array):
    '''
	Getting extremum of a 2D array
    	- Inputs:
    	   array     2D array
    	- Outputs:
    	   extremum  extremum value of 2D array
    '''

    extremum = array[0][0]

    for j in range(len(array)):
        for k in range(len(array[j])):
            if (extremum < abs(array[j][k])):
                extremum = abs(array[j][k])

    return extremum

def multiplyArrayByScalar(inputArray, scalar):
    '''
	Multiplying each component of a 2D array by a scalar
    	- Inputs:
          inputArray  input array
          scal        input scalar
    	- Outputs:
    		  outputArray output array
    '''

    outputArray = []

    for j in range(len(inputArray)):
        outputArray.append(inputArray[j])
        for k in range(len(inputArray[j])):
            outputArray[j][k] = inputArray[j][k] * scalar

    return outputArray

def planarProjection(inputCoord, alpha, beta):
    '''
	Projecting an array of 3D coordinates on a plane
    - Inputs:
	     inputCoord 	 array of 3D coordinates
		 beta	         rotation of the projection plane around Oz
		 alpha         rotation of the p.p. around intersection of itself and Oxy plane
	- Outputs:
	   outputCoord   2D projection of the input array on the plane
    '''

    outputCoord = []

    for j in range(len(inputCoord)):
        outputCoord.append([
            math.cos(beta) * inputCoord[j][0] - math.sin(beta) * inputCoord[j][1],
            -math.sin(beta) * math.sin(alpha) * inputCoord[j][0] -
                math.cos(beta) * math.sin(alpha) * inputCoord[j][1] +
                math.cos(alpha) * inputCoord[j][2]
        ])

    return outputCoord

def computePositions(a, e, i, RAAN, om, MU_PLANET, ROT_PLANET, numberOfSteps, dt):
    '''
      This function computes the positions of the orbiting object in a ECEF reference frame
      - Inputs:
        - orbital parameters:
          a           semi-major axis
          e           eccentricity
          i           inclination
          RAAN        right ascension of the ascending node
          om          argument of periapsis
        - planetary constants:
          MU_PLANET   planet gravitational parameter
          ROT_PLANET  planet rotational velocity
        numberOfSteps number of steps
        dt            step time
      Ouputs:
        positionsArray array of positions
    '''

    positionsArray = [];

    for j in range(numberOfSteps):
        positionsArray.append(
            fromOrbitalToCartesianCoordinates(
                a,
                e,
                i,
                RAAN,
                om,
                j * dt,
                MU_PLANET
            )
        )
        positionsArray[j] = fromJ2000ToECEF(positionsArray[j], j * dt, ROT_PLANET)

    return positionsArray

def adaptCoordinatesToCanvasFrame(inputCoord, canvasWidth, canvasHeight):
    '''
      The reference frame of the canvas is centered in the upper left corner, with the X axis
      pointing to the right, and the Y axis pointing down.
      This function transforms 2D coordinates so that the drawing is centered at the center of the
      canvas, with the Y axis pointing up.
      - Input:
        inputCoord   array of 2D coordinates expressed in a "classical" cartesian reference frame
        canvasWidth  canvas width
        canvasHeight canvas height
      - Output:
        outputCoord  array of 2D coordinates, adapted to the HTML5 canvas' reference frame
    '''

    outputCoord = []

    for j in range(len(inputCoord)):
        outputCoord.append([
            inputCoord[j][0] + canvasWidth / 2,
            -inputCoord[j][1] + canvasHeight / 2
        ])

    return outputCoord

def resizeDrawingToFitCanvas(inputCoord, canvasHeight):
    '''
    This function resizes the drawing to fit the canvas size
      - Input:
          inputCoord    input coordinates
          canvasHeight  canvas height
      - Output:
          outputCoord   output coordinates
    '''

    drawingSize = getExtremum(inputCoord) * 2

    outputCoord = multiplyArrayByScalar(
        inputCoord,
        canvasHeight / drawingSize * DRAWING_SIZE_FACTOR
    )

    return outputCoord

def drawOrbit(pos, alpha, beta, canvas, canvasWidth, canvasHeight):
    '''
      This function draws the trajectory of an orbiting object on the canvas
      This function calls several previously defined functions
      - Input:
          pos   array of 3D cartesian positions the orbiting object
          alpha         projection plane angle
          beta          projection plane angle
          canvas        canvas
          canvasWidth   canvas width
          canvasHeight  canvas height
    '''

    canvas.delete('all')

    pos2D = planarProjection(pos, alpha, beta)

    pos2D = adaptCoordinatesToCanvasFrame(pos2D, canvasWidth, canvasHeight)

    line = canvas.create_line(pos2D)

def getInputParameters():
    '''
    This function assigns parameters values entered by the user in the input forms to the
    calculation variables
    '''

    global pos, alpha, beta, a, e, i, RAAN, om, numberOfSteps, dt, dur

    # semi-major axis
    a = random.uniform(float(aMinInput.get()), float(aMaxInput.get()))
    # eccentricity
    e = random.uniform(float(eMinInput.get()), float(eMaxInput.get()))
    # inclination
    i = random.uniform(float(iMinInput.get()), float(iMaxInput.get()))
    # argument of periapsis
    om = random.uniform(float(omMinInput.get()), float(omMaxInput.get()))
    # duration
    dur = NB_OF_SECONDS_PER_DAY * random.uniform(float(durMinInput.get()), float(durMaxInput.get()))
    # number of steps
    numberOfSteps = int(numberOfStepsInput.get())
    # step time
    dt = dur / numberOfSteps

    # projection plane angles
    alpha = random.uniform(float(alphaMinInput.get()), float(alphaMaxInput.get())) * math.pi / 180
    beta = random.uniform(float(betaMinInput.get()), float(betaMaxInput.get())) * math.pi / 180

def displayInputParameters():
    '''
    This function displays the values used for the input parameters
    '''

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

'''defining user-interface related functions--------------------------------------------------------
'''

def computeButtonCallBack():
    '''
    Callback function of the compute button
    '''

    global pos, alpha, beta, a, e, i, RAAN, om, numberOfSteps, dt

    getInputParameters()

    pos = computePositions(a, e, i, RAAN, om, MU_PLANET, ROT_PLANET, numberOfSteps, dt)

    pos = resizeDrawingToFitCanvas(pos, canvasHeight)

    drawOrbit(pos, alpha, beta, canvas, canvasWidth, canvasHeight)

    displayInputParameters()

def upButtonCallBack():
    '''
    Callback function of the "up" button
    '''

    global pos, alpha, beta, a, e, i, RAAN, om, numberOfSteps, dt

    # increase projection plane angle
    alpha = alpha + 5 * math.pi / 180

    drawOrbit(pos, alpha, beta, canvas, canvasWidth, canvasHeight)

def downButtonCallBack():
    '''
    Callback function of the "down" button
    '''

    global pos, alpha, beta, a, e, i, RAAN, om, numberOfSteps, dt

    # increase projection plane angle
    alpha = alpha - 5 * math.pi / 180

    drawOrbit(pos, alpha, beta, canvas, canvasWidth, canvasHeight)

'''defining user-interface elements-----------------------------------------------------------------
'''

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
headerLabel = Label(header, text = 'Erdorbit_random', font = 'Helvetica 20 bold')
headerLabel.pack()

# creating and positioning simulation frame
simulationFrame = Frame(root)
simulationFrame.grid(row = 1, column = 0)

# creating and positioning footer
footer = Frame(root)
footer.grid(row = 2, column = 0)
footerLabel = Label(footer, text = 'Erdorbit_random v0.1.0')
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
aLabel = Label(orbitalParametersFrame, text = 'Semi-major axis (km, must be >0)')
aLabel.grid(row = 1, column = 0)
aMinInput = Entry(orbitalParametersFrame, width = 6)
aMinInput.insert(0, '6471')
aMinInput.grid(row = 1, column = 1)
aMaxInput = Entry(orbitalParametersFrame, width = 6)
aMaxInput.insert(0, '50000')
aMaxInput.grid(row = 1, column = 2)

# eccentricity
eLabel = Label(orbitalParametersFrame, text = 'Eccentricity (must be between 0 and 1)')
eLabel.grid(row = 2, column = 0)
eMinInput = Entry(orbitalParametersFrame, width = 6)
eMinInput.insert(0, '0')
eMinInput.grid(row = 2, column = 1)
eMaxInput = Entry(orbitalParametersFrame, width = 6)
eMaxInput.insert(0, '0.99')
eMaxInput.grid(row = 2, column = 2)

# inclination
iLabel = Label(orbitalParametersFrame, text = 'Inclination (deg)')
iLabel.grid(row = 3, column = 0)
iMinInput = Entry(orbitalParametersFrame, width = 6)
iMinInput.insert(0, '0')
iMinInput.grid(row = 3, column = 1)
iMaxInput = Entry(orbitalParametersFrame, width = 6)
iMaxInput.insert(0, '360')
iMaxInput.grid(row = 3, column = 2)

# argument of periapsis
omLabel = Label(orbitalParametersFrame, text = 'Argument of periapsis (deg)')
omLabel.grid(row = 4, column = 0)
omMinInput = Entry(orbitalParametersFrame, width = 6)
omMinInput.insert(0, '0')
omMinInput.grid(row = 4, column = 1)
omMaxInput = Entry(orbitalParametersFrame, width = 6)
omMaxInput.insert(0, '360')
omMaxInput.grid(row = 4, column = 2)

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
durLabel = Label(simulationParametersFrame, text = 'Duration (days)')
durLabel.grid(row = 1, column = 0)
durMinInput = Entry(simulationParametersFrame, width = 6)
durMinInput.insert(0, '0.1')
durMinInput.grid(row = 1, column = 1)
durMaxInput = Entry(simulationParametersFrame, width = 6)
durMaxInput.insert(0, '100')
durMaxInput.grid(row = 1, column = 2)

# number of steps
numberOfStepsLabel = Label(simulationParametersFrame, text = 'Number of steps (-)')
numberOfStepsLabel.grid(row = 2, column = 0)
numberOfStepsInput = Entry(simulationParametersFrame, width = 6)
numberOfStepsInput.insert(0, '2000')
numberOfStepsInput.grid(row = 2, column = 1)

#plane angles
alphaLabel = Label(simulationParametersFrame, text='Projection angle alpha (deg)')
alphaLabel.grid(row = 3, column = 0)
alphaMinInput = Entry(simulationParametersFrame, width = 6)
alphaMinInput.insert(0, '0')
alphaMinInput.grid(row = 3, column = 1)
alphaMaxInput = Entry(simulationParametersFrame, width = 6)
alphaMaxInput.insert(0, '90')
alphaMaxInput.grid(row = 3, column = 2)
betaLabel = Label(simulationParametersFrame, text='Projection angle beta (deg)')
betaLabel.grid(row = 4, column = 0)
betaMinInput = Entry(simulationParametersFrame, width = 6)
betaMinInput.insert(0, '0')
betaMinInput.grid(row = 4, column = 1)
betaMaxInput = Entry(simulationParametersFrame, width = 6)
betaMaxInput.insert(0, '90')
betaMaxInput.grid(row = 4, column = 2)

# creating results frame
resultsFrame = Frame(simulationFrame)
resultsFrame.grid(row = 1, column = 0)
resultsLabel = Label(resultsFrame, text = 'Results', font = 'Helvetica 16 bold')
resultsLabel.grid(row = 0, column = 0)

# creating canvas
canvasWidth = 300
canvasHeight = 300
canvas = Canvas(resultsFrame, width = str(canvasWidth), height = str(canvasHeight))
canvas.grid(row = 2, column = 0)

# creating 'compute' button
computeButton = Button(resultsFrame, text='compute', command = computeButtonCallBack)
computeButton.grid(row = 1, column = 0)

# creating 'up' button
upButton = Button(resultsFrame, text='up', command = upButtonCallBack)
upButton.grid(row = 1, column = 1)

# creating 'down' button
downButton = Button(resultsFrame, text='down', command = downButtonCallBack)
downButton.grid(row = 1, column = 2)

root.mainloop()
