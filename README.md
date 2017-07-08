Erdorbit_random V0.1.0 <br>
Date: 07.07.2017 <br>
Author: Thibaut Voirand <br>

Uses the same principle as Erdorbit (see erdorbit.com). <br>
Generates random orbital parameters, and draws the corresponding orbit in a ECEF frame (Earth fixed, Earth centered). <br>

I basically rewrote the Erdorbit script in Python3, to get to know the language, and added a random parameters generation feature, to find new ideas of interesting orbital parameters. <br>

Added in this version:<br>
    - refactoring:<br>
        - names improved: more explicit and use CamelCase
        - functions added for more abstraction
            - compute positions
            - adaptCoordinatesToCanvasFrame
            - resizeDrawingToFitCanvas
            - drawOrbit
            - getInputParameters
            - displayInputParameters
