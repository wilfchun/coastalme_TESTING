# coastalme
CoastalME (Coastal Modelling Environment) simulates the long-term behaviour of a coast. This initial version is a prototype which considers simple soft cliff cross-shore effects only.

By David Favis-Mortlock (Environmental Change Institute, University of Oxford) and Andres Payo (British Geological Survey)

See <a href="https://github.com/coastalme/coastalme" target="_blank">https://github.com/coastalme/coastalme</a> for the latest version of the source code.

INSTRUCTIONS

1. Create a local copy of the github repository, for example by downloading a zipfile, then unpacking it.

2. At a command-line prompt, change to the coastalme-master folder, then to the src folder.

3. If using Linux: copy run_cmake.sh.LINUX to run_cmake.sh, then run run_cmake.sh. If using Windows: copy run_cmake.bat.MSCV2013 to run_cmake.bat, then run run_cmake.bat. If you see error messages re. missing software (if, for example, it tells you that CMake cannot be found, or GDAL cannot be found) then you need to install the missing software.

4. Run make install. This should create an executable file called cme (on Linux) or cme.exe (on Windows) in the coastalme-master folder.

5. Edit cme.ini to tell CoastalME which input file you wish to use (for example, in/CliffFineBays/CliffFineBays.dat)

6. Run cme (or cme.exe on Windows). The output will appear in the out/ folder.

7. If using QGIS to visualize the output, look for a QGIS project file (e.g. CliffFineBays.qgs) somewhere in the in/ folder.

8. Enjoy!

Dave Favis-Mortlock and Andres Payo






