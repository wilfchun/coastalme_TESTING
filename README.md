CoastalME (Coastal Modelling Environment) simulates the long-term behaviour of a coast. This initial version is a prototype which considers simple soft cliff cross-shore effects only.

By Andres Payo and Dave Favis-Mortlock (both British Geological Survey)

This is the TESTING version of CoastalME. See <a href="https://github.com/coastalme/coastalme" target="_blank">https://github.com/coastalme/coastalme</a> for the most recent release version of the source code.

INSTRUCTIONS

CoastalME builds easily using Linux. If you wish to run CoastalME on Windows, then we currently recommend using the Cygwin pseudo-Linux software to do this.

1. Create a local copy of the github repository, for example by downloading a zipfile, then unpacking it. For example, unpack it to /home/<YOUR NAME>/coast/CoastalME.

2. At a command-line prompt, change to the CoastalME folder, then to the src folder

3. If using Linux: copy run_cmake.sh.LINUX to run_cmake.sh OR if using Cygwin under Windows: copy run_cmake.CYGWIN to run_cmake.sh. Then run run_cmake.sh. If you see error messages re. missing software (for example, telling you that CMake cannot be found or is too old, or GDAL cannot be found or is too old) then you need to install or update the software that is causing the problem

4. Run make install. This will create an executable file called cme in the coastalme-master folder

5. Edit cme.ini to tell CoastalME which input file you wish to use (for example, in/CliffFineBays/CliffFineBays.dat)

6. Run cme. Output will appear in the out/ folder

8. Enjoy!

Dave Favis-Mortlock and Andres Payo

PS We are working on getting CoastalME to build on Windows using the Visual Studio compiler... watch this space!




