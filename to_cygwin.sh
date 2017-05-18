#!/bin/sh

scp src/*.cpp 		davefm@Windows-PC:~/coast/CoastalME/src/
scp src/*.h   		davefm@Windows-PC:~/coast/CoastalME/src/
scp src/run_cmake* 	davefm@Windows-PC:~/coast/CoastalME/src/
scp README*		davefm@Windows-PC:~/coast/CoastalME/
scp src/*.txt		davefm@Windows-PC:~/coast/CoastalME/src/
scp src/cmake/Modules/* davefm@Windows-PC:~/coast/CoastalME/src/cmake/Modules/
scp COPY*		davefm@Windows-PC:~/coast/CoastalME/
scp src/do_cppcheck.sh	davefm@Windows-PC:~/coast/CoastalME/src/
scp cme.ini             davefm@Windows-PC:~/coast/CoastalME/
scp -r cshore/*         davefm@Windows-PC:~/coast/CoastalME/cshore/
#scp -r in/*             davefm@Windows-PC:~/coast/CoastalME/in/
scp -r scape/*          davefm@Windows-PC:~/coast/CoastalME/scape/

