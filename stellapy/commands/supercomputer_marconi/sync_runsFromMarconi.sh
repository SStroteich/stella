#!/bin/bash                                                                                                                                              

# ==================================================                                                                                                     
# Script to synchronize the stella version that is on xula to marconi                                                                                    
# ==================================================                                                                                                     

# Location of the runs on the local computer                                                                                                  
here=$RUNS

# Location of the runs on MARCONI                                                                                                                 
there=$RUNS_MARCONI

# Output message to show the process                                                                                                                     
echo "Syncronization between the directories:"
echo " "$here
echo " "$there

# Syncronize the stella versions                                                                                                                         
rsync -avzh -L --progress --exclude 'not_important*' --exclude 'Finished' --exclude 'stella' --exclude '*.nc' --exclude 'restart*' --exclude '*.out.nc'  --exclude '*.error' --exclude '*.vmec_geo' --exclude '*.jacob' --exclude '*.dat' --exclude '*.scratch' --exclude 'mat' $there $here

exit

