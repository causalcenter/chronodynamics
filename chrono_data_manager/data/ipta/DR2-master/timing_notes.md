## Some notes for setting up parfiles

Use TEMPO2, but try to keep it backwards compatible as much as
possible with TEMPO
  * UNITS(TDB) is a crucial one (note that IPTA DR1 parfiles are in TCB!)
  * Use all explicit flags (see NANOGrav parfiles, for example)
```
SOLARN0               0.00
EPHEM               DE421
ECL                 IERS2010
CLK                 TT(BIPM)    
UNITS               TDB
TIMEEPH             FB90
T2CMETHOD           TEMPO
CORRECT_TROPOSPHERE N
PLANET_SHAPIRO      N
DILATEFREQ          N
```

Since IPTA DR1 parfiles are in TCB format, the command to convert them is:
```
> tempo2 -gr transform parFile outputFile back
```
Where "back" is used to go from TCB to TDB, which is what we want.
Note that the remaining file will have the values converted, but it will still
say ```UNITS TCB```.  So just change that to TDB (and if T2 binary model is 
used, change that to DD as well).

**Use the latest version of TEMPO2 from here:  https://bitbucket.org/psrsoft/tempo2**

## More detailed processing notes

  - Assemble appropriate .tim files using the ```assemble_tims.py``` script
  - As an initial parfile, if NANoGrav data are involved, use the NANOGrav 9yr parfile, but remove the ECORR params and DMX lines (keep T2EFAC and T2EQUAD)
  - Use the python script ```check_groups_and_JUMPSs.py``` to find the required JUMPs (including the reference data set)
  - Add PPTA fixed jumps if any PPTA data are present ```utils/PPTA_fixed_JUMPs.list```
  - Use autoDM plugging in tempo2 to find the DMMODEL mjd grid 
    - e.g.: ```tempo2 -gr autoDM -f parfile timfile -dt 60```
	- Here, I used a 60d (-dt 60) grid. This will check whether there are multi-frequency data available for these grid points, and send out a new par file called (XXX_dm.par) with two additional files, for the DM and common-mode (CM) signal, respestively. Note that the CM signal is not useful for us. This new parfile has the DM grid values in it.
  - Run tempo2 with this new par file and then you can be able to fit the DMMODEL.
    - **NOTE:** make sure you run tempo2 with ```-qrfit``` option, otherwise you will not get a stable solution! (we discussed this before)
  - Now the fit will modify the .dm file, which you can use gnuplot to see the DM variation.
    - This will give a better DM grid, but you may need to tweak some grid points if you wish.
    - Running autoDM may take a while for some pulsars

