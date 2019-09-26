This repository contains three C-codes to create histograms.

dist_bin.c creates a histogram of the distance between two groups of atoms.
water_bin.c creates a histogram of the number of waters in the groove of nhTMEM16
lipid_tail_bin.c creats a histogram of the lipid tail atom count in the groove of nhTMEM16

dist_bin.c and water_bin.c are identical except of the bin size used. Currently, for dist_bin.c, the bin size is set to 0.2A; for water_bin.c the bin size is 2. To run these routines, first use these command to compile the codes:
gcc -lm -O2 -o dist_bin dist_bin.c
gcc -lm -O2 -o water_bin water_bin.c
These will create executables dist_bin and water_bin. To run these executables, use will need a parameters_bin.dat file which contains number of lines in the data that needs to be binned. Then the binning can be done by executing:

./dist_bin parameters_bin.dat distance_313_432.dat distance_313_432_binned.dat

./water_bin parameters_bin.dat waters_groove.dat waters_groove_binned.dat

In the above, distance_313_432.dat is the data file containing distance values between residues 313 and 432 in nhTMEM16, and water_groove.dat contains water counts in the groove of nhTMEM16. The binned data will be in distance_313_432_binned.dat and waters_groove_binned.dat files. In these binned files the first column contains bin intervals, the second - normalized probabilities, and the third - counts. 

To bin the lipid tail atom count data, we follow the same set of commands:
gcc -lm -O2 -o lipid_tail_bin lipid_tail_bin.c
./lipid_tail_bin parameters_bin.dat lipid_tail_groove.dat lipid_tail_groove_binned.dat
lipid_tail_groove.dat file contains sequence of counts of POPE and POPG lipid tail atom count in the groove, arranged in two columns.  


## Repositories

