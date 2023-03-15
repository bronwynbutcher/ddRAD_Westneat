##Run 1
mkdir denovo_run1
denovo_map.pl --samples /workdir/bgb27/demultfilter --popmap /workdir/bgb27/popmap_bushtit.txt -m 3 -M 1 -n 1 -T 15 -o /workdir/bgb27/denovo_run1
populations -P /workdir/bgb27/denovo_run1 -M /workdir/bgb27/popmap_bushtit.txt -t 15 -r 0.8

# note -T or -t = # of threads/CPUs to use
# for populations: -P = path to directory containing the Stacks files (output from the denovo run)
#                  -M = path to the population map (in thsi case same as used in the denovo trial)
#                  -r = minimum percentage of individuals in a population req to have a locus to process the locus

##Run 2
mkdir denovo_run2
denovo_map.pl --samples /workdir/bgb27/demultfilter --popmap /workdir/bgb27/popmap_bushtit.txt -m 3 -M 2 -n 2 -T 15 -o /workdir/bgb27/denovo_run2
populations -P /workdir/bgb27/denovo_run2 -M /workdir/bgb27/popmap_bushtit.txt -t 15 -r 0.8

##Run 3
mkdir denovo_run3
denovo_map.pl --samples /workdir/bgb27/demultfilter --popmap /workdir/bgb27/popmap_bushtit.txt -m 3 -M 3 -n 3 -T 15 -o /workdir/bgb27/denovo_run3
populations -P /workdir/bgb27/denovo_run3 -M /workdir/bgb27/popmap_bushtit.txt -t 15 -r 0.8

##Run 4
mkdir denovo_run4
denovo_map.pl --samples /workdir/bgb27/demultfilter --popmap /workdir/bgb27/popmap_bushtit.txt -m 3 -M 4 -n 4 -T 15 -o /workdir/bgb27/denovo_run4
populations -P /workdir/bgb27/denovo_run4 -M /workdir/bgb27/popmap_bushtit.txt -t 15 -r 0.8

##Run 5
mkdir denovo_run5
denovo_map.pl --samples /workdir/bgb27/demultfilter --popmap /workdir/bgb27/popmap_bushtit.txt -m 3 -M 5 -n 5 -T 15 -o /workdir/bgb27/denovo_run5
populations -P /workdir/bgb27/denovo_run5 -M /workdir/bgb27/popmap_bushtit.txt -t 15 -r 0.8

##Run 6
mkdir denovo_run6
denovo_map.pl --samples /workdir/bgb27/demultfilter --popmap /workdir/bgb27/popmap_bushtit.txt -m 3 -M 6 -n 6 -T 15 -o /workdir/bgb27/denovo_run6
populations -P /workdir/bgb27/denovo_run6 -M /workdir/bgb27/popmap_bushtit.txt -t 15 -r 0.8

##Run 7
mkdir denovo_run7
denovo_map.pl --samples /workdir/bgb27/demultfilter --popmap /workdir/bgb27/popmap_bushtit.txt -m 3 -M 7 -n 7 -T 15 -o /workdir/bgb27/denovo_run7
populations -P /workdir/bgb27/denovo_run7 -M /workdir/bgb27/popmap_bushtit.txt -t 15 -r 0.8

##Run 8
mkdir denovo_run8
denovo_map.pl --samples /workdir/bgb27/demultfilter --popmap /workdir/bgb27/popmap_bushtit.txt -m 3 -M 8 -n 8 -T 15 -o /workdir/bgb27/denovo_run8
populations -P /workdir/bgb27/denovo_run8 -M /workdir/bgb27/popmap_bushtit.txt -t 15 -r 0.8

##Run 9
mkdir denovo_run9
denovo_map.pl --samples /workdir/bgb27/demultfilter --popmap /workdir/bgb27/popmap_bushtit.txt -m 3 -M 9 -n 9 -T 15 -o /workdir/bgb27/denovo_run9
populations -P /workdir/bgb27/denovo_run9 -M /workdir/bgb27/popmap_bushtit.txt -t 15 -r 0.8