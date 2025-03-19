# Generate_Plot_times.py 
# This python script can be used to create the target_plot_times input file
# used by CO2GraVISim.

import numpy as np

#This script builds up the plot times by specifying particular time intervals (e.g. 
#between the times when the injection flux changes), and the corresponding number of 
#plot times to create within each interval

#Set the key times between which to linearly space plot times in t_intervals.
#In nplot_intervals, specify the number of plot times generated within each interval.
#These should each count the first and last element of the interval. The length of 
#nplot_intervals should be one less than that of t_intervals.
# e.g. t_intervals = [0.,1.,3.,10.] and nplot_intervals = [3,3,2] will generate
# t = [ 0., 0.5, 1.0, 2.0, 3.0, 10. ]

# t_intervals     = [ 0., 1., 3., 4., 10.]
# nplot_intervals = [   5,  5,  5,  8    ] 

# t_intervals     = [ 0., 1., 2., 3., 6., 7., 10.]
# nplot_intervals = [   5,  3,  5,  3,  5,  4    ] 

# t_intervals     = [ 0., 676.463, 266557.21 ]
# # t_intervals     = [ 0., 676.463, 30218. ]
# nplot_intervals = [   11, 10   ] 

t_intervals = [ 0., 2.] #[ 0., 5.]
nplot_intervals = [ 20 ]


#Create the first interval
Plot_times = np.linspace(t_intervals[0], t_intervals[1], nplot_intervals[0])

#Build the subsequent intervals and add them onto what exists already,
#skipping the first element that would otherwise be repeated
for i,n in enumerate(nplot_intervals[1:]):
    t_new = np.linspace(t_intervals[i+1],t_intervals[i+2],n)
    Plot_times = np.hstack((Plot_times,t_new[1:]))

n_plot = len(Plot_times)


# Save the generated plot times to file

text_file = open("./Input/target_plot_times.txt", "w")
# First line: write comment about input structure
text_file.write("-- number of plot times followed by the plot times, each on a new row -- \n")
# Second line: number of plot times
text_file.write("{0} \n".format(n_plot))
# Remaining lines: Plot times on individual lines
np.savetxt(text_file, Plot_times, fmt='%.8f')
text_file.close()
