# Generate_Plot_times.py
# This python script can be used to create the target_plot_times input file
# used by CO2GraVISim.

# This script builds up the plot times by specifying particular time intervals (e.g.
# between the times when the injection flux changes), and the corresponding number of
# plot times to create within each interval

# Set the key times between which to linearly space plot times in t_intervals.
# In nplot_intervals, specify the number of plot times generated within each interval.
# These should each count the first and last element of the interval. The length of
# nplot_intervals should be one less than that of t_intervals.
# e.g. t_intervals = [0.,1.,3.,10.] and nplot_intervals = [3,3,2] will generate
# t = [ 0., 0.5, 1.0, 2.0, 3.0, 10. ]


import numpy as np
import argparse
from pathlib import Path


def main():
    # Parse command-line input arguments
    parser = argparse.ArgumentParser(
        description="Target folder, Interval endtimes, and number of plots per interval."
    )
    parser.add_argument(
        "--input", type=str, default="./Input/", help="Input data folder path"
    )
    parser.add_argument(
        "--times",
        type=float,
        nargs="+",
        default=[0, 10],
        help="specified plot interval times",
    )
    parser.add_argument(
        "--intervals",
        type=int,
        nargs="+",
        default=[20],
        help="number of plot times within each interval",
    )
    # Flag for whether to print output to command line
    parser.add_argument("--print", action="store_true")

    args = parser.parse_args()

    # Extract the input and output folder paths, and remove any trailing quotation marks and slashes
    input_folder = Path(args.input)

    # Ensure that the folders exist
    if not input_folder.is_dir():
        print(f"Error: The input folder {input_folder} does not exist.")
        exit(1)

    print(f" Input folder : {input_folder}")

    Print_flag = args.print

    t_intervals = args.times
    nplot_intervals = args.intervals

    Plot_times, n_plot, dt = pt_generate(t_intervals, nplot_intervals)
    pt_save(input_folder, Print_flag, Plot_times, n_plot, dt)


def pt_generate(t_intervals, nplot_intervals):

    # Create the first interval
    Plot_times = np.linspace(t_intervals[0], t_intervals[1], nplot_intervals[0])

    # Build the subsequent intervals and add them onto what exists already,
    # skipping the first element that would otherwise be repeated
    for i, n in enumerate(nplot_intervals[1:]):
        t_new = np.linspace(t_intervals[i + 1], t_intervals[i + 2], n)
        Plot_times = np.hstack((Plot_times, t_new[1:]))

    n_plot = len(Plot_times)

    dt = np.zeros([len(Plot_times + 1)])
    dt[1:] = np.diff(Plot_times)

    return Plot_times, n_plot, dt

def pt_save(input_folder, Print_flag, Plot_times, n_plot, dt):

    if Print_flag:
        print(f"Plot number   Plot time   Time step")
        # print(f'{i:3d},      {t:.3f},   {dt[i]:.3f}')
        for i, t in enumerate(Plot_times):
            print(f"{i:3d}           {t:4.3f}       {dt[i]:4.3f}")

    # Save the generated plot times to file

    text_file = open( input_folder / "target_plot_times.txt", "w")
    # First line: write comment about input structure
    text_file.write(
        "-- number of plot times followed by the plot times, each on a new row -- \n"
    )
    # Second line: number of plot times
    text_file.write("{0} \n".format(n_plot))
    # Remaining lines: Plot times on individual lines
    np.savetxt(text_file, Plot_times, fmt="%.8f")
    text_file.close()


if __name__ == "__main__":
    main()
