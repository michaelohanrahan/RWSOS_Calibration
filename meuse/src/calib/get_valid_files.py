from glob import glob
import os
import argparse
import traceback

def main(basin, level):
    os.chdir(fr"p:\11209265-grade2023\wflow\RWSOS_Calibration\{basin}")
    # Get the list of directories containing run.done
    run_done_files = glob("data/2-interim/calib_data/{level}/**/run.done", recursive=True)
    run_done_dirs = {os.path.dirname(file) for file in run_done_files}

    # Get the list of output_scalar.nc files
    output_files = glob("data/2-interim/calib_data/{level}/**/output_scalar.nc", recursive=True)

    # Find the output_scalar.nc files where run.done also exists
    valid_outputs = [file for file in output_files if os.path.dirname(file) in run_done_dirs]

    # Write the list to a comma-separated text file
    output_file_path = f"{level}_valid_outputs.txt"
    with open(output_file_path, "w") as f:
        f.write(",".join(valid_outputs))

    print(f"List of valid output files written to {output_file_path}")
    print(f"Number of run.done directories: {len(run_done_dirs)}")
    print(f"Number of output files: {len(output_files)}")
    print(f"Number of valid output files: {len(valid_outputs)}")

if __name__ == "__main__":
    try:
        arg_parser = argparse.ArgumentParser()
        arg_parser.add_argument("basin", help="Name of the basin")
        arg_parser.add_argument("level", help="Calibration level")
        args = arg_parser.parse_args()
        basin = args.basin
        level = args.level
        main(basin, level)
    except Exception as e:
        print(f"Error: {e}")
        traceback.print_exc()
        input("Press Enter to exit")