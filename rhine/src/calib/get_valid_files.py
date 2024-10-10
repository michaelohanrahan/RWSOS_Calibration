from glob import glob
import os
import sys 
import argparse
import traceback

def main(basin, level):
    if sys.platform == "win32":
        DRIVE = "p:"
    elif sys.platform == "linux":
        DRIVE = "/p"
    else:
        raise ValueError(f"Unsupported platform: {sys.platform}")
    
    os.chdir(f"{DRIVE}/11209265-grade2023/wflow/RWSOS_Calibration/{basin}")
    # Get the list of directories containing run.done
    run_done_files = glob(f"data/2-interim/calib_data/{level}/**/run.done", recursive=True)
    run_done_dirs = {os.path.dirname(file) for file in run_done_files}
    # Get the list of output_scalar.nc files
    output_files = glob(f"data/2-interim/calib_data/{level}/**/output_scalar.nc", recursive=True)

    # Find the output_scalar.nc files where run.done also exists
    valid_outputs = [file for file in output_files if os.path.dirname(file) in run_done_dirs]
    non_valid_outputs = [file for file in output_files if os.path.dirname(file) not in run_done_dirs]
    
    print(f"Number of valid output files: {len(valid_outputs)}")
    first_file_size = None
    if valid_outputs:
        first_file_size = os.path.getsize(valid_outputs[0])
        print(f"First file size: {first_file_size}")
        valid_outputs = [file for file in valid_outputs if os.path.getsize(file) == first_file_size]
    print(f"Number of valid output files after filtering by size: {len(valid_outputs)}")
    
    print(f"Number of non-valid output files: {len(non_valid_outputs)}")
    
    invalid_with_size = []
    if non_valid_outputs and first_file_size is not None:
        invalid_with_size = [file for file in non_valid_outputs if os.path.getsize(file) == first_file_size]
    
    if invalid_with_size:
        print(f"Number of invalid output files with the same size as the first valid file: {len(invalid_with_size)}")
        valid_outputs.extend(invalid_with_size)
        non_valid_outputs = [file for file in non_valid_outputs if file not in invalid_with_size]
    else:
        print("No invalid output files with the same size as the first valid file found")
    
    print(f"Newly extended valid outputs: {len(valid_outputs)}")
    
    count = 0
    for file in valid_outputs:
        done = os.path.join(os.path.dirname(file), "run.done")
        if not os.path.exists(done):
            with open(done, "w") as f:
                f.write("")
                count += 1
    print(f"Number of run.done files created: {count}")
    
    #FULL PATH
    valid_outputs = [f"{os.getcwd()}/{file}" for file in valid_outputs]
    
    # Write the list to a comma-separated text file
    output_file_path = f"{level}_valid_outputs.txt"
    with open(output_file_path, "w") as f:
        f.write(",".join(valid_outputs))

    non_valid_outputs_file_path = f"{level}_non_valid_outputs.txt"
    with open(non_valid_outputs_file_path, "w") as f:
        f.write(",".join(non_valid_outputs))
    
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