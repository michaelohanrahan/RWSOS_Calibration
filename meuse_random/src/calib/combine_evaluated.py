from pathlib import Path
from setuplog import setup_logging
import traceback
import pandas as pd

def read_results_to_dataframe(results_file: Path) -> pd.DataFrame:
    """
    Read the results txt file into a DataFrame.

    Args:
    results_file (Path): Path to the results txt file.

    Returns:
    pd.DataFrame: DataFrame with gauge measurements as columns and parameter sets as index.
    """
    with open(results_file, 'r') as f:
        header = f.readline().strip().split(',')
        #first element is string the rest are float
        header = [header[0]] + [float(h) for h in header[1:]]
        data = []
        for line in f:
            # Split the line at the closing curly brace
            parts = line.strip().split('}', 1)
            if len(parts) == 2:
                param_set = parts[0] + '}'  # Add the closing brace back
                values = parts[1].strip().split(',')
                data.append((param_set, values))

    # Extract parameter sets and values
    param_sets = [row[0] for row in data]
    values = []
    for row in data:
        row_values = []
        for val in row[1]:
            try:
                row_values.append(float(val))
            except ValueError:
                row_values.append(float('nan'))  # Convert empty strings or invalid floats to NaN
        #remove the nan values
        row_values = [val for val in row_values if not pd.isna(val)]
        values.append(row_values)

    # Create DataFrame
    df = pd.DataFrame(values, columns=header[1:], index=param_sets)
    df.index.name = 'params'

    return df

def main(
    l,
    results_file: str|Path,
    out: str|Path,
    best_n: int = 10,
    level: str = None,
):
    df = read_results_to_dataframe(results_file)
    # ic(df)

    best_params = []
    for gauge in df.columns:
        #gauge is a float column
        valid_euclidean = df[gauge]
        valid_ranked = valid_euclidean.sort_values(ascending=True)
        # ic(valid_ranked)
        top_n = min(best_n, len(valid_euclidean))
        best_for_gauge = {'gauge': gauge}
        for idx in range(top_n):
            param_dict =valid_ranked.index[idx]
            best_for_gauge[f'Top_{idx+1}'] = param_dict
        best_params.append(best_for_gauge)


    #TODO: Add the level column to best_params csv
    # Convert the results to a DataFrame and save as CSV
    level_int = int(level.split("level")[-1])
    if level_int == 0:
        final_df = pd.DataFrame(best_params)
        final_df['level'] = level_int
    else:
        inter_df = pd.DataFrame(best_params)
        inter_df['level'] = level_int
        prev_df = pd.read_csv(Path(out).parent.parent / f"level{level_int-1}" / "best_params.csv")
        # append the new rows, there will be an extension of gauges and level,
        # but the columns will be the same Top_1 ... Top_n
        final_df = pd.concat([prev_df, inter_df], ignore_index=True)
        # this way the next phase has access to the best 10 of all levels

    # ic(final_df)
    final_df.to_csv(out, index=False)

    # Create an 'eval.done' file to indicate completion
    with open(Path(out).parent / "level.done", "w") as f:
        f.write("")
    l.info(f"Saved {out}")


if __name__ == "__main__":
    l = setup_logging("data/0-log", "05-combine_evaluated.log")
    try:
        if 'snakemake' in globals():
            main(
                l=l,
                results_file=snakemake.input.results_file,
                out=snakemake.output.best_params,
                level=snakemake.params.level,
            )
        else:
            base_dir = r"/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level0"
            l.info(f"Base dir exists: {Path(base_dir).exists()}\n{base_dir}")

            results_file = Path(base_dir) / "results_level0.txt"
            main(
                l=l,
                results_file=results_file,
                out="/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level0/best_params.csv",
                level="level0",
            )
    except Exception as e:
        l.exception(e)
        l.error(traceback.format_exc())
        raise