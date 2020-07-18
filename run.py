import sys
import time
import pkg_resources
from multiprocessing import Pool
from pathlib import Path
import argparse
import csv

from rd_filters.rd_filters import read_rules, RDFilters
import pandas as pd


def filter(rf: RDFilters, rule_dict:dict, input_path: str, output_path: str):
    file_name = Path(input_path).stem
    print(f"Reading {input_path}")
    orig_df = pd.read_feather(input_path)
    input_data = list(orig_df[['SMILE','ID']].to_records(index=False))
 
    start_time = time.time()
    pool = Pool(16)
    res = list(pool.map(rf.evaluate, input_data))
    df = pd.DataFrame(res, columns=["SMILES", "NAME", "FILTER", "MW", "LogP", "HBD", "HBA", "TPSA", "Rot"])
    df_ok = df[
        (df.FILTER == "OK") &
        df.MW.between(*rule_dict["MW"]) &
        df.LogP.between(*rule_dict["LogP"]) &
        df.HBD.between(*rule_dict["HBD"]) &
        df.HBA.between(*rule_dict["HBA"]) &
        df.TPSA.between(*rule_dict["TPSA"]) &
        df.Rot.between(*rule_dict["Rot"])
    ]
    # refine original dataframe
    df_final = orig_df[orig_df['ID'].isin(df_ok['NAME'])].reset_index(drop=True)
    save_path = str(Path(output_path, f'{file_name}.feather'))
    df_final.to_feather(save_path)

    # stat
    num_input_rows = df.shape[0]
    num_output_rows = df_ok.shape[0]
    fraction_passed = "%.1f" % (num_output_rows / num_input_rows * 100.0)
    print(f"{num_output_rows} of {num_input_rows} passed filters {fraction_passed}%", file=sys.stderr)
    elapsed_time = "%.2f" % (time.time() - start_time)
    print(f"Elapsed time {elapsed_time} seconds", file=sys.stderr)


def load_file_list(file_list_path:str):
    if not Path(file_list_path).exists():
        print(f"Cannot file {file_list_path}")
        sys.exit(1)

    with open(file_list_path) as f:
        file_list = []
        reader = csv.reader(f, delimiter=" ")
        for line in reader:
            file_list.append(line[0])

    return file_list


def parse_args():
    psr = argparse.ArgumentParser()
    psr.add_argument('--input', type=str, default=None, help='input file list')
    psr.add_argument('--out', default='./filtered/', help='output folder')
    args, _ = psr.parse_known_args()

    return args


def main():
    args = parse_args()

    # initialize filter class
    alert_file_name = pkg_resources.resource_filename('rd_filters', "data/alert_collection.csv")
    rd_filter = RDFilters(alert_file_name)
    rd_filter.build_rule_list(['PAINS'])
    rule_dict = {
        "MW": [0, 500],
        "LogP": [-5, 5],
        "HBD": [0, 5],
        "HBA": [0, 10],
        "TPSA": [0, 200],
        "Rot": [0, 10],
        "TPSA": [0, 200]
    }

    input_files = load_file_list(args.input)
    for input_file in input_files:
        filter(rd_filter, rule_dict, input_file, args.out)


if __name__ == '__main__':
    main()

