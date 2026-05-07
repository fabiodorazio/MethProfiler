# import custom packages
import utils as utl
import functions as fct
import plots

# import packages
import argparse
import logging

log = logging.getLogger()

if __name__ == "__main__":
    log = logging.getLogger("main")
    usage = ""

    parser = argparse.ArgumentParser(description='Some description')

    parser.add_argument(
        "--config",
        type=str,
        required=True,
        dest="Config", # name of the attribute to be added to arg -> check args.Config
        help="some text help around config file",
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        dest="Output_dir",
        help="some text help around output dir",
    )

    parser.set_defaults(isForceOut=False)
    args = parser.parse_args()
    print(args)

    # use arguments: 
    # read_input being a function in utils.py
    config_file = read_input(args.Config)

    for index, row in config_file.iterrows():
        if not row["path_to_file"].endswith(
            ".gz"
        ):  # os.path.exists(row['path_to_vcf']+'.tbi')
            errMsg = "Please gzip the file, see README file for more info"
            sys.stderr.write(errMsg + "\nProgram exits.")
            sys.exit(1)
    
    # check_output_dir and create_log_dir being functions in utils.py
    output_dir = check_output_dir(args.Output_dir)
    create_log_dir(output_dir)
