#!/usr/bin/env python3

# -----------------------------------------------------

# Download FASTQ files from ENA using aria2c.
# This script is a modified version of fastq-dl, a tool
# to download FASTQ files from SRA or ENA. The original
# script can be found here: https://github.com/rpetit3/fastq-dl.
#
# The purpose of this script is to borrow some of the
# functionality of fastq-dl to download FASTQ files but
# make it portable, lightweight, and easy to run.

# Usage:
#
# ./get.py -a <accession> -o <outdir> -m <max_attempts> \\
# -s <sleep> -f <force> -g <group_by_experiment> -G <group_by_sample> -p <prefix>
#
# The base use case could be achieve just by:
#
# ./get.py -a <accession> -o <outdir>
#
# -----------------------------------------------------

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "alejandrogzi"
__version__ = "0.0.1"
__credits__ = ["Robert A. Petit III"]

# -----------------------------------------------------

import time
import sys
import re
import requests
import hashlib
from pathlib import Path
from executor import ExternalCommand, ExternalCommandFailed
import csv
import argparse


ENA_FAILED = "ENA_NOT_FOUND"
ENA = "ENA"
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv"


def get_run_info(
    query: str,
    max_attempts: int = 10,
    sleep: int = 10,
) -> tuple:
    """Retrieve a list of samples available from ENA.

    Args:
        query (str): A formatted query for ENA searches.
        max_attempts (int, optional): Maximum number of download attempts
        sleep (int): Minimum amount of time to sleep before retry

    Returns:
        tuple: Records associated with the accession.
    """
    attempt = 1
    ena_attempt = 1
    while True:
        success, ena_data = get_ena_metadata(query)
        if success:
            return ENA, ena_data
        elif attempt >= max_attempts:
            print("There was an issue querying ENA, exiting...")
            print(f"STATUS: {ena_data[0]}")
            print(f"TEXT: {ena_data[1]}")
            sys.exit(1)
        time.sleep(sleep)


def validate_query(query: str) -> str:
    """
    Check that query is an accepted accession type and return the accession type. Current
    accepted types are:

        Projects - PRJEB, PRJNA, PRJDA
        Studies - ERP, DRP, SRP
        BioSamples - SAMD, SAME, SAMN
        Samples - ERS, DRS, SRS
        Experiments - ERX, DRX, SRX
        Runs - ERR, DRR, SRR

    Parameters:
        query (str): A string containing an accession.

    Returns:
        str: A string containing the query for ENA search.

    https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html
    """
    if re.match(r"^PRJ[EDN][A-Z][0-9]+$|^[EDS]RP[0-9]{6,}$", query):
        # Is a project or study accession
        return f"(study_accession={query} OR secondary_study_accession={query})"
    elif re.match(r"^SAM[EDN][A-Z]?[0-9]+$|^[EDS]RS[0-9]{6,}$", query):
        # Is a sample or biosample accession
        return f"(sample_accession={query} OR secondary_sample_accession={query})"
    elif re.match(r"^[EDS]RX[0-9]{6,}$", query):
        # Is an experiment accession
        return f"experiment_accession={query}"
    elif re.match(r"^[EDS]RR[0-9]{6,}$", query):
        # Is a run accession
        return f"run_accession={query}"
    else:
        print(
            f"{query} is not a Study, Sample, Experiment, or Run accession. See https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html for valid options"
        )
        sys.exit(1)


def get_ena_metadata(query: str) -> list:
    """Fetch metadata from ENA.
    https://docs.google.com/document/d/1CwoY84MuZ3SdKYocqssumghBF88PWxUZ/edit#heading=h.ag0eqy2wfin5

    Args:
        query (str): The query to search for.

    Returns:
        list: Records associated with the accession.
    """
    url = f'{ENA_URL}&query="{query}"&fields=all'
    headers = {"Content-type": "application/x-www-form-urlencoded"}
    r = requests.get(url, headers=headers)
    if r.status_code == requests.codes.ok:
        data = []
        col_names = None
        for line in r.text.split("\n"):
            cols = line.rstrip().split("\t")
            if line:
                if col_names:
                    data.append(dict(zip(col_names, cols)))
                else:
                    col_names = cols
        if data:
            return [True, data]
        else:
            return [
                False,
                [r.status_code, "Query was successful, but received an empty response"],
            ]
    else:
        return [False, [r.status_code, r.text]]


def ena_download(
    run: str, outdir: str, max_attempts: int = 10, force: bool = False, sleep: int = 10
) -> dict:
    """Download FASTQs from ENA FTP using wget.

    Args:
        accession (str): The accession to download associated FASTQs.
        outdir (str): Directory to write FASTQs to.
        max_attempts (int, optional): Maximum number of download attempts. Defaults to 10.
        force: (bool, optional): Whether to overwrite existing files if the MD5's do not match
        sleep (int): Minimum amount of time to sleep before retry

    Returns:
        dict: A dictionary of the FASTQs and their paired status.
    """
    fastqs = {"r1": "", "r2": "", "single_end": True}
    ftp = run["fastq_ftp"]
    if not ftp:
        return ENA_FAILED

    ftp = ftp.split(";")
    md5 = run["fastq_md5"].split(";")
    for i in range(len(ftp)):
        is_r2 = False
        # If run is paired only include *_1.fastq and *_2.fastq, rarely a
        # run can have 3 files.
        # Example:ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR114/007/ERR1143237
        if run["library_layout"] == "PAIRED":
            if ftp[i].endswith("_2.fastq.gz"):
                # Example: ERR1143237_2.fastq.gz
                is_r2 = True
            elif ftp[i].endswith("_1.fastq.gz"):
                # Example: ERR1143237_1.fastq.gz
                pass
            else:
                # Example: ERR1143237.fastq.gz
                # Not a part of the paired end read, so skip this file. Or,
                # its the only fastq file, and its not a paired
                obs_fq = Path(ftp[i]).name
                exp_fq = f'{run["run_accession"]}.fastq.gz'
                if len(ftp) != 1 and obs_fq != exp_fq:
                    continue

        # Download Run
        if md5[i]:
            fastq = download_ena_fastq(
                ftp[i],
                outdir,
                md5[i],
                max_attempts=max_attempts,
                force=force,
                sleep=sleep,
            )
            if fastq == ENA_FAILED:
                return ENA_FAILED

            if is_r2:
                fastqs["r2"] = fastq
                fastqs["single_end"] = False
            else:
                fastqs["r1"] = fastq

    return fastqs


def download_ena_fastq(
    ftp: str,
    outdir: str,
    md5: str,
    max_attempts: int = 10,
    force: bool = False,
    sleep: int = 10,
) -> str:
    """Download FASTQs from ENA using FTP.

    Args:
        ftp (str): The FTP address of the FASTQ file.
        outdir (str): Directory to download the FASTQ to.
        md5 (str): Expected MD5 checksum of the FASTQ.
        max_attempts (int, optional): Maximum number of download attempts. Defaults to 10.
        force: (bool, optional): Whether to overwrite existing files if the MD5's do not match
        sleep (int): Minimum amount of time to sleep before retry

    Returns:
        str: Path to the downloaded FASTQ.
    """
    success = False
    attempt = 0
    fastq = f"{outdir}/{Path(ftp).name}"

    if Path(fastq).exists() and force:
        print(f"Overwriting existing file: {fastq}")
        Path(fastq).unlink()

    if not Path(fastq).exists():
        Path(outdir).mkdir(parents=True, exist_ok=True)

        while not success:
            print(f"\t\t{Path(ftp).name} FTP download attempt {attempt + 1}")
            outcome = execute(
                f"aria2c -c -x 4 -o {fastq} http://{ftp}",
                max_attempts=max_attempts,
                sleep=sleep,
            )
            if outcome == ENA_FAILED:
                return ENA_FAILED

            if force:
                print(f"--force used, skipping MD5sum check for {fastq}")
                success = True
            else:
                fastq_md5 = md5sum(fastq)
                if fastq_md5 != md5:
                    print(f"MD5s, Observed: {fastq_md5}, Expected: {md5}")
                    attempt += 1
                    if Path(fastq).exists():
                        Path(fastq).unlink()
                    if attempt > max_attempts:
                        print(
                            f"Download failed after {max_attempts} attempts. "
                            "Please try again later or manually from SRA/ENA."
                        )
                        sys.exit(1)
                else:
                    success = True
    else:
        print(f"Skipping re-download of existing file: {fastq}")

    return fastq


def execute(
    cmd: str,
    directory: str = str(Path.cwd()),
    capture_stdout: bool = False,
    stdout_file: str = None,
    stderr_file: str = None,
    max_attempts: int = 1,
    sleep: int = 10,
) -> str:
    """A simple wrapper around executor.

    Args:
        cmd (str): A command to execute.
        directory (str, optional): Set the working directory for command. Defaults to str(Path.cwd()).
        capture_stdout (bool, optional): Capture and return the STDOUT of a command. Defaults to False.
        stdout_file (str, optional): File to write STDOUT to. Defaults to None.
        stderr_file (str, optional): File to write STDERR to. Defaults to None.
        max_attempts (int, optional): Maximum times to attempt command execution. Defaults to 1.
        is_sra (bool, optional): The command is from SRA. Defaults to False.
        sleep (int): Minimum amount of time to sleep before retry

    Raises:
        error: An unexpected error occurred.

    Returns:
        str: Exit code, accepted error message, or STDOUT of command.
    """
    attempt = 0
    while attempt < max_attempts:
        attempt += 1
        try:
            command = ExternalCommand(
                cmd,
                directory=directory,
                capture=True,
                capture_stderr=True,
                stdout_file=stdout_file,
                stderr_file=stderr_file,
            )

            command.start()
            # print(command.decoded_stdout)
            # print(command.decoded_stderr)

            if capture_stdout:
                return command.decoded_stdout
            else:
                return command.returncode
        except ExternalCommandFailed:
            print(f'"{cmd}" return exit code {command.returncode}')

            if attempt < max_attempts:
                print(f"Retry execution ({attempt} of {max_attempts})")
                time.sleep(sleep)
            else:
                return ENA_FAILED


def merge_runs(runs: list, output: str) -> None:
    """Merge runs from an experiment or sample.

    Args:
        runs (list): A list of FASTQs to merge.
        output (str): The final merged FASTQ.
    """
    if len(runs) > 1:
        run_fqs = " ".join(runs)
        execute(f"cat {run_fqs} > {output}")
        for p in runs:
            Path(p).unlink()
    else:
        Path(runs[0]).rename(output)


def write_tsv(data: dict, output: str) -> None:
    """Write a TSV file.

    Args:
        data (dict): Data to be written to TSV.
        output (str): File to write the TSV to.
    """
    with open(output, "w") as fh:
        if output.endswith("-run-mergers.tsv"):
            writer = csv.DictWriter(
                fh, fieldnames=["accession", "r1", "r2"], delimiter="\t"
            )
            writer.writeheader()
            for accession, vals in data.items():
                writer.writerow(
                    {
                        "accession": accession,
                        "r1": ";".join(vals["r1"]),
                        "r2": ";".join(vals["r2"]),
                    }
                )
        else:
            writer = csv.DictWriter(fh, fieldnames=data[0].keys(), delimiter="\t")
            writer.writeheader()
            for row in data:
                writer.writerow(row)


def md5sum(fastq: str) -> str:
    """Calculate the MD5 checksum of a file.

    Source: https://stackoverflow.com/a/3431838/5299417

    Args:
        fastq (str): Input FASTQ to calculate MD5 checksum for.

    Returns:
        str: Calculated MD5 checksum.
    """
    MB = 1_048_576
    BUFFER_SIZE = 10 * MB
    if Path(fastq).exists():
        hash_md5 = hashlib.md5()
        with open(fastq, "rb") as fp:
            for chunk in iter(lambda: fp.read(BUFFER_SIZE), b""):
                hash_md5.update(chunk)

        return hash_md5.hexdigest()
    else:
        return None


def get_fastqs(
    accession,
    group_by_experiment,
    group_by_sample,
    outdir,
    max_attempts,
    sleep,
    force,
    prefix,
):
    """Download FASTQ files from ENA or SRA."""
    query = validate_query(accession)
    data_from, ena_data = get_run_info(
        query,
        max_attempts=max_attempts,
        sleep=sleep,
    )

    print(f"Total Runs Found: {len(ena_data)}")
    downloaded = {}
    runs = {} if group_by_experiment or group_by_sample else None
    outdir = "DOWNLOAD" if outdir == "./" else f"{outdir}"

    for i, run_info in enumerate(ena_data):
        run_acc = run_info["run_accession"]
        if run_acc not in downloaded:
            downloaded[run_acc] = True
        else:
            print(f"Duplicate run {run_acc} found, skipping re-download...")
            continue
        print(f"\tWorking on run {run_acc}...")
        fastqs = None

        fastqs = ena_download(
            run_info,
            outdir,
            max_attempts=max_attempts,
            force=force,
            sleep=sleep,
        )

        if fastqs == ENA_FAILED:
            print(f"\tNo fastqs found in ENA for {run_acc}")
            ena_data[i]["error"] = ENA_FAILED
            fastqs = None

        if fastqs:
            if group_by_experiment or group_by_sample:
                name = run_info["sample_accession"]
                if group_by_experiment:
                    name = run_info["experiment_accession"]

                if name not in runs:
                    runs[name] = {"r1": [], "r2": []}

                if fastqs["single_end"]:
                    runs[name]["r1"].append(fastqs["r1"])
                else:
                    runs[name]["r1"].append(fastqs["r1"])
                    runs[name]["r2"].append(fastqs["r2"])

    # If applicable, merge runs
    if runs:
        for name, vals in runs.items():
            if len(vals["r1"]) and len(vals["r2"]):
                # Not all runs labeled as paired are actually paired.
                if len(vals["r1"]) == len(vals["r2"]):
                    print(f"\tMerging paired end runs to {name}...")
                    merge_runs(vals["r1"], f"{outdir}/{name}_R1.fastq.gz")
                    merge_runs(vals["r2"], f"{outdir}/{name}_R2.fastq.gz")
                else:
                    print("\tMerging single end runs to experiment...")
                    merge_runs(vals["r1"], f"{outdir}/{name}.fastq.gz")
            else:
                print("\tMerging single end runs to experiment...")
                merge_runs(vals["r1"], f"{outdir}/{name}.fastq.gz")
        print(f"Writing merged run info to {outdir}/{prefix}-run-mergers.tsv")
        write_tsv(runs, f"{outdir}/{prefix}-run-mergers.tsv")
    print(f"Writing metadata to {outdir}/{prefix}-run-info.tsv")
    write_tsv(ena_data, f"{outdir}/{prefix}-run-info.tsv")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Download FASTQ files from ENA or SRA."
    )
    parser.add_argument(
        "-a",
        "--accession",
        type=str,
        required=True,
        help="A valid ENA or SRA accession.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        default="./",
        help="Directory to write FASTQs to.",
    )
    parser.add_argument(
        "-m",
        "--max-attempts",
        type=int,
        default=10,
        help="Maximum number of download attempts.",
    )
    parser.add_argument(
        "-s",
        "--sleep",
        type=int,
        default=10,
        help="Minimum amount of time to sleep before retry.",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force overwrite of existing files.",
    )
    parser.add_argument(
        "-g",
        "--group-by-experiment",
        action="store_true",
        help="Group FASTQs by experiment.",
    )
    parser.add_argument(
        "-G",
        "--group-by-sample",
        action="store_true",
        help="Group FASTQs by sample.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default="fastq",
        help="Prefix for output files.",
    )
    args = parser.parse_args()

    return args


def main():
    """Main function."""
    args = parse_args()

    if args.accession.endswith(".txt"):
        with open(args.accession, "r") as fh:
            for line in fh:
                accession = line.strip()
                if accession:
                    print(f"Working on accession {accession}...")
                    get_fastqs(
                        accession,
                        args.group_by_experiment,
                        args.group_by_sample,
                        args.outdir,
                        args.max_attempts,
                        args.sleep,
                        args.force,
                        args.prefix,
                    )
    else:
        get_fastqs(
            args.accession,
            args.group_by_experiment,
            args.group_by_sample,
            args.outdir,
            args.max_attempts,
            args.sleep,
            args.force,
            args.prefix,
        )

    return None


if __name__ == "__main__":
    main()
