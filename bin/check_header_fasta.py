#!/usr/bin/env python

import argparse
import logging
import sys
from pathlib import Path
from Bio import SeqIO

logger = logging.getLogger()

def check_fasta(file_in, file_out):
	# print('\033[94m The header of assembly {} is not compatible with REAPR.\nPlease, execute reapr facheck and modify the file assemblies.tsv to point to assembly with corrected header! \033[0m \n\n\n To this y using command line: \033[96m reapr facheck <in.fa> [out_prefix] \033[0m'.format(input.genome))
	sequences = []
	with open(file_in, "r") as ingenome:
			for record in SeqIO.parse(ingenome, "fasta"):
					header = str(record.description)
					# Things like trailing whitespace or characters |':- could break the pipeline.
					replace_dict = {"|": "_", "'": "_", ":": "_", "-": "_", " ": "_"}
					record_id = record.id
					description = record.description
					for old, new in replace_dict.items():
						record_id = record_id.replace(old, new)
					record.id = record_id
					record.description=""
					sequences.append(record)
	destination = Path(file_out)
	# destination = Path(input.genome+".bkp")
	# with destination.open(mode="xb") as fid:
	# 		fid.write(source.read_bytes())
	with open(str(destination), "w") as outgenome:
			SeqIO.write(sequences, outgenome, "fasta")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Transform record fasta header if it contains things like trailing whitespace or characters |':- that could break the REAPR execution.",
        epilog="Example: python check_header_fasta.py sample.fasta sample_valid.fasta",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Fasta input.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Fasta output.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_fasta(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())