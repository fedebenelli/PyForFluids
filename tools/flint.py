# This tool only exits because fontran-linter dont support directories as
# inputs

# This tools is part of PyForFluids
# https://github.com/fedebenelli/PyForFluids

# =============================================================================
# IMPORTS
# =============================================================================

import argparse
import itertools as it
import os
import pathlib
import sys

from fortran_linter import cli

# =============================================================================
# CONSTANTS
# =============================================================================
GLOBS = ["*.f90", "*.f95"]

# =============================================================================
# FUNCTIONS
# =============================================================================


def _expand_files(file_or_dir):
    files = []
    if os.path.isdir(file_or_dir):
        path = pathlib.Path(file_or_dir)
        for glob in GLOBS:
            files.extend([str(p) for p in path.glob(glob)])
    else:
        files.append(file_or_dir)  # always return a collection
    return files


def main(args=None):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "input",
        nargs="+",
        type=_expand_files,
        help=(
            "Input file(s) or directories.\n"
            "If the input is a directory all files with extension "
            "f90 and f95 are checked."
        ),
    )

    args = parser.parse_args(args)

    files = list(it.chain(*args.input))
    return cli.main(["--syntax-only"] + files)


if __name__ == "__main__":
    main(sys.argv[1:])
