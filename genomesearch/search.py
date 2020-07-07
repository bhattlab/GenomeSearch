from os import makedirs
from os.path import isdir, join
import shutil
import click
import sys
from genomesearch.prodigal import run_prodigal
from genomesearch import *

def _search(fasta, num_markers, outdir, prefix, force):

    tmpdir = join(outdir, 'tmp')
    if force and isdir(outdir):
        shutil.rmtree(outdir)
    try:
        makedirs(tmpdir)
    except FileExistsError:
        click.echo("Output directory exists, please delete or overwrite with --force")
        sys.exit(1)

    run_prodigal(PRODIGAL_PATH, fasta, tmpdir, meta=False)