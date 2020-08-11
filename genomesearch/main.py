import click
from genomesearch import *
from genomesearch.help import CustomHelp
from genomesearch.isolate import _refbank as isolate_refbank
from genomesearch.isolate import _uhgg as isolate_uhgg
from genomesearch.meta import _refbank as meta_refbank
from genomesearch.meta import _uhgg as meta_uhgg

from genomesearch.install import _install

@click.group(cls=CustomHelp)
def cli():
    """A command line tool to quickly search for closely related microbial genomes using a marker-gene based approach."""
    pass

@cli.command(short_help='Download the GenomeSearch database')
@click.option('--threads', '-t', default=10)
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
def install(threads, force):
    log_params(threads=threads, force=force)
    _install(threads, force)


@cli.group(short_help='Run genomesearch on a complete or draft sequence of a single isolate.')
def isolate():
    """A click access point for the run module. This is used for creating the command line interface."""
    pass

@isolate.command(short_help='Run genomesearch isolate on the RefSeq/GenBank (RefBank) database.')
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--num-markers', '-m', default=40, help='The number of marker genes to use (default 40).')
@click.option('--outdir', '-o', default='genomesearch_isolate_refbank_output', help='The name of the output directory.')
@click.option('--prefix', '-prefix', default='genomesearch', help='The prefix of all files in the output directory.')
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
@click.option('--threads', '-t', default=16, help="Number of threads to use for diamond searches.")
@click.option('--max-target-seqs', '-k', default=200, help="The maximum number of target seqs returned by the diamond search.")
@click.option('--min-percent-identity', '-mpi', default=50, help="The minimum percent identity to keep a marker match from the diamond search.")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files.")
@click.option('--fasta-type', '-ft', type=click.Choice(['genome', 'proteome', 'markers']), default='genome', help="Select the type of fasta input.")
def refbank(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(fasta=fasta, num_markers=num_markers, outdir=outdir, prefix=prefix, force=force, threads=threads,
               max_target_seqs=max_target_seqs, min_percent_identity=min_percent_identity,
               keep_intermediate=keep_intermediate, fasta_type=fasta_type)
    isolate_refbank(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, min_percent_identity,
                    keep_intermediate, fasta_type)


@isolate.command(short_help='Run genomesearch isolate on the Unified Human Gastrointestinal Genome (UHGG) database.')
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--num-markers', '-m', default=40, help='The number of marker genes to use (default 40).')
@click.option('--outdir', '-o', default='genomesearch_isolate_uhgg_output', help='The name of the output directory.')
@click.option('--prefix', '-prefix', default='genomesearch', help='The prefix of all files in the output directory.')
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
@click.option('--threads', '-t', default=16, help="Number of threads to use for diamond searches.")
@click.option('--max-target-seqs', '-k', default=200, help="The maximum number of target seqs returned by the diamond search.")
@click.option('--min-percent-identity', '-mpi', default=50, help="The minimum percent identity to keep a marker match from the diamond search.")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files.")
@click.option('--fasta-type', '-ft', type=click.Choice(['genome', 'proteome', 'markers']), default='genome', help="Select the type of fasta input.")
def uhgg(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(fasta=fasta, num_markers=num_markers, outdir=outdir, prefix=prefix, force=force, threads=threads,
               max_target_seqs=max_target_seqs, min_percent_identity=min_percent_identity,
               keep_intermediate=keep_intermediate, fasta_type=fasta_type)
    isolate_uhgg(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type)


@cli.group(short_help='Run genomesearch on a metagenomic assembly with multiple unbinned sequences.')
def meta():
    """A click access point for the run module. This is used for creating the command line interface."""
    pass

@meta.command(short_help='Run genomesearch meta on the RefSeq/GenBank (RefBank) database.')
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--num-markers', '-m', default=400, help='The number of marker genes to use (default 40).')
@click.option('--outdir', '-o', default='genomesearch_meta_refbank_output', help='The name of the output directory.')
@click.option('--prefix', '-prefix', default='genomesearch', help='The prefix of all files in the output directory.')
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
@click.option('--threads', '-t', default=16, help="Number of threads to use for diamond searches.")
@click.option('--max-target-seqs', '-k', default=20, help="The maximum number of target seqs returned by the diamond search.")
@click.option('--min-percent-identity', '-mpi', default=80, help="The minimum percent identity to keep a marker match from the diamond search.")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files.")
@click.option('--fasta-type', '-ft', type=click.Choice(['genome', 'proteome', 'markers']), default='genome', help="Select the type of fasta input.")
def refbank(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(fasta=fasta, num_markers=num_markers, outdir=outdir, prefix=prefix, force=force, threads=threads,
               max_target_seqs=max_target_seqs, min_percent_identity=min_percent_identity, keep_intermediate=keep_intermediate, fasta_type=fasta_type)
    meta_refbank(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type)

@meta.command(short_help='Run genomesearch meta on the Unified Human Gastrointestinal Genome (UHGG) database.')
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--num-markers', '-m', default=400, help='The number of marker genes to use (default 40).')
@click.option('--outdir', '-o', default='genomesearch_meta_uhgg_output', help='The name of the output directory.')
@click.option('--prefix', '-prefix', default='genomesearch', help='The prefix of all files in the output directory.')
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
@click.option('--threads', '-t', default=16, help="Number of threads to use for diamond searches.")
@click.option('--max-target-seqs', '-k', default=20, help="The maximum number of target seqs returned by the diamond search.")
@click.option('--min-percent-identity', '-mpi', default=80, help="The minimum percent identity to keep a marker match from the diamond search.")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files.")
@click.option('--fasta-type', '-ft', type=click.Choice(['genome', 'proteome', 'markers']), default='genome', help="Select the type of fasta input.")
def uhgg(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(fasta=fasta, num_markers=num_markers, outdir=outdir, prefix=prefix, force=force, threads=threads,
               max_target_seqs=max_target_seqs, min_percent_identity=min_percent_identity, keep_intermediate=keep_intermediate, fasta_type=fasta_type)
    meta_uhgg(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type)

def log_params(**kwargs):
    click.echo("#### PARAMETERS ####")
    click.echo('\n'.join(list(map(lambda x: ': '.join(list(map(str, x))), kwargs.items()))))
    click.echo("####################")

if __name__ == '__main__':

    cli()