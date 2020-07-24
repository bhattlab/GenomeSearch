import click
from genomesearch import *
from genomesearch.help import CustomHelp
from genomesearch.search import _refbank, _meta
from genomesearch.download import _download

@click.group(cls=CustomHelp)
def cli():
    """A command line tool to quickly search for closely related microbial genomes using a marker-gene based approach."""
    pass

@cli.command(short_help='Download the GenomeSearch database', help_priority=1)
@click.option('--threads', '-t', default=10)
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
def download(threads, force):
    log_params(threads=threads, force=force)
    _download(threads, force)


@cli.command(short_help='Run genomesearch on a complete or draft sequence of a single species against refseq/genbank genomes.', help_priority=1)
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--num-markers', '-m', default=40, help='The number of marker genes to use (default 40).')
@click.option('--outdir', '-o', default='genomesearch_output', help='The name of the output directory.')
@click.option('--prefix', '-prefix', default='genomesearch', help='The prefix of all files in the output directory.')
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
@click.option('--threads', '-t', default=16, help="Number of threads to use for diamond searches.")
@click.option('--max-target-seqs', '-k', default=200, help="The maximum number of target seqs returned by the diamond search.")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files.")
@click.option('--fasta-type', '-ft', type=click.Choice(['genome', 'proteome', 'markers']), default='genome', help="Select the type of fasta input.")
def refbank(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, keep_intermediate, fasta_type):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(fasta=fasta, num_markers=num_markers, outdir=outdir, prefix=prefix, force=force, threads=threads,
               max_target_seqs=max_target_seqs, keep_intermediate=keep_intermediate, fasta_type=fasta_type)
    _refbank(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, keep_intermediate, fasta_type=fasta_type)


@cli.command(short_help='Run genomesearch on a complete or draft sequence of a single species against human gut metagenome-assembled genomes (MAGs).', help_priority=1)
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--num-markers', '-m', default=40, help='The number of marker genes to use (default 40).')
@click.option('--outdir', '-o', default='genomesearch_output', help='The name of the output directory.')
@click.option('--prefix', '-prefix', default='genomesearch', help='The prefix of all files in the output directory.')
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
@click.option('--threads', '-t', default=16, help="Number of threads to use for diamond searches.")
@click.option('--max-target-seqs', '-k', default=200, help="The maximum number of target seqs returned by the diamond search.")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files.")
@click.option('--fasta-type', '-ft', type=click.Choice(['genome', 'proteome', 'markers']), default='genome', help="Select the type of fasta input.")
def meta(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, keep_intermediate, fasta_type):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(fasta=fasta, num_markers=num_markers, outdir=outdir, prefix=prefix, force=force, threads=threads,
               max_target_seqs=max_target_seqs, keep_intermediate=keep_intermediate, fasta_type=fasta_type)
    _meta(fasta, num_markers, outdir, prefix, force, threads, max_target_seqs, keep_intermediate, fasta_type=fasta_type)


def log_params(**kwargs):
    click.echo("#### PARAMETERS ####")
    click.echo('\n'.join(list(map(lambda x: ': '.join(list(map(str, x))), kwargs.items()))))
    click.echo("####################")

if __name__ == '__main__':

    cli()