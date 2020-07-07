import click
from genomesearch import *
from genomesearch.help import CustomHelp
from genomesearch.search import _search
from genomesearch.download import _download

@click.group(cls=CustomHelp)
def cli():
    """A command line tool to quickly search for closely related microbial genomes using a marker-gene based approach."""
    pass

@cli.command(short_help='Download the GenomeSearch database', help_priority=1)
def download():
    _download()



@cli.command(short_help='Run deepsmorfnet on a complete or draft sequence of a single species.', help_priority=1)
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--outfile', '-o', default='genomesearch_output.tsv')
def search(fasta, outfile):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(fasta=fasta, outfile=outfile)

    _search(fasta, outfile)



def log_params(**kwargs):
    click.echo("#### PARAMETERS ####")
    click.echo('\n'.join(list(map(lambda x: ': '.join(list(map(str, x))), kwargs.items()))))
    click.echo("####################")

if __name__ == '__main__':

    cli()