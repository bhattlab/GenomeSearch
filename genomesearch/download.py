from genomesearch import *
import click
import wget
from os import makedirs

def _download():
    click.echo("#### INPUT PARAMETERS ####")
    try:
        num_markers = int(input("How many markers do you want to use? (This can be any number between 1 and 400)\n[default=150]>> ") or "150")
        if num_markers > 400 or num_markers < 1:
            raise Exception('wrong_markers')
    except:
        print("ERROR!")
        print("Please input a number between 1 and 400, the default is 150.")
        num_markers = int(input("How many markers do you want to use?]\n[default=150] >> ") or "150")
    click.echo("####################")

    if not isfile(SQLDB_PATH):
        print("Downloading DSN1 model...")
        makedirs(dirname(SQLDB_PATH), exist_ok=True)
        wget.download('https://storage.googleapis.com/genomesearch/downloads/genomesearch.db', SQLDB_PATH)
        print()