from os import makedirs
from os.path import isdir, join
import shutil
import click
import sys
import os
from genomesearch.prodigal import run_prodigal
from genomesearch import *
from Bio import SeqIO
from subprocess import run

def _search(fasta, num_markers, outdir, prefix, force):

    tmpdir = join(outdir, 'tmp')
    if force and isdir(outdir):
        shutil.rmtree(outdir)
    try:
        makedirs(tmpdir)
    except FileExistsError:
        click.echo("Output directory exists, please delete or overwrite with --force")
        sys.exit(1)

    click.echo("Running prodigal...")
    run_prodigal(PRODIGAL_PATH, fasta, tmpdir, meta=False)


def get_marker_genes(protein_fasta_path):
    run('diamond blastp --query {0} --out {1}.dmd.tsv --outfmt 6 --db {2}'.format(query, outfile, database).split())

    top_markers = dict()
    with open(outfile + '.dmd.tsv') as infile:
        for line in infile:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore =  line.strip().split('\t')
            length, evalue, bitscore = int(length), float(evalue), float(bitscore)
            finding = (qseqid, length, evalue, bitscore)
            marker = sseqid.split('_')[1]

            if finding[-2] >= 1e-4:
                continue

            if marker not in top_markers:
                top_markers[marker] = finding
            else:
                if finding[-1] > top_markers[marker][-1]:
                    top_markers[marker] = finding

    marker2gene = dict()
    gene2marker = dict()
    for rec in top_markers:
        marker2gene[rec] = top_markers[rec][0]
        gene2marker[top_markers[rec][0]] = rec


    records = []
    for rec in SeqIO.parse(query, 'fasta'):
        if rec.id in gene2marker:
            rec.id = gene2marker[rec.id] + '__' + rec.id + '__' + record_suffix
            rec.description = rec.id
            records.append(rec)

    SeqIO.write(records, outfile, 'fasta')
    os.remove(outfile+'.dmd.tsv')