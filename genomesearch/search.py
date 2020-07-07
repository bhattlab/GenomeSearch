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
from collections import defaultdict
import sqlite3
from glob import glob
import numpy as np

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

    marker_output = join(outdir, prefix+'.markers.faa')
    click.echo("Identifying marker genes...")
    get_marker_genes(join(tmpdir, 'prodigal.faa'), marker_output, prefix)

    click.echo("Searching for closest genomes in database...")
    get_closest_genomes(marker_output, num_markers, tmpdir)


def get_marker_genes(protein_fasta_path, outfile, prefix):
    run('diamond blastp --query {0} --out {1}.dmd.tsv --outfmt 6 --db {2}'.format(protein_fasta_path, outfile, PHYLOPHLAN_MARKER_PATH).split())

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
    for rec in SeqIO.parse(protein_fasta_path, 'fasta'):
        if rec.id in gene2marker:
            rec.id = gene2marker[rec.id] + '__' + rec.id + '__' + prefix
            rec.description = rec.id
            records.append(rec)

    SeqIO.write(records, outfile, 'fasta')
    os.remove(outfile+'.dmd.tsv')

def get_closest_genomes(marker_genes_fasta, num_markers, outdir):

    markers = []
    with open(MARKER_RANKS_PATH) as infile:
        for line in infile:
            marker = line.strip()
            markers.append(marker)

    markers = markers[:num_markers]

    conn = sqlite3.connect(SQLDB_PATH)
    c = conn.cursor()

    c.execute("SELECT genome_id, taxon_id FROM genome;")

    genome2taxid = dict()
    for line in c.fetchall():
        genome2taxid[line[0]] = line[1]

    c.execute("SELECT taxon_id,phylum,species FROM taxon;")
    taxon2species = dict()
    for line in c.fetchall():
        taxon2species[line[0]] = (line[1], line[2])

    split_markers_dir = os.path.join(outdir, 'markers')
    diamond_dir = os.path.join(outdir, 'diamond')

    os.makedirs(split_markers_dir, exist_ok=True)
    os.makedirs(diamond_dir, exist_ok=True)

    for rec in SeqIO.parse(marker_genes_fasta, 'fasta'):
        marker = rec.id.split('__')[0]
        SeqIO.write([rec], os.path.join(split_markers_dir, marker + '.faa'), 'fasta')

    for marker in markers:
        db = join(UNIQUE_MARKERS_PATH, marker + '.unique.dmnd')
        marker = os.path.basename(db).split('.')[0]

        run('diamond blastp -k 1000 --query {0} --out {1}.dmd.tsv --outfmt 6 --db {2}'.format(
            os.path.join(split_markers_dir, marker + '.faa'), os.path.join(diamond_dir, marker), db).split())

    total_markers = len(glob(diamond_dir + '/*tsv'))

    all_markers = set()
    all_pident = defaultdict(list)
    for f1 in glob(diamond_dir + '/*tsv'):

        marker = os.path.basename(f1).split('.')[0]
        all_markers.add(marker)

        seq_mapping = defaultdict(list)
        with open('INPUT/unique_markers/' + marker + '.unique.tsv') as infile:
            infile.readline()
            for line in infile:
                line = line.strip().split('\t')
                seq_mapping[line[0]].append(line[1])

        with open(diamond_dir + '/' + marker + '.dmd.tsv') as infile:
            for line in infile:
                qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split(
                    '\t')
                pident, evalue = float(pident), float(evalue)
                if evalue >= 1e-4:
                    continue
                for genome in seq_mapping[sseqid]:
                    all_pident[genome].append((marker, pident))

    outfile = open(os.path.join(outdir, 'closest_genomes.tsv'), 'w')

    print(*(['genome', 'num_markers', 'total_markers', 'avg_pident'] + [marker for marker in all_markers]), sep='\t',
          file=outfile)
    closest_genomes = []
    for genome in all_pident:

        taxid = genome2taxid[int(genome)]

        if len(all_pident[genome]) / float(total_markers) < 0.25:
            continue

        marker_pident = dict(all_pident[genome])
        pidents = []
        for marker in all_markers:
            try:
                pidents.append(marker_pident[marker])
            except:
                pidents.append(None)
        closest_genomes.append(
            [genome, taxid, taxon2species[taxid][0], taxon2species[taxid][1], len(all_pident[genome]), total_markers,
             np.mean(list(marker_pident.values()))] + pidents)

    closest_genomes = list(reversed(sorted(closest_genomes, key=lambda x: x[3])))

    for res in closest_genomes:
        print(*res, sep='\t', file=outfile)

    outfile.close()