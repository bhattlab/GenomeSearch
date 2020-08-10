import subprocess
from os.path import join
from os import remove, devnull
from Bio import SeqIO
import gzip
import math
from multiprocessing import Pool
from itertools import cycle


def run_prodigal(prodigal_path, infile, outdir, meta):
    bashCommand = '{prodigal} -c -f gff -o {gff} -a {faa} -d {ffn} -i {input}'.format(
        prodigal=prodigal_path, gff=join(outdir, 'prodigal.gff'), faa=join(outdir, 'prodigal.faa'),
        ffn=join(outdir, 'prodigal.ffn'), input=infile
    )
    if meta == True:
        bashCommand += ' -p meta'
    print('prodigal command:', bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    output, error = process.communicate()

def run_prodigal_simple(prodigal_path, infile, outprefix):
    FNULL = open(devnull, 'w')
    bashCommand = '{prodigal} -p meta -c -f gff -o {gff} -a {faa} -d {ffn} -i {input}'.format(
        prodigal=prodigal_path, gff=outprefix + '.gff', faa=outprefix + '.faa',
        ffn=outprefix + '.ffn', input=infile
    )

    print('prodigal command:', bashCommand)
    process = subprocess.Popen(bashCommand.split(), stdout=FNULL, stderr=FNULL)
    output, error = process.communicate()

def run_prodigal_multithread(prodigal_path, infile, outdir, threads):

    print("Counting records in FASTA file...")
    record_count = count_records_in_fasta(infile)
    print("The FASTA file contains %d records..." % record_count)

    print("Writing FASTA file to batches for multithreading..." )
    record_count_per_file = math.ceil(record_count / threads)
    filecount = 0
    outrecs = []
    outfiles = []

    if infile.endswith('.gz'):
        infile = gzip.open(infile, "rt")

    for record in SeqIO.parse(infile, 'fasta'):
        outrecs.append(record)
        if len(outrecs) == record_count_per_file:
            filecount += 1
            outfile = join(outdir, 'input%d.fna' % filecount)
            SeqIO.write(outrecs, outfile, 'fasta')
            outfiles.append(outfile)
            outrecs = []

    if len(outrecs) > 0:
        filecount += 1
        outfile = join(outdir, 'input%d.fna' % filecount)
        SeqIO.write(outrecs, outfile, 'fasta')
        outfiles.append(outfile)
    del outrecs

    prodigal_files = [join(outdir, 'prodigal%d' % (i+1)) for i in range(len(outfiles))]
    with Pool(processes=threads) as pool:
        pool.starmap(run_prodigal_simple, zip(cycle([prodigal_path]), outfiles, prodigal_files))

    combine_files(outfiles, join(outdir, 'input.fna'))
    combine_files([f+'.faa' for f in prodigal_files], join(outdir, 'prodigal.faa'))
    combine_files([f+'.ffn' for f in prodigal_files], join(outdir, 'prodigal.ffn'))
    combine_files([f+'.gff' for f in prodigal_files], join(outdir, 'prodigal.gff'), ['#'])


def combine_files(files, outfile, exclude_startswith=[]):
    with open(outfile, 'w') as out:

        for f in files:

            with open(f) as infile:
                for line in infile:

                    skip = False
                    for exclude in exclude_startswith:
                        if line.startswith(exclude):
                            skip = True
                            break
                    if skip is True:
                        continue
                    out.write(line)

    for f in files:
        remove(f)

def count_records_in_fasta(fasta):
    records = 0
    if fasta.endswith('.gz'):
        with gzip.open(fasta, "rt") as infile:
            for line in infile:
                if line.startswith('>'):
                    records += 1
    else:
        with open(fasta) as infile:
            for line in infile:
                if line.startswith('>'):
                    records += 1
    return records