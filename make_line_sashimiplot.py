#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2019.11.28

chr2    HAVANA  transcript      218781749       218815288

python main.py line \
    -b all_bam.txt \
    -e chr2:218780000-218786749:+ \
    -g gencode.v32.annotation.gtf  \
    --color-factor 3 \
    -o ENST00000494263_ATAC_not_y.pdf \
    -p 20 \
    --share-y \
    --plot-by 4 \
    --sep-by-color
"""
import os
from multiprocessing import Pool
from subprocess import check_call, CalledProcessError


SCRIPT = "pysashimi/main.py"


class Gene(object):

    def __init__(self, chrom, start, end, strand, gene):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.gene = gene
        self.strand = strand

    def __hash__(self):
        return hash(self.gene)

    def promoter(self, up = 2000, down = 500):
        u"""
        return a promoter region
        """

        return "{}:{}-{}:{}".format(
            self.chrom,
            self.start - up if self.strand == "+" else self.end + up,
            self.start + down if self.strand == "+" else self.end - down,
            self.strand
        )


def load_reference(bed):
    u"""
    Load reference
    :param bed: path to a bed file
    """
    print("Loading %s" % bed)
    data = {}
    with open(bed) as r:
        for line in r:
            if line.startswith("M"):
                continue

            lines = line.split()

            data[lines[3]] = Gene(lines[0], lines[1], lines[2], lines[5], lines[3])

    return data


def call(cmd):
    try:
        with open(os.devnull, "w+") as w:
            check_call(cmd, shell = True, stdout = w, stderr = w)
    except CalledProcessError as err:
        print(err)


def main(bed, sample, bam_list, gtf, output, n_jobs = 10):
    u"""
    Main
    """
    if not os.path.exists(output):
        os.makedirs(output)
    print("Reading sample")

    genes = sample.split(",")

    reference = load_reference(bed)

    tasks = []
    for i in genes:
        if i not in reference.keys():
            continue
        print(i)
        tasks.append("python pysashimi/main.py normal -e {} -b {} -g {} --color-factor 3 -o {} -p {} --remove-empty-gene --no-gene --title {} --share-y".format(  # --plot-by 4  --sep-by-color --distance-ratio 0.35
            reference[i].promoter(),
            bam_list, gtf,
            os.path.join(output, "{}.txt".format(i)), n_jobs,
            i
        ))


    with Pool(min(n_jobs, 6)) as p:
        p.map(call, tasks)


if __name__ == '__main__':
    from fire import Fire
    Fire(main)

   
