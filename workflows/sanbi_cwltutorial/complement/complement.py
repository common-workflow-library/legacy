#!/usr/bin/env python

import sys
import click

@click.command()
@click.argument('infile', type=click.File('rb'))
@click.argument('outfile', type=click.File('wb'), required=False, default=sys.stdout)
def complement(infile, outfile):
    dnastring = infile.read().upper()
    c_bases = dict(A='T', T='A', C='G', G='C')
    comp_dnastring = ''.join([c_bases.get(base, base) for base in dnastring])
    outfile.write(comp_dnastring)

if __name__ == '__main__':
    complement()
