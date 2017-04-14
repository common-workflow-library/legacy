#!/usr/bin/env python

import sys
import click

@click.command()
@click.argument('infile', type=click.File('rb'))
@click.argument('outfile', type=click.File('wb'), required=False, default=sys.stdout)
def reverse(infile, outfile):
    dnastring = infile.read()
    if dnastring.endswith('\n'):
        add_eol = True
        dnastring = dnastring.rstrip()
    else:
        add_eol = False
    reverse_string = dnastring[::-1]
    if add_eol:
        reverse_string += '\n'
    outfile.write(reverse_string)

if __name__ == '__main__':
    reverse()
