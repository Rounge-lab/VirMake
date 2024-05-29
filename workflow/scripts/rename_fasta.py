#!/usr/bin/env python

import argparse


def main():
    # Parse arguments.
    parser = argparse.ArgumentParser(description='Rename FASTA files.', epilog='Work out those contigs.')
    parser.add_argument('-i', '--input', help='specify input FASTA file', required=True)
    parser.add_argument('--pre', help='string pre contig count', type=str, default='')
    parser.add_argument('--pos', help='string post contig count', type=str, default='')
    parser.add_argument('-o', '--output', help='specify output FASTA file', required=True)
    parser.add_argument('--int_length', help = 'length of number to make output of uniform length', type = int, default = 5)
    parser.add_argument('--contig_id_file', help = 'old to new conversion file', required = True)
    parser.add_argument('-s', help = 'first number for vOTU rename', default=1)

    args=parser.parse_args()

    # Open FASTA.
    fasta_in = open(args.input)

    # Create FASTA output file.
    fasta_out = open(args.output, 'w')

    # Start counter.
    count = int(args.s)

    old_names, new_names = [], []

    if args.pos == "None":
        pos = ""
    else:
        pos = args.pos

    # Parse file and write to output
    for line in fasta_in.readlines():
        if line.startswith('>'):
            contig_id = '>' + args.pre + str(count).zfill(args.int_length) + pos + '\n'
            fasta_out.write(contig_id)
            old_names.append(line[1:].strip())
            new_names.append(contig_id[1:].strip())
            count += 1
        else:
            fasta_out.write(line)

     # Finish.
    fasta_out.close()
    fasta_in.close()

    # Write file linking new to old ids
    with open(args.contig_id_file, 'w') as f:
            f.write("old_id\tnew_id\n")
            for i in range(0,len(old_names)):
                f.write(old_names[i]+"\t"+new_names[i]+"\n")
            f.close()


if __name__ == '__main__':
    main()
