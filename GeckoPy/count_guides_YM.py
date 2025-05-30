from Bio import SeqIO
import csv
from collections import OrderedDict
import numpy as np
import argparse

KEY_REGION_START = 30  # start index of key region
KEY_REGION_END = 55    # end index of key region
KEY = "CGAAACACC"      # identifies sequence before guide to determine guide position

def count_spacers(input_file, fastq_file, output_file, guide_g): 
    """
    Parses a FASTQ file to count occurrences of guide sequences defined in a library file.

    Parameters:
    - input_file: CSV file with guide sequences (one per line).
    - fastq_file: FASTQ file with sequencing reads.
    - output_file: Output CSV file to store counts.
    - guide_g: Boolean indicating if a G is appended to the key before the spacer.
    """
    num_reads = 0
    perfect_matches = 0
    non_perfect_matches = 0
    key_not_found = 0

    key_seq = KEY + ("G" if guide_g else "")

    try:
        with open(input_file, newline='') as infile:
            reader = csv.reader(infile)
            dictionary = {rows[0]: 0 for rows in reader if rows}
    except Exception as e:
        print(f'Could not open input file: {input_file}')
        print(str(e))
        return

    try:
        handle = open(fastq_file, "r")
    except Exception as e:
        print(f"Could not open FASTQ file: {fastq_file}")
        print(str(e))
        return

    readiter = SeqIO.parse(handle, "fastq")
    for record in readiter:
        num_reads += 1
        read_sequence = str(record.seq).upper()
        key_region = read_sequence[KEY_REGION_START:KEY_REGION_END]
        key_index = key_region.find(key_seq)
        if key_index >= 0:
            start_index = key_index + KEY_REGION_START + len(key_seq)
            guide = read_sequence[start_index:(start_index + 20)]
            if guide in dictionary:
                dictionary[guide] += 1
                perfect_matches += 1
            else:
                non_perfect_matches += 1
        else:
            key_not_found += 1

    dict_sorted = OrderedDict(sorted(dictionary.items()))
    with open(output_file, 'w', newline='') as csvfile:
        mywriter = csv.writer(csvfile)
        mywriter.writerow(["Guide", "Count"])
        for guide, count in dict_sorted.items():
            mywriter.writerow([guide, count])

    values = list(dictionary.values())
    percent_matched = round(perfect_matches / float(perfect_matches + non_perfect_matches) * 100, 1) \
                      if (perfect_matches + non_perfect_matches) > 0 else 0
    guides_with_reads = np.count_nonzero(values)
    guides_no_reads = len(values) - guides_with_reads
    percent_no_reads = round(guides_no_reads / float(len(values)) * 100, 1) if values else 0
    try:
        top_10 = np.percentile(values, 90)
        bottom_10 = np.percentile(values, 10)
        skew_ratio = round(top_10 / bottom_10, 2) if bottom_10 > 0 else 'Inf'
    except:
        skew_ratio = 'Not enough data to calculate skew ratio'

    with open('statistics.txt', 'w') as statsfile:
        statsfile.write(f'Number of perfect guide matches: {perfect_matches}\n')
        statsfile.write(f'Number of nonperfect guide matches: {non_perfect_matches}\n')
        statsfile.write(f'Number of reads where key was not found: {key_not_found}\n')
        statsfile.write(f'Number of reads processed: {num_reads}\n')
        statsfile.write(f'Percentage of guides that matched perfectly: {percent_matched}%\n')
        statsfile.write(f'Percentage of undetected guides: {percent_no_reads}%\n')
        statsfile.write(f'Skew ratio of top 10% to bottom 10%: {skew_ratio}\n')

    handle.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze sequencing data for sgRNA library distribution')
    parser.add_argument('-f', '--fastq', type=str, dest='fastq_file',
                        help='FASTQ file name', default='NGS.fastq')
    parser.add_argument('-o', '--output', type=str, dest='output_file',
                        help='Output CSV file name', default='library_count.csv')
    parser.add_argument('-i', '--input', type=str, dest='input_file',
                        help='Input guide library CSV file', default='library_sequences.csv')
    parser.add_argument('-no-g', dest='guide_g', help='Disable guanine before spacer', action='store_false')
    parser.set_defaults(guide_g=True)
    args = parser.parse_args()

    count_spacers(args.input_file, args.fastq_file, args.output_file, args.guide_g)
