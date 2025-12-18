import argparse
import re

# Load the DeepTE.fa information into a dictionary
def load_deepte_fa(fasta_file):
    annotations = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                # Extract sequence range and annotation from the header
                header = line.strip()
                match = re.match(r'>[^:]+:(\d+)-\d+#(.+)', header)
                if match:
                    start = int(match.group(1)) + 1  # Add 1 to the start position
                    annotation = match.group(2)
                    annotations[start] = annotation  # Store annotation by the start position
    return annotations

# Update the test.out file based on DeepTE.fa annotations
def update_test_out(out_file, annotations, output_file):
    with open(out_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) > 10 and 'Unknown' in fields[10]:
                try:
                    # The sixth column corresponds to the start position in test.out
                    start = int(fields[5])  # Sixth column in test.out
                    if start in annotations:
                        # Replace "Unknown" in the eleventh column with the annotation
                        fields[10] = annotations[start]
                        print(f"Updated line with start {start}: {annotations[start]}")
                    else:
                        print(f"No match found for start {start}")
                except ValueError:
                    pass
            outfile.write('\t'.join(fields) + '\n')

# Main function to handle command-line arguments
def main():
    parser = argparse.ArgumentParser(description="Update annotations in test.out based on DeepTE.fa")
    parser.add_argument('-inf', '--input_fasta', required=True, help='Path to the DeepTE.fa file')
    parser.add_argument('-i', '--input_out', required=True, help='Path to the test.out file')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file where updates will be saved')

    args = parser.parse_args()

    # Load annotations from DeepTE.fa
    annotations = load_deepte_fa(args.input_fasta)

    # Update test.out based on annotations
    update_test_out(args.input_out, annotations, args.output)

    print(f"Updated {args.input_out} saved to {args.output}")

if __name__ == "__main__":
    main()
