import os, sys, csv, re


def parse_hisat2_summary(filepath):
    """
    Parse relevant information from a HISAT2 summary file (_h2report.txt).
    Handles paired-end format:
        12345 reads; of these:
          12345 (100.00%) were paired; of these:
            123 (1.00%) aligned concordantly 0 times
            11111 (90.00%) aligned concordantly exactly 1 time
            1111 (9.00%) aligned concordantly >1 times
            ----
            123 pairs aligned concordantly 0 times; of these:
              12 (9.76%) aligned discordantly 1 time
            ----
            111 pairs aligned 0 times concordantly or discordantly; of these:
              222 mates make up the pairs; of these:
                111 (50.00%) aligned 0 times
                56 (25.23%) aligned exactly 1 time
                55 (24.77%) aligned >1 times
        95.50% overall alignment rate
    """
    data = dict()
    with open(filepath, 'r') as f:
        content = f.read()

    # Total reads/pairs
    match = re.search(r'^(\d+) reads; of these:', content, re.MULTILINE)
    if match:
        data["Total reads"] = match.group(1)

    # Were paired
    match = re.search(r'(\d+) \(([\d.]+)%\) were paired', content)
    if match:
        data["Paired reads"] = match.group(1)
        data["Paired reads %"] = match.group(2) + "%"

    # Concordantly aligned 0 times
    match = re.search(r'(\d+) \(([\d.]+)%\) aligned concordantly 0 times', content)
    if match:
        data["Aligned concordantly 0 times"] = match.group(1)
        data["Aligned concordantly 0 times %"] = match.group(2) + "%"

    # Concordantly aligned exactly 1 time
    match = re.search(r'(\d+) \(([\d.]+)%\) aligned concordantly exactly 1 time', content)
    if match:
        data["Aligned concordantly exactly 1 time"] = match.group(1)
        data["Aligned concordantly exactly 1 time %"] = match.group(2) + "%"

    # Concordantly aligned >1 times
    match = re.search(r'(\d+) \(([\d.]+)%\) aligned concordantly >1 times', content)
    if match:
        data["Aligned concordantly >1 times"] = match.group(1)
        data["Aligned concordantly >1 times %"] = match.group(2) + "%"

    # Aligned discordantly 1 time
    match = re.search(r'(\d+) \(([\d.]+)%\) aligned discordantly 1 time', content)
    if match:
        data["Aligned discordantly 1 time"] = match.group(1)
        data["Aligned discordantly 1 time %"] = match.group(2) + "%"

    # Overall alignment rate
    match = re.search(r'([\d.]+)% overall alignment rate', content)
    if match:
        data["Overall alignment rate"] = match.group(1) + "%"

    return data


def parse_all_summaries(directory):
    """
    Process all HISAT2 summary files (*_h2report.txt) in the directory and collect results.
    """
    final_data = list()
    header = None

    for filename in sorted(os.listdir(directory)):
        if filename.endswith("_h2report.txt"):
            filepath = os.path.join(directory, filename)
            parsed_data = parse_hisat2_summary(filepath)
            sample_name = filename.replace("_h2report.txt", "")
            parsed_data["Sample"] = sample_name
            final_data.append(parsed_data)

            if header is None:
                header = list(parsed_data.keys())

    # Reorder header to move "Sample" to the front
    if header is not None:
        header = ["Sample"] + [col for col in header if col != "Sample"]

    return final_data, header


def write_to_csv(final_data, header, output_csv):
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        for row in final_data:
            writer.writerow(row)


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python generate_map_metrics.py <config_directory> <directory_with_hisat2_summaries>")
        sys.exit(1)

    config_directory = sys.argv[1]
    directory = sys.argv[2]

    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        sys.exit(1)

    final_data, header = parse_all_summaries(directory)

    os.makedirs(os.path.join(config_directory, "3_1_map_metrics_output_qc"), exist_ok=True)
    output_csv = os.path.join(config_directory, "3_1_map_metrics_output_qc", "hisat2_summary.csv")
    write_to_csv(final_data, header, output_csv)
    print(f"Summary CSV file saved as {output_csv}")
