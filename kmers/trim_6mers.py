def trim_file(input_filename, output_filename):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            columns = line.split()  # Assuming columns are separated by spaces
            if len(columns) >= 2:
                outfile.write(columns[0] + ' ' + columns[1] + '\n')

input_filename = '9.4_6mers_450bps.txt'  # Replace with your input file
output_filename = '9.4_6mers_450bps_trimmed.txt'  # Replace with your desired output file name
trim_file(input_filename, output_filename)

