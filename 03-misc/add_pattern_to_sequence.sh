#!/bin/bash
## use this script to build an index/barcode file from a fastq file.
## it will replace the seq & qual rows with the pattern provided, keeping sequence header unchanged
## use R1 for I1 index and R2 for I2 index

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <pattern> <action> <input_fastq_file>"
    echo "Actions: replace, prefix, suffix"
    exit 1
fi

PATTERN=$1
ACTION=$2
INPUT_FILE=$3


# Use AWK to modify the sequence based on the specified action
awk -v pattern="$PATTERN" -v action="$ACTION" '

BEGIN {
    # Seed the random number generator
    srand();

    # Define IUPAC ambiguity codes
    iupac["A"] = "A";
    iupac["C"] = "C";
    iupac["G"] = "G";
    iupac["T"] = "T";
    iupac["R"] = "AG";    # G or A
    iupac["Y"] = "CT";    # C or T
    iupac["S"] = "GC";    # G or C
    iupac["W"] = "AT";    # A or T
    iupac["K"] = "GT";    # G or T
    iupac["M"] = "AC";    # A or C
    iupac["B"] = "CGT";   # C, G, or T
    iupac["D"] = "AGT";   # A, G, or T
    iupac["H"] = "ACT";   # A, C, or T
    iupac["V"] = "ACG";   # A, C, or G
    iupac["N"] = "ACGT";  # Any base
}

{
    if (NR % 4 == 2) {
        # This is a sequence line
        original_seq = $0;
        new_seq = "";
        for (i = 1; i <= length(pattern); i++) {
            char = substr(pattern, i, 1);
            if (char in iupac) {
                # Choose a random character from the IUPAC set
                options = iupac[char];
                chosen_char = substr(options, int(rand() * length(options)) + 1, 1);
                new_seq = new_seq chosen_char;
            } else {
                # If the character is not an IUPAC code, use it as is
                new_seq = new_seq char;
            }
        }

        # Apply the specified action
        if (action == "replace") {
            $0 = new_seq;
        } else if (action == "prefix") {
            $0 = new_seq original_seq;
        } else if (action == "suffix") {
            $0 = original_seq new_seq;
        }
        # Print the line
        print $0;

    } else if (NR % 4 == 0) {
        # This is a quality line, replace it with "I" repeated for the length of the modified sequence
        if (action == "replace") {
            print substr(quality, 1, length(pattern));
        } else if (action == "prefix" || action == "suffix") {
            original_quality = $0;
            quality_for_pattern = substr(quality, 1, length(pattern));
            if (action == "prefix") {
                $0 = quality_for_pattern original_quality;
            } else {
                $0 = original_quality quality_for_pattern;
            }
            # Print the line
            print $0;
        }
    } else{
        # Print the line
        print $0;
    }

}

' quality=$(printf "I%.0s" $(seq 1 ${#PATTERN})) "$INPUT_FILE" 
