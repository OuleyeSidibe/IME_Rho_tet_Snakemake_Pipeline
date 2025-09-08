## Librairies
import csv
import pandas as pd 
import os


#############################################################################################################
# ## Notes 

"""  Parse FIMO outputs to extract TIR which match with IME_Rho_tet positions and orientations """

# ## Load script  in work directory with appropriate input files

#############################################################################################################


# Functions
def process_files():
    # Open the CSV file "RefseqFINAL3.csv)"
    try:
        with open("RefseqFINAL3.csv", newline='') as csvfile, open("parse_outIME_Rho_tet.csv", "a") as parse_out:
            csv_reader = csv.reader(csvfile)
            
            
            # print("ID strand TIRamont_s TIRamont_e TIRamont_seq TIRaval_s TIRaval_e TIRaval_seq mge_size TIRamont_s2 TIRamont_e2 TIRamont_seq2")
            for row in csv_reader:
                # Check if row has enough columns (at least 23 columns for index 22)
                if len(row) < 23:
                    continue
                # Only process rows where column 23 (index 22) equals "TetMGE"
                if row[22] != 'TetMGE':
                    continue

                # Extract required columns from CSV (1-indexed columns used in instructions)
                # Column 1: ID, Column 2: ACC, Column 9: NUCL, Column 13: REL,
                # Column 16: REC, Column 21: MGE
                ID = row[0]
                ACC = row[1]
                NUCL = row[8]
           
                REL_str = row[12]
                REC_str = row[15]

                MGE_str = row[20]

                # Initialize variables that will be set based on conditions in TSV file.
                amont_s = None
                amont_e = None
                amont_seq = None
                amont_s2 = None
                amont_e2 = None
                amont_seq2 = None
                aval_s = None
                aval_e = None
                aval_seq = None
                mge_size = None

                # Parse REC, REL, and MGE ranges (format 'start..end')
                try:
                    rec_parts = REC_str.split("..")
                    rec_start = int(rec_parts[0])
                    rec_end   = int(rec_parts[1])

                    rel_parts = REL_str.split("..")
                    rel_start = int(rel_parts[0])
                    rel_end   = int(rel_parts[1])

                    mge_parts = MGE_str.split("..")
                    mge_start = int(mge_parts[0])
                    mge_end   = int(mge_parts[1])
                    # Some MGE positions are reverse (when the TET gene is in the opposite strand)
                    if mge_start > mge_end:
                     mge_start = int(mge_parts[1])
                     mge_end   = int(mge_parts[0])
                     
                    
                except (IndexError, ValueError):
                    # If parsing fails, skip this row
                    print(ID, strand, amont_s, amont_e, amont_seq, aval_s, aval_e, aval_seq, mge_size, amont_s2, amont_e2, amont_seq2, "error in parsing positions")
                    continue

                # Determine STRAND based on comparison of REC start and REL start
                ##for TetMGE groupe
                if rec_start > rel_start:
                    strand = "+"
                else:
                    strand = "-"
                

                # Construct the filename for the TSV file based on the ACC variable
                acc_filename = f"/home/sbleclercq/pgba_work/sebastien/tetPPR/TIRdetections/results/{ACC}.tsv"
              
                if not os.path.isfile(acc_filename):
                    # Skip processing if the file does not exist
                    print(ID, strand, amont_s, amont_e, amont_seq, aval_s, aval_e, aval_seq, mge_size, amont_s2, amont_e2, amont_seq2, "FIMO file does not exist")
                    continue

                try:
                    with open(acc_filename, newline='') as tsvfile:
                        tsv_reader = csv.reader(tsvfile, delimiter='\t')
                        for tsv_row in tsv_reader:
                            # Ensure the row has enough columns (we require at least 10 columns)
                            if len(tsv_row) < 10:
                                continue

                            # Extract values from specific columns:
                            # Column 3 => index 2: remove the last two characters
                            n2_raw = tsv_row[2]
                            n2 = n2_raw[:-2] if len(n2_raw) >= 2 else ""

                            try:
                                # Column 4 and Column 5 => T_START and T_END
                                t_start = int(tsv_row[3])
                                t_end   = int(tsv_row[4])
                            except ValueError:
                                # If conversion fails, skip this row
                                continue

                            # Column 6 => T_STR
                            t_str = tsv_row[5]
                            # Column 10 => T_SEQ (index 9)
                            t_seq = tsv_row[9]

                            # Process conditions only if n2 matches the NUCL value from CSV
                            if n2 == NUCL:
                                # First conditional block (for T_START conditions)
                                if (t_str == '+' and t_start < mge_start and t_start > (mge_start - 6000)):
                                    if strand == '+':
                                      if amont_seq == None: # we may have more than one correct TIR in amont
                                        amont_s = t_start
                                        amont_e = t_end
                                        amont_seq = t_seq
                                      else:
                                        amont_s2 = t_start
                                        amont_e2 = t_end
                                        amont_seq2 = t_seq
                                    else:
                                        aval_s = t_start
                                        aval_e = t_end
                                        aval_seq = t_seq

                                # Second conditional block (for T_END conditions)
                                if (t_str == '-' and t_end > mge_end and t_end < (mge_end + 6000)):
                                    if strand == '-':
                                      if amont_seq == None: # we may have more than one correct TIR in amont
                                        amont_s = t_start
                                        amont_e = t_end
                                        amont_seq = t_seq
                                      else: 
                                        # When we are in reverse strand, we want the last TIR found to be the correct one, so we need to shift the values
                                        amont_s2 = amont_s
                                        amont_e2 = amont_e
                                        amont_seq2 = amont_seq
                                        amont_s = t_start
                                        amont_e = t_end
                                        amont_seq = t_seq
                                          
                                    else:
                                        aval_s = t_start
                                        aval_e = t_end
                                        aval_seq = t_seq
                except Exception as e:
                    # If any error occurs processing the TSV file, skip to the next CSV row.
                    print(ID, strand, amont_s, amont_e, amont_seq, aval_s, aval_e, aval_seq, mge_size, amont_s2, amont_e2, amont_seq2, str(e))
                    continue

                # Check if both AMONT and AVAL values were assigned before computing MGE_SIZE
                if amont_s is None or aval_e is None:
                    # If either assignment is missing, skip this entry.
                    mge_size = "Unknown"
                else:
                 if strand == '+':
                    mge_size = aval_e - amont_s
                 else:    
                    mge_size = amont_e - aval_s

                # Print the required output: ID, STRAND, AMONT_S, AMONT_E, AMONT_SEQ, 
                # AVAL_S, AVAL_E, AVAL_SEQ, MGE_SIZE
                parse_out.write(f"{ID}, {NUCL}, {MGE_str}, {amont_s}, {amont_e}, {amont_seq}, {aval_s}, {aval_e}, {aval_seq}, {mge_size}, {amont_s2}, {amont_e2}, {amont_seq2}\n")
                
                

    except FileNotFoundError:
        print("The file 'RefseqFINAL3.csv' was not found.")


print("Process Done !!")
if __name__ == "__main__":
    process_files()

