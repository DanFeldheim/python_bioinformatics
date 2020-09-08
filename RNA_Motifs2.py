# This program imports merged and trimmed fasta files and performs sequence alignment
# The alignment includes a bead database such that bead binders are aligned with the sequences of interest
# Output is an html file with families and sequences that appear in families that did not co-align with bead binders
# Sequences are checked against the bead database one more to find matches


# Import packages
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import regex
import os
import glob
import csv



    
class Motif:
    
    def __init__(self):
        
        # Enter path to import merged and trimmed fasta file for alignment
        self.file_path = '/Users/danfeldheim/Documents/test_files/'
        
        # Enter path for alignment results html file
        self.results_path = '/Users/danfeldheim/Documents/'
        
        # Enter path and filename to export csv file of most abundant sequences that were found in a family with no relation to bead binders
        self.most_abundant_from_py_aligment = '/Users/danfeldheim/Documents/most_abundant_from_py_aligment.csv'
        
        # Enter True on the right hand side if you wish to export a list of motifs to the database of all motifs for all selections
        # If not enter False
        # This will also export a list of motifs for each individual sequencing file
        self.add_motifs = True
        # Enter path and filename of master motif database
        self.motif_db_filename = '/Users/danfeldheim/Documents/motif_master_list.csv'
        
        # Enter length of motif to search for
        self.motif_length = 13
        
        # Enter allowed number of substitutions in motif
        self.subs = 1
        
        # Enter shortest length sequence allowed (note average length = 40)
        self.shortest_length = 40
        
        # Enter longest length sequence allowed (note average length = 40)
        self.longest_length = 40
    
        # Input minimum desired family size
        self.min_fam_size = 3
        
        # Enter a bead database to blast in EITHER fasta or csv format
        # Enter path and file for master bead database in csv format (preferred)
        self.beaddb = '/Users/danfeldheim/Documents/Bead_Database/master_bead_binder_database.csv'
        # self.beaddb = '/Users/danfeldheim/Documents/Bead_Database/bead_test2.csv'
        
        # Enter path and file to bead database for bead database in fasta format
        # self.beaddb = '/Users/danfeldheim/Documents/Bead_Database/NP08-9-8_R1_001_041520_mergedandtrimmed_allsense.fasta'
        # self.beaddb = '/Users/danfeldheim/Documents/bead_test2.fasta'
        
    def import_files(self):
        
        # Call glob to place all files a list for processing
        # self.import_files = sorted(glob.glob(self.file_path + '*.fasta'))
        self.import_files = sorted(glob.glob(self.file_path + 'NP01-3-A01_R1_001.fastq_merged_and_trimmed.fasta'))
        # self.import_files = sorted(glob.glob(self.file_path + 'traveling_motif_test.fasta'))
        
        # Loop through files and create dataframe
        
        for file in self.import_files:
            
            # Use path and filename to create an output filename for motifs
            # Separate filename and path
            self.extract_filename = os.path.basename(file)
                
            # Strip out the extension
            self.base_filename = self.extract_filename.rpartition('.')[0]
            self.new_name, self.sep, self.tail = self.base_filename.partition('.')
            
            # Add the motif length to the filename
            self.new_name = self.new_name + "_" + str(self.motif_length) + "mer"
            
            # Create a path and filename for results
            self.results_path_name = os.path.join(self.results_path, self.new_name + ".html")
            
            # Read fasta file
            with open(file, mode ='r') as handle:
        
                # Create lists for sequences 
                self.seq_list = []
                
                # Create lists for sequence IDs
                self.seqID_list = []
            
                # Use Biopython's parse function to process individual FASTA records
                for record in SeqIO.parse(handle, 'fasta'):
                
                    # Extract sequence
                    self.sequence = str(record.seq)
                    
                    # Extract sequence ID
                    self.sequence_id = str(record.id)
        
                    # Append sequence to list
                    self.seq_list.append(self.sequence)
                    
                    # Append sequence ID to list
                    self.seqID_list.append(self.sequence_id)
        
            # Create dataframe with identifier and sequence
            self.sequence_df = pd.DataFrame(columns = ['Sequence', 'Sequence_ID'], dtype = 'str')
        
            # Append seq_list to dataframe
            self.sequence_df['Sequence'] = self.seq_list
            self.sequence_df['Sequence_ID'] = self.seqID_list
            
            # Get frequency of every sequence in df if desired
            # self.sequence_df['freq'] = self.sequence_df.Sequence.map(self.sequence_df.Sequence.value_counts())
            # print ('freq:')
            # print(self.sequence_df.head())
            
            # Call other functions so that loop can be run on many files
            
            # Send to self.clean to remove replicates and sequences with Ns
            x = self.sequence_df
            self.clean(x)
            
            # Call bead_csv OR import_beaddb to load the bead database
            # bead_csv imports bead database from csv file, import_beaddb imports from fasta file
            self.bead_csv()
            # self.import_beaddb
            
            # Join sequence and bead dataframes
            self.join_df()
            
            # Create motifs
            self.motifs()
            
            # Align into families
            self.families()
            
            # Export html file
            self.create_html()
            
            # Delete bead families
            self.no_bead_families()
            
            # Determine motif travel parameter
            self.motif_travel()
            
            # Calculate abundance of top sequences
            self.abundance()
            
            # Confirm that the top sequences are not bead sequences
            self.bead_check()
            
            # Append most abundant sequences to html file
            self.append_top_seqs()
            
            # Export top sequences to csv file
            self.export_csv()
            
            # Export motifs to motif database
            if self.add_motifs:
                self.export_motifs()
            
 
    # Function to read in a bead database in fasta format     
    # Should be merged and trimmed since it will be combined with merged and trimmed sequence file
    # Can also read in a csv file of bead binders using self.bead_csv below                  
    def import_beaddb(self):
        
        # Read bead fasta file
            with open(self.beaddb, mode ='r') as handle:
        
                # Create lists for sequences 
                self.bead_seq_list = []
                
                # Create lists for sequence IDs
                self.bead_seqID_list = []
            
                # Use Biopython's parse function to process individual FASTA records
                for record in SeqIO.parse(handle, 'fasta'):
                
                    # Extract sequence
                    self.bead_seq = str(record.seq)
                    
                    # Extract sequence ID
                    self.bead_seq_id = str(record.id)
        
                    # Append sequence to list
                    self.bead_seq_list.append(self.bead_seq)
                    
                    # Append sequence ID to list
                    self.bead_seqID_list.append(self.bead_seq_id)
        
            # Create dataframe with identifier and sequence
            self.bead_sequence_df = pd.DataFrame(columns = ['Sequence', 'Sequence_ID'], dtype = 'str')
        
            # Append seq_list to dataframe
            self.bead_sequence_df['Sequence'] = self.bead_seq_list
            self.bead_sequence_df['Sequence_ID'] = self.bead_seqID_list
            
            # Send to self.clean to remove replicate sequences
            # Create the generic x so self.clean can be called from other functions
            y = self.bead_sequence_df
            self.clean(y)
            
    # Function to read in a bead database in csv format
    # The database contains only unique bead binding sequences
    def bead_csv(self):
        # Read csv
        self.bead_sequence_df = pd.read_csv(self.beaddb, dtype = 'str')
    
    # Function to remove duplicates from self.sequence_df 
    def clean(self, arg):
        # Find sequences with N and delete
        self.cleaned_sequence_df = self.sequence_df[~self.sequence_df.Sequence.str.contains("N")]
        
        # Delete sequences outside the range of lengths specified by user above
        # If shortest and longest length desired are the same (e.g., 40)
        if self.shortest_length == self.longest_length:
            
            self.cleaned_sequence_df = self.cleaned_sequence_df[(self.cleaned_sequence_df.Sequence.str.len() == self.shortest_length)]
        
        else:
            
            self.cleaned_sequence_df = self.cleaned_sequence_df[~(self.cleaned_sequence_df.Sequence.str.len() < self.shortest_length) & ~(self.cleaned_sequence_df.Sequence.str.len() > self.longest_length)]
        
        
        # Get abundance of each sequence (perfect matches) and then re-order from highest to lowest abundance
        # This prevents the most abundant sequences from being filtered out when close matches are removed below
        self.cleaned_sequence_df['freq'] = self.cleaned_sequence_df.Sequence.map(self.cleaned_sequence_df.Sequence.value_counts())
        # print ('freq:')
        # print(self.cleaned_sequence_df.head())
        
        # Order from highest to lowest abundance
        self.cleaned_sequence_df = self.cleaned_sequence_df.sort_values(by = 'freq', ascending = False)
        
        # Remove ALL duplicate sequences
        self.cleaned_sequence_df.drop_duplicates(subset = "Sequence", keep = 'first', inplace = True) 
        
        # print ('cleaned:')
        # print(self.cleaned_sequence_df.head())
        
        # Remove all close matches
        # Create new dataframe
        self.sequence_df_noReps = pd.DataFrame(columns = ['Sequence', 'Sequence_ID', 'freq'], dtype = 'str')
    
        # Move top row from self.cleaned_sequence_df to bottom row of self.sequence_df_noReps until 
        # self.cleaned_sequence_df is empty
        # First, reset index of self.cleaned_sequence_df so that rows are numbered 0,1,2... again.
        self.cleaned_sequence_df = self.cleaned_sequence_df.reset_index(drop = True)
    
        while not self.cleaned_sequence_df.empty:
            # Create a list for matching row index values
            # Create inside the loop so it 'refreshes' as each sequence is checked for a match
            self.match_list = []
            
            # Append row of sequence_df_nodups to sequence_df_cleaned
            # Copy the top row out of sequence_df_nodups
            self.top_row = self.cleaned_sequence_df.iloc[0].copy()
            
            # Append top row to self.sequence_df_noReps
            self.sequence_df_noReps = self.sequence_df_noReps.append(self.top_row)
    
            # Reset row index of cleaned to start at 0
            self.sequence_df_noReps = self.sequence_df_noReps.reset_index(drop = True)
        
            # Get sequence from top of self.cleaned_sequence_df
            self.check_seq = self.top_row['Sequence']
            
            # Set up regex search sequence and parameters (e.g., s <= 4 allows for 4 subs and 0 indels; e for indelsubs)
            self.close_match = ('(%s){s<=2}' % (self.check_seq))
            
            # Use regex.search to find matches with the allowed # of substitutions
            # Loop through self.cleaned_sequence_df comparing search sequence with each sequence in self.sequence_df_noReps
            for row in self.cleaned_sequence_df.itertuples():
                self.comp_seq = row.Sequence
                self.found_one = regex.search(self.close_match, self.comp_seq)
            
                if self.found_one:
                    # Append matching sequences to match_list
                    self.match_list.append(row.Sequence)
        
            # Subset self.cleaned_sequence_df by sequences that do not appear in match_list
            # Note that df.drop did not work here for a number of reasons!
            self.cleaned_sequence_df = self.cleaned_sequence_df[~self.cleaned_sequence_df['Sequence'].isin(self.match_list)]
            
            # Reset row index of self.cleaned_sequence_df so that it starts at 0 again
            self.cleaned_sequence_df = self.cleaned_sequence_df.reset_index(drop = True)
        
        # Get length of self.sequence_df_noReps. This is the number of unique sequences fed into alignment
        self.seqs_for_alignment = len(self.sequence_df_noReps)
        print ('Number of sequences for alignment:')
        print (self.seqs_for_alignment)
        
        # print ('noReps')
        # print (self.sequence_df_noReps.head())
            

    def join_df(self):
        
        # Join sequence and bead dataframes
        self.sequence_df_merged = self.sequence_df_noReps.append(self.bead_sequence_df)
        # Make sure all columns are strings
        self.sequence_df_merged.to_string()
        
        # print ('merged')
        # print (self.sequence_df_merged)
        
            
    def motifs(self):
    
    # Create list for all possible motifs of lenght motif_length
        
        self.motif_list = []
            
        # Iterate down the sequence column in dataframe and split every sequence up into the desired motif length
        for row in self.sequence_df_merged.itertuples():
            self.seq = row.Sequence
            for i in range(0, len(self.seq) - self.motif_length + 1):
                self.motif = self.seq[i:self.motif_length + i]
                self.motif_list.append(self.motif)
            
        # Remove duplicate motifs from motif_list (must use set() and convert result to another list)
        self.motif_list = list(set(self.motif_list))  
        self.motifs_nodups = list(self.motif_list) 
          
        
    def families(self):
       
       # Create a dataframe for families
       self.family_df = pd.DataFrame(columns = ['Motif', 'Sequence', 'Sequence_ID'], dtype = 'str')
       
       # Add self.motifs_nodups to self.family_df
       self.family_df['Motif'] = self.motifs_nodups
       
       # Create lists for all sequences and IDs that match each motif (this will become a list of lists)
       # self.matching_seq_list will house the sequence after coloring the motif
       # This is useful in dataframes where the coloring code behaves as a string instead of a color designation
       self.matching_seq_list = []
       self.matching_seqID_list = []
       self.matching_seq_notColored_list = []
       self.matching_motif_start_index = []
       self.matching_motif_end_index = []
       
       # Compare every motif to every sequence in self.sequence_df_merged
       for row in self.family_df.itertuples():
           # Get motif from self.family_df
           self.motif = row.Motif
           
           # Set up partial match criterion with self.subs = # of allowed substitutions (s for subs only, e for subs, insertions, and deletions)
           self.pattern = ('(%s){s<=%s}' % (self.motif, self.subs))
           
           # Create a list that will hold the sequences that match the current motif
           self.current_matching_seq_list = []
           self.current_matching_seqID_list = []
           self.current_matching_seq_notColored_list = []
           self.current_matching_motif_start_index = []
           self.current_matching_motif_end_index = []
           
           # Create a list of lists to place sequences that contain self.motif
           # Compare motif with sequences in test_sequence_df['Sequence']:
           for row in self.sequence_df_merged.itertuples():
                self.seq_id = row.Sequence_ID
                self.sequence = row.Sequence
                
                # Use regex.search to find matches with the allowed # of substitutions
                self.match = regex.search(self.pattern, self.sequence)
                
                # Extract sequence containing motif, motif start and end positions, and change color of motif within sequence
                if self.match:
                    # Extract starting and ending indeces where the motif appears in the sequence
                    self.match_start = (self.match.start())
                    self.match_end = (self.match.end())   
                    
                    # Change color of motif in the sequence
                    self.seq_color = self.sequence[:self.match_start] + "<span style = \"color:Blue\">" + self.match.group() + "</span>" + self.sequence[self.match_end:]
                    
                    # Append match to self.current_matching_seq_list
                    self.current_matching_seq_list.append(self.seq_color)
                    
                    # Append matching sequence ID to self.current_matching_seqID_list
                    self.current_matching_seqID_list.append(self.seq_id)
                    
                    # Append original non-colored version of matching sequence to self.current_matching_seq_notColored_list
                    self.current_matching_seq_notColored_list.append(self.sequence)
                    
                    # Append motif starting and ending indeces to current lists
                    self.current_matching_motif_start_index.append(self.match_start)
                    self.current_matching_motif_end_index.append(self.match_end)
                    
                 
           # Append all of the matches for the current motif to the master list
           self.matching_seq_list.append(self.current_matching_seq_list)
           
           # Append all of the matching IDs for the current motif to the master list
           self.matching_seqID_list.append(self.current_matching_seqID_list)
           
           # Append all of the matches for the current motif without motif coloring to the master list
           self.matching_seq_notColored_list.append(self.current_matching_seq_notColored_list)
           
           # Append current starting and ending indeces to their master lists
           self.matching_motif_start_index.append(self.current_matching_motif_start_index)
           self.matching_motif_end_index.append(self.current_matching_motif_end_index)
                
       # Append all of the sequence matches (families) for all motifs to self.family_df
       self.family_df['Sequence'] = self.matching_seq_list
       
       # Append all of the sequence matching ID for all motifs to self.family_df
       self.family_df['Sequence_ID'] = self.matching_seqID_list
       
       # Append all of the sequence matching ID for all motifs to self.family_df
       self.family_df['Sequence_no_color'] = self.matching_seq_notColored_list
       
       # Append all of the starting and ending index values for all motifs to self.family_df
       self.family_df['Starting_Index'] = self.matching_motif_start_index
       self.family_df['Ending_Index'] = self.matching_motif_end_index
       
       # Could drop duplicate families based upon sequences in the family
       # If two motifs produce exactly the same family, drop all but one
       # Since 'Sequence' is a list of lists, must recast as a tuple first
       # Would need to create a new column with the original (pre-motif-colored) sequence
       # self.family_df['Sequence_tuple'] = self.family_df['Original_Sequence'].apply(lambda x : tuple(x) if type(x) is list else x)
       # self.family_df.drop_duplicates('Sequence_tuple')
       
       # Delete families less than the desired size (self.min_fam_size)
       self.family_df = self.family_df[self.family_df['Sequence'].str.len() >= self.min_fam_size]
       
       self.family_df.reset_index(drop = True, inplace = True)
       
       # print('family_df:')
       # print(self.family_df)
       
       # Delete families containing only bead sequences
       # Create list to tag all rows as delete or dont delete
       self.all_beads_list = []
       
       # Loop through all sequence IDs and determine if all are beads or not
       for row in self.family_df.itertuples():
           self.all_beads_temp_list = []
  
           for i in range(len(row.Sequence_ID)):
               
               if 'Bead' in row.Sequence_ID[i]:
                   self.all_beads_temp_list.append('yes')
                    
               else:
                   self.all_beads_temp_list.append('no')
                    
           if 'no' in self.all_beads_temp_list:
                self.all_beads_list.append('keep')
                
           else:
                self.all_beads_list.append('delete')
                
           # print(self.all_beads_temp_list)
           
       # Append self.all_beads_list to self.family_df
       self.family_df['All_beads'] = self.all_beads_list
       # print('family_df:')
       # print(self.family_df)
       
       # Delete all rows with delete in All_beads column
       self.family_df = self.family_df[~self.family_df.All_beads.str.contains("delete", na = False)] 
       # print('family_df:')
       # print(self.family_df) 
            
       
    # Function to export families to html
    def create_html(self):
       
       # Start a counter for family number
       self.fam_no = 1 
       # Create a variable for html line break
       self.cr = '<br>'
       self.crlb = '<br><br>'
       
       with open(self.results_path_name, "a") as file:
            
            # Write out number of unique sequences that went into alignment
            self.no_aligned = 'Number of unique sequences aligned: ' + str(self.seqs_for_alignment)
            file.write(self.no_aligned)
            file.write(self.crlb)
            
            file.write('Sequence Families:')
            file.write(self.crlb)
            
            if len(self.family_df .index.values) == 0:
                file.write('All families contained nothing but bead sequences!')
                
            else:
                
                # Indent to continue writing to html file
                for row in self.family_df.itertuples():
                    # Get motif, sequences, and sequence ID
                    self.motif = row.Motif
                    self.seq_family = row.Sequence
                    self.seq_familyID = row.Sequence_ID
                  
                    # Write family number and motif followed by a line break to html
                    self.motif_header = 'Family ' + str(self.fam_no) + ': ' + self.motif
                    file.write(self.motif_header)
                    file.write(self.cr)
                    
                    for i in range(len(self.seq_family)):
                        # Write sequence and sequence ID to html
                        # self.entry = self.seq_family[i] + '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp' +  '  Sequence ID: ' + self.seq_familyID[i]
                        self.entry = self.seq_family[i] + '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp' + self.seq_familyID[i]
                        file.write(self.entry)
                        file.write(self.cr)
                        
                    file.write(self.cr)
                        
                    # Increment the family number
                    self.fam_no += 1
            
    # Function to remove every family that contains a bead binder 
    def no_bead_families(self):
       
       # It would be nice to do this with apply
       self.bead_list = []
       
       for row in self.family_df.itertuples():
           self.temp_list = []
           for j in row.Sequence_ID:
               if 'Bead' in j:
                   self.temp_list.append('bead')
                   
               else:
                   self.temp_list.append('no')
                
           if 'bead' in self.temp_list:
               self.bead_list.append('bead')
               
           else:
               self.bead_list.append('no')
               
       self.family_df['Bead'] = self.bead_list
                
       # Delete rows with bead in self.family_df['Bead']
       self.no_bead_families_df = self.family_df[~self.family_df.Bead.str.contains("bead", na = False)] 
       
       # print('no bead families:')
       # print(self.no_bead_families_df)
       
    # Define a motif 'travelers parameter', which is the how far the motif moves between sequences in a family
    # Travelers parameter = max starting index value - min starting index value
    def motif_travel(self):
        if len(self.no_bead_families_df.index.values) == 0:
            pass
            
        else:
            # Create list for max motif starting index value - min motif starting index value
            self.no_bead_families_df['max_starting_index'] = self.no_bead_families_df.Starting_Index.apply(lambda x: max(x))
            self.no_bead_families_df['min_starting_index'] = self.no_bead_families_df.Starting_Index.apply(lambda x: min(x))
            
            self.no_bead_families_df['max_starting_index'] = self.no_bead_families_df['max_starting_index'].astype(int)
            self.no_bead_families_df['min_starting_index'] = self.no_bead_families_df['min_starting_index'].astype(int)
            
            # Subtract max from min
            self.no_bead_families_df['travel_param'] = self.no_bead_families_df['max_starting_index'] - self.no_bead_families_df['min_starting_index']
            
            # Create a copy of self.no_bead_families_df that is ranked from largest travel_param to smallest
            self.travel_df = self.no_bead_families_df.sort_values(by = 'travel_param', ascending = False) 
            
            # Truncate self.travel_df to the top 6 
            self.travel_df = self.travel_df.head(6)
            
            with open(self.results_path_name, "a") as file:
                file.write('Top 6 Traveling Motifs and Motif Travel Parameter')
                file.write(self.cr)
            
                # Write Motif and travel_param to csv
                if len(self.no_bead_families_df.index.values) == 0:
                    file.write('Could not be calculated.')
                        
                else:
                    for row in self.travel_df.itertuples():
                        # Write motif to html
                        self.motif = row.Motif
                        self.travel_param = row.travel_param
                        self.motif_travel_param = self.motif + ': ' + str(self.travel_param)
                        file.write(self.motif_travel_param)
                        file.write(self.crlb)
                 
            
        # print('travel_df:')
        # print(self.travel_df)
       

    # Function to count the number of copies of each sequence that didn't align with a bead sequence
    def abundance(self):
        
        # Remove duplicate sequences from self.no_bead_families_df
        
        if len(self.no_bead_families_df.index.values) == 0:
            with open(self.results_path_name, "a") as file:
                file.write('All sequences aligned into families containing a bead binder.')
            
        else:    
            # Use the Sequence_no_color column here because the color code makes every sequence in the Sequence column unique!
            self.no_duplicate_seq_list = (list(set([a for b in self.no_bead_families_df.Sequence_no_color.tolist() for a in b])))
                     
            # Create dictionary for each sequence and its abundance
            self.abundance_dict = {}
            
            for seq in self.no_duplicate_seq_list:
                
                # Create counter for abundance of each sequence
                self.abundance_counter = 0
                
                # Set up partial match criterion with N allowed substitutions
                self.match_pattern = ('(%s){s<=%s}' % (seq, 2))
                    
                for j in self.sequence_df.itertuples():
                    self.comparison_seq = j.Sequence
                        
                    # Use regex.search to find matches with the allowed # of substitutions
                    self.find_match = regex.search(self.match_pattern, self.comparison_seq)
                    
                    if self.find_match:
                        self.abundance_counter += 1
       	
                # Append seq and self.abundance_counter to self.abundance_dict
                self.abundance_dict.update({seq:self.abundance_counter})
                
            # Append dictionary to dataframe
            self.abundance_df = pd.DataFrame(self.abundance_dict.items(), columns=['Sequence', 'Abundance'])  
            
            # Order from highest to lowest abundance
            self.abundance_df = self.abundance_df.sort_values(by = 'Abundance', ascending = False) 
            
            # Convert abundance column to string (can't print to html without this)
            self.abundance_df['Abundance'] = self.abundance_df['Abundance'].astype(str) 
            # print ('abundance df sorted:')
            # print(self.abundance_df)  
            
  
    def bead_check(self):
        
        # Create a final list of yes/no bead matches for each sequence in top ten
        self.bead_matches_found_list = []
        
        # Check every sequence in self.abundance_df against every sequence in self.bead_sequences_df for close matches
        for row in self.abundance_df.itertuples():
            
            # Define search conditions (allows for 4 subs)
            self.match_params = ('(%s){s<=%s}' % (row.Sequence, 4))
        
            # Create a temp list of yes/no for match or no match
            self.bead_match_list = []
        
            # Compare row.Sequence to every sequence in bead database according to conditions in match_oligo_criterion
            for row in self.bead_sequence_df.itertuples():
        
                # Search for close matches
                self.bead_match = regex.search(self.match_params, row.Sequence)
        
                if self.bead_match:
                    # Add yes to bead_match_list
                    self.bead_match_list.append('yes')
        
                else:
                    # Add no to bead_match_list
                    self.bead_match_list.append('no')
        
            # If bead_match_list contains a single yes add yes to new bead match column in self.abundance_df
            # Else add no
        
            if 'yes' in self.bead_match_list:
                self.bead_matches_found_list.append('Yes')
        
            else:
                self.bead_matches_found_list.append('No')
    
            
        # Append matches_found list to sequence_df
        # self.abundance_df.insert(loc = 3, column = 'Bead Match (y/n)', value = self.bead_matches_found_list)
        self.abundance_df['Bead_Match'] = self.bead_matches_found_list
        # print(self.abundance_df)          
                           
    # Function to append top sequences and their abundance to HTML file  
    def append_top_seqs(self):
        
        # Append to html
        with open(self.results_path_name, "a") as file:
            file.write("Sequences that appeared in families with no relation to bead sequences:")
            file.write(self.crlb)  
                   
            for row in self.abundance_df.itertuples():
                # Write sequence to html
                self.seq = row.Sequence
                file.write(self.seq)
                file.write(self.cr)
                    
                # Write Abundance to html
                self.abundance = 'Abundance: ' + row.Abundance 
                file.write(self.abundance)
                file.write(self.cr)    
               
                # Write bead matches to html
                self.bead_matches = 'Matched to bead sequence (Y/N): ' + row.Bead_Match
                file.write(self.bead_matches)
                file.write(self.crlb)
            file.write(self.crlb)
    
    # Export scores to csv file        
    def export_csv(self):
       
        # Get first 6 rows of self.abundance_df
        self.export_df = self.abundance_df.head(6)
     
        with open(self.most_abundant_from_py_aligment, "a") as f:
            wr = csv.writer(f, dialect = 'excel')
            wr.writerow([])
            
            # Convert each row to a list
            for row in self.export_df.itertuples():
                self.row_list = [self.new_name, row.Sequence, row.Abundance]
                wr.writerow(self.row_list)
                
            # Insert a space between data for each sequencing file    
            wr.writerow([])
        
                
    def export_motifs(self):
        
        # Export a csv file containing all the motifs from the current file
        # Use the motifs from self.no_bead_families_df to exclude bead motifs
        # Create new df from motif column of self.no_bead_families_df
        self.motif_export_df = self.no_bead_families_df[['Motif']].copy()
        # Add a column containing the current filename to self.no_bead_families_df
        # This will throw a warning if done on the original df self.no_bead_families_df
        self.motif_export_df['Filename'] = self.new_name
        
        # Create a path and filename for motifs
        # Create a path and filename for results
        self.motif_path_name = os.path.join(self.results_path, self.new_name + "_motifs" + ".csv")
        
        # Export filename and motif columns to csv
        with open(self.motif_path_name, "a") as f:
            wr = csv.writer(f, dialect = 'excel')
            # wr.writerow([])
            
            # Convert each row to a list
            for row in self.motif_export_df.itertuples():
                self.row_list = [row.Filename, row.Motif]
                wr.writerow(self.row_list)
            
        # Export motifs to the database of motifs for all selections    
        with open(self.motif_db_filename, "a") as f:
            wr = csv.writer(f, dialect = 'excel')
            wr.writerow([])
            
            # Convert each row to a list
            for row in self.motif_export_df.itertuples():
                self.row_list = [row.Filename, row.Motif]
                wr.writerow(self.row_list)
                
            # Insert a space between data for each sequencing file    
            wr.writerow([])
        
            
        
               
            
            
            
                                                      
                                                            
# Call class
obj1 = Motif()

# Call functions
importer = obj1.import_files()

# The remainder of these calls have been moved to the import_files function so for loop could be run
# Use one of the next two to import a bead database
# For fasta format
# fasta_bead_import = obj1.import_beaddb()
# For csv format
# csv_bead_import = obj1.bead_csv()

# merge = obj1.join_df()
# motifs = obj1.motifs()
# families = obj1.families()
# html = obj1.create_html()
# no_bead_families = obj1.no_bead_families()
# abundance = obj1.abundance()
# bead_check = obj1.bead_check()
# top_seqs = obj1.append_top_seqs()
# export_csv = obj1.export_csv()


    
    
    
    
    
    
    
    
    
    
    
    
    