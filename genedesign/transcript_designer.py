from genedesign.seq_utils.translate import Translate
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.seq_utils.hairpin_counter import hairpin_counter
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.codon_checker import CodonChecker

import random
from itertools import product
import re

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence.
    Attempts to minimize hairpins, rare codons, forbidden sequences, and internal promoters.
    Attempts to maximmize CAI and codon diversity.
    Chooses an RBS to best satisfy above criteria.
    """

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        # Codon table with highest CAI codon for each amino acid (for E. coli)
        self.aminoAcidToCodon = {
            'A': "GCG", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTC",
            'G': "GGT", 'H': "CAC", 'I': "ATC", 'K': "AAA", 'L': "CTG",
            'M': "ATG", 'N': "AAC", 'P': "CCG", 'Q': "CAG", 'R': "CGT",
            'S': "TCT", 'T': "ACC", 'V': "GTT", 'W': "TGG", 'Y': "TAC"
        } 

        self.codon_checker = CodonChecker()
        self.codon_checker.initiate()
        for freq_list in self.codon_checker.aa_map.values():
            freq_list.sort(reverse=True, key=(lambda x: x[1]))
        
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.forbidden_checker.initiate()

        self.promoter_checker = PromoterChecker()
        self.promoter_checker.initiate()

        self.translator = Translate()
        self.translator.initiate()

    def count_forbidden(self, construct: str) -> int:
        """
        Counts number of forbidden sequences in a DNA construct.

        Parameters:
            construct (str): DNA sequence (all characters in set "ATCG"). 

        Returns:
            int: number of forbidden sequences in DNA construct,
        """
        score = 0
        for seq in self.forbidden_checker.forbidden:
            score += len(re.findall(seq, construct))
        return score

    def pick_rscu(self, aa: str) -> str:
        """
        Picks DNA codon for an amino acid based on relative synonymous codon usage.
        Also attempts to improve codon diversity by incorporating some lower RSCU codons.

        Paramters:
            aa (str): Standard IUPAC code representing one amino acid (e.g. "M", "K").

        Returns:
            str: String representing chosen codon.
        """
        freq_list = self.codon_checker.aa_map[aa]
        random_float = random.random()
        codon_choice = freq_list[0][0]
        for codon in freq_list:
            random_float -= codon[1]
            if random_float < 0.0:
                if codon[0] in self.codon_checker.rare_codons:
                    codon_choice = freq_list[0][0]
                else:
                    codon_choice = codon[0]
        # Add diversity
        random_float = random.random()
        if random_float < 0.5 and len(freq_list) > 1 and codon_choice == freq_list[0][0]:
            index = random.randint(1, len(freq_list) - 1)
            if freq_list[index][0] in self.codon_checker.rare_codons:
                codon_choice = freq_list[0][0]

        return codon_choice

    def pick_random(self, aa: str) -> str: 
        """
        Picks a DNA codon for an amino acid randomly, avoiding rare codons.

        Parameters:
            aa (str): Standard IUPAC code representing one amino acid (e.g. "M", "K").

        Returns:
            str: String representing chosen codon.
        """
        freq_list = self.codon_checker.aa_map[aa]
        index = random.randint(0, len(freq_list) - 1)
        codon = freq_list[index][0]
        if codon in self.codon_checker.rare_codons:
            codon = freq_list[0][0]
        return codon 

    def enumerate_codons(self, peptide: str) -> list[list[str]]:
        """
        Enumerates all possible combinations of codons for a peptide sequence.

        Parameters:
            peptide (str): Peptide sequence with amino acids represented in standard IUPAC codes.

        Returns:
            list[list[str]]: List of all possible codon combinations for this peptide.
                             Each codon combination is a list, e.g. ["TTC", "TTG"].
        """
        aa_map = self.codon_checker.aa_map
        codon_lists = [aa_map[aa] for aa in peptide]

        combinations = []
        for codon_tuple in product(*codon_lists):
            combinations.append([codon[0] for codon in codon_tuple])

        return combinations

    def score_construct(self, construct: str) -> int:
        """
        Score a DNA construct based on presence of hairpins, forbidden sequences, and internal promoters.
        Lower score indicates a "better" construct.

        Parameters:
            construct (str): DNA sequence (all characters in set "ATCG"). 
        
        Returns:
            int: Calculated score for the construct.   
        """
        score = 0

        # Check for hairpins, add 1 to score for each hairpin
        hairpin_result = hairpin_counter(construct)
        if hairpin_result:
            score += hairpin_result[0]
        
        # Check for forbidden sequences, add 1 to score for each sequence
        score += self.count_forbidden(construct)

        # Add 15 to score if an internal promoter is found within sequence
        passed_promoter, _ = self.promoter_checker.run(construct)
        if not passed_promoter:
            score += 15

        return score

    def score_codons(self, codons: list[str]) -> int:
        """
        Score a list of codons based on codon diversity, rare codon count, and codon adaptiveness index.
        Lower score indicates a "better" construct.

        Parameters:
            codons (list[str]): List of three base-pair DNA codons in the construct. 
  
        Returns:
            int: Calculated score for the list of codons. 
        """
        score = 0
        codons_above_board, codon_diversity, rare_codon_count, cai_value = self.codon_checker.run(codons)
        if rare_codon_count > 3:
            score += rare_codon_count
        if codon_diversity < 0.5:
            score += 6
        if cai_value < 0.2:
            score += 6
        return score, [codon_diversity, rare_codon_count, cai_value]

    def eliminate_forbidden(self, codons: list[str], utr: str, seq: str, peptide: str) -> list[str]:
        """
        Attempts to eliminate all instances of a forbidden sequence from a construct.
        Uses a "sliding window" style approach by enumerating all possible codon combinations.
        Modifies sequence only if overall construct score also decreases.

        Parameters:
            codons (list[str]): List of three base-pair DNA codons in the construct. 
            utr (str): The untranslated region (UTR) sequence of the construct.
            seq (str): The forbidden sequence (e.g. "ACTAGT" for SpeI).
            peptide (str): The peptide sequence corresponding to the codons.

        Returns:
            list[str]: Modified list of codons with forbidden sequence eliminated (where possible).
        """
        start = 0
        construct = utr + ''.join(codons)
        while True:
            start = construct.find(seq, start)
            
            # Break out of loop if no more instances of the forbidden sequence
            if start == -1:
                break

            # Find codons containing start and end of forbidden sequence
            utr_length = len(utr)
            start_codon = (start - utr_length) // 3
            start_codon = max(start_codon, 0)
            end_codon = (start - utr_length + len(seq)) // 3
            end_codon = min(end_codon, len(codons) - 1)

            # Calculate initial construct score
            initial_forbidden = self.count_forbidden(construct)
            initial_score = self.score_construct(construct)
            initial_codon_score, _ = self.score_codons(codons)
            initial_score += initial_codon_score

            for codon_set in self.enumerate_codons(peptide[start_codon:end_codon + 1]):
                test_codons = codons[:]
                test_codons[start_codon:start_codon + len(codon_set)] = codon_set
                test_construct = utr + ''.join(test_codons)
                
                # Check if the new construct eliminates the forbidden sequence
                test_forbidden = self.count_forbidden(test_construct)
                if test_forbidden < initial_forbidden:
                    # Check overall scores to ensure no increase
                    test_score = self.score_construct(test_construct)
                    test_codon_score, _ = self.score_codons(test_codons)
                    test_score += test_codon_score
                    if test_score < initial_score:
                        codons = test_codons
            
            # Increment by 1 to detect overlapping forbidden sequences
            start += 1

        return codons

    def eliminate_hairpin(self, codons: list[str], utr: str, hairpin: str, peptide: str) -> list[str]:
        """
        Attempts to eliminate a hairpin by modifying codons in its stem regions.
        Uses a "sliding window" style approach by enumerating all possible codon combinations.
        Modifies sequence only if overall construct score also decreases.
    
        Parameters:
            codons (list[str]): List of three base-pair DNA codons in the construct. 
            utr (str): The untranslated region (UTR) sequence of the construct.
            hairpin (str): The hairpin sequence in the format "ATC(CGTCG)GAT".
            peptide (str): The peptide sequence corresponding to the codons.
    
        Returns:
            list[str]: Modified list of codons with hairpin eliminated (if possible).
        """
        # Assemble the full construct
        construct = utr + ''.join(codons)
        initial_hairpins = hairpin_counter(construct)
        
        # Return if no hairpins, otherwise count hairpins in initial construct
        if initial_hairpins == None:
            return codons
        else:
            initial_hairpins = initial_hairpins[0]

        # Extract the stem regions from the hairpin
        stem_start = hairpin.split("(")[0]
        stem_end = hairpin.split(")")[1]

        # Locate stem regions in the construct
        clean_hairpin = hairpin.replace("(", "").replace(")", "")
        stem1_loc = construct.find(clean_hairpin)
        stem2_loc = stem1_loc + len(clean_hairpin) - len(stem_end)

        # Ensure both stems are located
        if stem1_loc == -1 or stem2_loc == -1:
            return codons

        # Map stem locations to codons
        utr_length = len(utr)
        stem1_start_codon = (stem1_loc - utr_length) // 3
        stem2_start_codon = (stem2_loc - utr_length) // 3
        stem1_start_codon = max(stem1_start_codon, 0)
        stem2_start_codon = max(stem2_start_codon, 0)

        stem1_end_codon = (stem1_loc - utr_length + len(stem_start)) // 3
        stem2_end_codon = (stem2_loc - utr_length + len(stem_end)) // 3
        stem1_end_codon = min(stem1_end_codon, len(codons) - 1)
        stem2_end_codon = min(stem2_end_codon, len(codons) - 1)
        stem1_end_codon = max(stem1_end_codon, 0)
        stem2_end_codon = max(stem2_end_codon, 0)

        # Initial construct scores
        initial_score = self.score_construct(construct)
        initial_codon_score, _ = self.score_codons(codons)
        initial_score += initial_codon_score

        for codon1 in self.enumerate_codons(peptide[stem1_start_codon:stem1_end_codon + 1]):
            for codon2 in self.enumerate_codons(peptide[stem2_start_codon:stem2_end_codon + 1]):
                # Test the replacement codons
                test_codons = codons[:]
                test_codons[stem1_start_codon:stem1_start_codon + len(codon1)] = codon1
                test_codons[stem2_start_codon:stem2_start_codon + len(codon2)] = codon2
                test_construct = utr + ''.join(test_codons)
                
                # Check if the new construct eliminates the hairpin
                test_hairpins = hairpin_counter(test_construct)
                if test_hairpins == None or test_hairpins[0] < initial_hairpins: 
                    # Check overall scores to ensure no increase
                    test_score = self.score_construct(test_construct)
                    test_codon_score, _ = self.score_codons(test_codons)
                    test_score += test_codon_score
                    if test_score < initial_score:
                        return test_codons
       
        return codons

    def get_hairpin_locus(self, construct: str) -> tuple[list[str], list[int]]:
        """
        Given a DNA construct, locate all  hairpins and return hairpin sequences and locations.

        Parameters:
            construct (str): DNA sequence (all characters in set "ATCG"). 
        
        Returns:
            tuple[list[str], list[int]]: List of hairpins in the format "ATC(CGTCG)GAT".
                                         List of start locations for each hairpin.
        """
        hairpin_result = hairpin_counter(construct)
        if hairpin_result == None:
            return [], []
        hairpin_locations = []
        hairpin_string = hairpin_result[1]
        if hairpin_string == None:
            return [], []
        hairpin_list = hairpin_string.split("\n")
        hairpin_list.pop()
        for i in range(len(hairpin_list)):
            hairpin = hairpin_list[i].split(": ")[1]
            hairpin_list[i] = hairpin 
            hairpin = hairpin.replace("(", "").replace(")", "")
            locus = construct.find(hairpin)
            hairpin_locations.append(locus)
        return hairpin_list, hairpin_locations    

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """

        # Original implementation of run (preserved for testing)
        """
        # Translate peptide to codons
        codons = [self.aminoAcidToCodon[aa] for aa in peptide]

        # Append the stop codon (TAA in this case)
        codons.append("TAA")

        # Build the CDS from the codons
        cds = ''.join(codons)

        # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, codons)
        #"""
        

        filtered_rbs_options = []
        best_rbs = None
        best_score = 9999
        for rbs in self.rbsChooser.rbsOptions:
            if rbs in ignores:
                continue
            utr = rbs.utr.upper()
            score = self.score_construct(utr)
            if score < best_score:
                best_rbs = rbs
                best_score = score
            if score == 0:
                filtered_rbs_options.append(rbs)

        codons = []

        # If no RBS determined above, default to RBSChooser
        if best_rbs == None:
            best_rbs = self.rbsChooser.run(cds, ignores)

        best_score = 9999
        best_scores_list = [(999, []) for _ in range(2)]
        
        for i in range(300):
            candidate_codons = [self.pick_random(aa) for aa in peptide]
            candidate_codons.append(random.choice(["TAA", "TAG", "TGA"]))
            if filtered_rbs_options:
                rbs = random.choice(filtered_rbs_options)
            else:
                rbs = best_rbs
            utr = rbs.utr.upper()
            candidate_seq = utr + ''.join(candidate_codons)
            candidate_score = self.score_construct(candidate_seq)
            codon_score, _ = self.score_codons(candidate_codons)
            candidate_score += codon_score
            if candidate_score == 0:
                codons = candidate_codons
                best_score = candidate_score
                break
            if candidate_score < best_scores_list[-1][0]:
                best_scores_list[-1] = (candidate_score, candidate_codons, rbs)
                best_scores_list.sort(key=(lambda x: x[0]))

        best_score = best_scores_list[0][0]
        codons = best_scores_list[0][1]
        best_rbs = best_scores_list[0][2]

        # Use variation on sliding window approach to attempt to optimize targeted regions
        for i in range(len(best_scores_list)):
            candidate = best_scores_list[i]
            candidate_score = candidate[0]
            candidate_codons = candidate[1]
            candidate_rbs = candidate[2]

            # Attempt to eliminate hairpins
            for j in range(2):
                candidate_construct = candidate_rbs.utr.upper() + ''.join(candidate_codons)
                hairpin_list, hairpin_locs = self.get_hairpin_locus(candidate_construct)
                for k in range(len(hairpin_list)):
                    candidate_codons = self.eliminate_hairpin(candidate_codons, 
                                                             candidate_rbs.utr.upper(), 
                                                              hairpin_list[k], 
                                                              peptide)

            # Attempt to eliminate forbidden sequences
            for seq in self.forbidden_checker.forbidden:
                candidate_codons = self.eliminate_forbidden(candidate_codons, candidate_rbs.utr.upper(), seq, peptide)

            new_score, _ = self.score_codons(candidate_codons)
            new_score += self.score_construct(candidate_rbs.utr.upper() + ''.join(candidate_codons))
            best_scores_list[i] = (new_score, candidate_codons, candidate_rbs)
                        
        best_scores_list.sort(key=(lambda x: x[0]))
        best_score = best_scores_list[0][0]
        codons = best_scores_list[0][1]
        best_rbs = best_scores_list[0][2]

        return Transcript(best_rbs, peptide, codons)

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)
