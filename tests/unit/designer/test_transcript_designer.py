"""
test_transcript_designer.py:

Purpose: Tests the TranscriptDesigner class to ensure mRNA sequences are correctly generated from protein sequences.
Coverage:
** Test the reverse-translation of protein sequences into optimized DNA.
** Ensure proper RBS assignment for each mRNA.
** Validate handling of codon optimization and RNA folding requirements.
"""

import pytest
from genedesign.transcript_designer import TranscriptDesigner
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.seq_utils.hairpin_counter import hairpin_counter

@pytest.fixture
def transcript_designer():
    designer = TranscriptDesigner()
    designer.initiate()
    return designer

def test_enumerate_codons(transcript_designer):
    peptide = "FF"
    correct = {"TTTTTT", "TTTTTC", "TTCTTT", "TTCTTC"}

    enumerations = transcript_designer.enumerate_codons(peptide)
    answers = set()
    for possibility in enumerations:
        answers.add(''.join(possibility))
    assert correct == set(answers)

def test_rbs_ignores(transcript_designer):
    rbs_options = transcript_designer.rbsChooser.rbsOptions
    ignore = rbs_options[0]
    peptide = "MFFFKTAGMSDIDPFACEEEFLRFLRRRRY"
    transcript = transcript_designer.run(peptide, {ignore})
    assert transcript.rbs != ignore

def test_translation_accuracy(transcript_designer):
    peptide = "MFFFKTAGMSDIDPFACEEEFLRFLRRRRY"
    transcript = transcript_designer.run(peptide, {})
    cds = ''.join(transcript.codons)
    translated = transcript_designer.translator.run(cds)
    assert translated == peptide

def test_eliminate_hairpin(transcript_designer):
    hairpin = ["TAT", "CCC", "CCC", "ATA"] 
    filler = ["CCC" for _ in range(15)]
    codons = filler + hairpin + filler
    peptide = ("P" * 15) + "YPPI" + ("P" * 15)
    utr = transcript_designer.rbsChooser.rbsOptions[0].utr.upper()
    new_codons = transcript_designer.eliminate_hairpin(codons, utr, "TAT(CCCCCC)ATA", peptide)
    result = hairpin_counter(''.join(new_codons))
    assert result == (0, None)

def test_eliminate_forbidden(transcript_designer):
    ecoRI = "GAATTC"
    filler = ["TAT", "CCC", "ATA"] * 3
    filler_peptide = "YPI" * 3
    codons = filler + ["GAA", "TTC"] + filler
    peptide = filler_peptide + "EF" + filler_peptide
    utr = transcript_designer.rbsChooser.rbsOptions[0].utr.upper()
    new_codons = transcript_designer.eliminate_forbidden(codons, utr, ecoRI, peptide)
    result = transcript_designer.forbidden_checker.run(''.join(new_codons))
    assert result == (True, None)
