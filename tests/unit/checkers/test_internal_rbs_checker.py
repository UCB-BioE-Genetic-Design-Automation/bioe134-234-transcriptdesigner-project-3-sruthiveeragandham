import pytest
from genedesign.checkers.internal_rbs_checker import internal_rbs_checker

def test_no_shine_delgarno():
    """Test with a DNA sequence that has no Shine-Dalgarno sequence."""
    dna_seq = "ATGCATCATCGTACGTAGCTAGCATCGTATG"
    status, results = internal_rbs_checker(dna_seq)
    assert status is True
    assert results is None

def test_no_start_codon():
    """Test with a Shine-Dalgarno sequence but no valid start codon in spacing range."""
    dna_seq = "TTTAAGGAGGTTTTTTTACATTGGG"
    status, results = internal_rbs_checker(dna_seq)
    assert status is True
    assert results is None

def test_single_rbs():
    """Test with a DNA sequence that contains just one Shine-Dalgarno sequence and ATG."""
    dna_seq = "TTACCGCTAGGACACGCCATG"
    status, results = internal_rbs_checker(dna_seq)
    assert status is False
    assert len(results) == 1  # Only one RBS
    assert results[0]["Shine-Dalgarno Sequence"] == "GGA"
    assert results[0]["Start Codon Position"] - results[0]["SD End"] == 6

def test_multiple_rbs():
    """Test with a DNA sequence that contains multiple valid ribosome binding sites."""
    dna_seq = "ATCAAGGAGGCCCCCCCCCCCCATG"
    status, results = internal_rbs_checker(dna_seq)
    assert status is False
    assert len(results) == 2
    sequences = set([seq["Shine-Dalgarno Sequence"] for seq in results])
    print(sequences)
    assert "AAGG" in sequences
    assert "AGGA" in sequences


def test_edge_case_spacing():
    """Test with valid Shine-Dalgarno and an ATG start codon just within the max spacing range."""
    dna_seq = "TTT" + "AAGG" + ("C" * 23) + "ATG"
    status, results = internal_rbs_checker(dna_seq)
    assert status is False
    assert len(results) == 1
    assert results[0]["Shine-Dalgarno Sequence"] == "AAGG"
    assert results[0]["Start Codon Position"] - results[0]["SD End"] == 23
    """Test with valid Shine-Dalgarno and an ATG start codon just within the min spacing range."""
    dna_seq = "TTT" + "AAGGG" + "CC" + "ATG"
    status, results = internal_rbs_checker(dna_seq)
    assert status is False
    assert len(results) == 1
    assert results[0]["Shine-Dalgarno Sequence"] == "AAGGG"
    assert results[0]["Start Codon Position"] - results[0]["SD End"] == 2

def test_multiple_start_codons():
    """Test with multiple ATG start codons in the spacing range."""
    dna_seq = "AAGG" + ("C" * 10) + "ATG" + "CC" + "ATG"
    status, results = internal_rbs_checker(dna_seq)
    assert status is False
    assert len(results) == 2
    assert results[0]["Shine-Dalgarno Sequence"] == "AAGG"
    assert results[1]["Shine-Dalgarno Sequence"] == "AAGG"
    start_codons = set([seq["Start Codon Position"] for seq in results])
    assert 14 in start_codons # AAGGG + first start codon
    assert 19 in start_codons # AAGGG + second start codon
