from genetic_variant_interpreter import analyze_SNP, analyze_RNA, analyze_protein, classify_mutations


def test_analyze_SNP():

    sequences = {
        "reference": "ATGCCG",
        "experimental": "ATGTCG"
    }

    result = analyze_SNP(sequences)

    expected = [[4, "C", "T"]]
    wrong_result = [[5, "C", "T"]]


    assert result == expected
    assert result != wrong_result


def test_analyze_RNA():

    sequences = {
        "reference": "ATGCCGTTA",
        "experimental": "ATGTCGTAA"
    }

    ref_rna, exp_rna = analyze_RNA(sequences)

    assert ref_rna == "AUGCCGUUA"
    assert exp_rna == "AUGUCGUAA"

    assert "T" not in ref_rna
    assert "T" not in exp_rna


def test_analyze_protein():

    sequences = {
        "reference": "ATGGCTAAATGA",
        "experimental": "ATGGCGAAGTGA"
    }

    ref_protein, exp_protein = analyze_protein(sequences, "1")

    assert ref_protein == "MAK"
    assert exp_protein == "MAK"

    assert "STOP" not in ref_protein
    assert "STOP" not in exp_protein


def test_classify_mutations():
    sequences = {
        "reference": "GCTTTGTGG",
        "experimental": "GCCATGTGA"
    }

    SNPs = analyze_SNP(sequences)

    mutations = classify_mutations(sequences, "1", SNPs)

    mutation_types = [mutation[4] for mutation in mutations]

    assert "Synonymous" in mutation_types
    assert "Missense" in mutation_types
    assert "Nonsense" in mutation_types
