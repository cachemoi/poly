package codon_test

import (
	"fmt"
	"os"

	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/TimothyStiles/poly/synthesis/codon"
)

func ExampleTranslationTable_Translate() {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"
	testTranslation, _ := codon.NewTranslationTable(11).Translate(gfpDnaSequence) // need to specify which codons map to which amino acids per NCBI table

	fmt.Println(gfpTranslation == testTranslation)
	// output: true
}

func ExampleTranslationTable_OptimizeSequence() {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	sequence, _ := genbank.Read("../../data/puc19.gbk")
	codonTable := codon.NewTranslationTable(11)
	stats, _ := codonTable.UpdateWeightsWithSequence(sequence)

	// Here, we double check if the number of genes is equal to the number of stop codons
	stopCodonCount := 0
	for _, aa := range codonTable.AminoAcids {
		if aa.Letter == "*" {
			for _, codon := range aa.Codons {
				stopCodonCount = stopCodonCount + codon.Weight
			}
		}
	}
	if stopCodonCount != stats.GeneCount {
		fmt.Println("Stop codons don't equal number of genes!")
	}

	optimizedSequence, _ := codonTable.OptimizeSequence(gfpTranslation)
	optimizedSequenceTranslation, _ := codonTable.Translate(optimizedSequence)

	fmt.Println(optimizedSequenceTranslation == gfpTranslation)
	// output: true
}

func ExampleReadCodonJSON() {
	codontable := codon.ReadCodonJSON("../../data/bsub_codon_test.json")

	fmt.Println(codontable.GetWeightedAminoAcids()[0].Codons[0].Weight)
	//output: 28327
}

func ExampleParseCodonJSON() {
	file, _ := os.ReadFile("../../data/bsub_codon_test.json")
	codontable := codon.ParseCodonJSON(file)

	fmt.Println(codontable.GetWeightedAminoAcids()[0].Codons[0].Weight)
	//output: 28327
}

func ExampleWriteCodonJSON() {
	codontable := codon.ReadCodonJSON("../../data/bsub_codon_test.json")
	codon.WriteCodonJSON(codontable, "../../data/codon_test.json")
	testCodonTable := codon.ReadCodonJSON("../../data/codon_test.json")

	// cleaning up test data
	os.Remove("../../data/codon_test.json")

	fmt.Println(testCodonTable.GetWeightedAminoAcids()[0].Codons[0].Weight)
	//output: 28327
}

func ExampleCompromiseCodonTable() {
	sequence, _ := genbank.Read("../../data/puc19.gbk")

	// weight our codon optimization table using the regions we collected from the genbank file above
	optimizationTable := codon.NewTranslationTable(11)
	optimizationTable.UpdateWeightsWithSequence(sequence)

	sequence2, _ := genbank.Read("../../data/phix174.gb")
	optimizationTable2 := codon.NewTranslationTable(11)
	optimizationTable2.UpdateWeightsWithSequence(sequence2)

	finalTable, _ := codon.CompromiseCodonTable(optimizationTable, optimizationTable2, 0.1)
	for _, aa := range finalTable.GetWeightedAminoAcids() {
		for _, codon := range aa.Codons {
			if codon.Triplet == "TAA" {
				fmt.Println(codon.Weight)
			}
		}
	}
	//output: 2727
}

func ExampleAddCodonTable() {
	sequence, _ := genbank.Read("../../data/puc19.gbk")

	// weight our codon optimization table using the regions we collected from the genbank file above
	optimizationTable := codon.NewTranslationTable(11)
	optimizationTable.UpdateWeightsWithSequence(sequence)

	sequence2, _ := genbank.Read("../../data/phix174.gb")
	optimizationTable2 := codon.NewTranslationTable(11)
	optimizationTable2.UpdateWeightsWithSequence(sequence2)

	finalTable := codon.AddCodonTable(optimizationTable, optimizationTable2)
	for _, aa := range finalTable.GetWeightedAminoAcids() {
		for _, codon := range aa.Codons {
			if codon.Triplet == "GGC" {
				fmt.Println(codon.Weight)
			}
		}
	}
	//output: 90
}
