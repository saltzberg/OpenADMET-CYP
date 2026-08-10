# Human CYP target sequences

The FASTA files in this folder contain the full-length canonical human sequences from reviewed UniProtKB/Swiss-Prot entries. These include the native N-terminal membrane-anchor regions; no construct trimming, mutation, tag, heme record, or redox-partner sequence has been introduced.

| Target | UniProtKB | Entry | Length | Sequence version | UniProt MD5 |
|---|---|---|---:|---:|---|
| CYP3A4 | [P08684](https://www.uniprot.org/uniprotkb/P08684/entry) | CP3A4_HUMAN | 503 aa | 4 | `0D1C7296209804664D2DA95C23EEA089` |
| CYP2C9 | [P11712](https://www.uniprot.org/uniprotkb/P11712/entry) | CP2C9_HUMAN | 490 aa | 3 | `189EDC2D0AD3DFF2B6D26A9E2C32ED86` |
| CYP2D6 | [P10635](https://www.uniprot.org/uniprotkb/P10635/entry) | CP2D6_HUMAN | 497 aa | 2 | `C71FFC3C831A845301011DF0B1EC5638` |
| CYP1A2 | [P05177](https://www.uniprot.org/uniprotkb/P05177/entry) | CP1A2_HUMAN | 516 aa | 4 | `BCEC75A2BE2D5C4654F8D58B240766C2` |

Organism for every entry: *Homo sapiens* (NCBI taxonomy 9606).

Files:

- `CYP3A4_HUMAN_P08684.fasta`
- `CYP2C9_HUMAN_P11712.fasta`
- `CYP2D6_HUMAN_P10635.fasta`
- `CYP1A2_HUMAN_P05177.fasta`
- `cyp_panel_human_canonical.fasta` (the same four records concatenated)

Source endpoint pattern: `https://rest.uniprot.org/uniprotkb/<accession>.fasta`.
