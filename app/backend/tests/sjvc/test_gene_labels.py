"""Gene-level matrix row labels come in three shapes the app must all handle:

- a bare gene symbol (older / hand-built matrices), e.g. "TP53"
- "<gene_name>:<gene_id>" (prepTCGAdata::count_novel_sjs() output — every row
  is exactly one gene, but the id disambiguates symbols that repeat across
  loci, e.g. the Rfam small RNAs; see annotate_sj())
- a comma-joined multi-gene row from an offline-condensed matrix, e.g.
  "CFH,CFHR1,CFHR3"
"""
import numpy as np

from sjvc.services.features import FeatureMatrix
from sjvc.services.geneset import match_matrix_labels
from sjvc.services.junctions import gene_name_of_label
from sjvc.services.mad import top_features_by_mad
from sjvc.services.rds import Matrix


def test_gene_name_of_label_shapes():
    assert gene_name_of_label("TP53") == "TP53"
    assert gene_name_of_label("TP53:ENSG00000141510.19") == "TP53"
    assert gene_name_of_label("TP53:ENSG00000141510") == "TP53"          # unversioned id
    assert gene_name_of_label("RF00019:ENSG00000200344.1") == "RF00019"
    # a comma row is left alone — callers split it themselves
    assert gene_name_of_label("CFH,CFHR1,CFHR3") == "CFH,CFHR1,CFHR3"
    # not a real gene_id suffix -> treated as an opaque (already-bare) label
    assert gene_name_of_label("SOME:THING") == "SOME:THING"


def test_match_matrix_labels_name_colon_id_rows():
    labels = ["TP53:ENSG00000141510.19", "KRAS:ENSG00000133703.11", "RF00019:ENSG00000200344.1"]
    matched, unmatched, warnings = match_matrix_labels(["TP53", "kras", "NOPE"], labels)
    assert matched == ["TP53:ENSG00000141510.19", "KRAS:ENSG00000133703.11"]
    assert unmatched == ["NOPE"]
    assert warnings == []  # each row is unambiguously one gene -> no "merged" warning


def test_match_matrix_labels_duplicate_symbol_different_ids_both_addressable():
    """Two distinct gene_ids sharing a symbol (e.g. RF00019 copies) are two
    separate rows; searching the symbol should pull in both."""
    labels = ["RF00019:ENSG00000200344.1", "RF00019:ENSG00000252254.1", "TP53:ENSG00000141510.19"]
    matched, unmatched, warnings = match_matrix_labels(["RF00019"], labels)
    assert set(matched) == {"RF00019:ENSG00000200344.1", "RF00019:ENSG00000252254.1"}


def test_match_matrix_labels_still_handles_bare_and_comma_rows():
    labels = ["TP53", "CFH,CFHR1,CFHR3"]
    matched, unmatched, warnings = match_matrix_labels(["TP53", "CFHR1"], labels)
    assert matched == ["TP53", "CFH,CFHR1,CFHR3"]
    assert unmatched == []
    assert "CFHR1" in warnings[0]  # merged-row fallback still warns


def test_mad_protein_coding_filter_sees_through_name_colon_id():
    labels = ["TP53:ENSG00000141510.19", "KRAS:ENSG00000133703.11", "NONCODING:ENSG00000000001.1"]
    m = Matrix(samples=["s1", "s2", "s3"], features=labels,
               values=np.array([[1.0, 2.0, 3.0], [4.0, 1.0, 0.0], [9.0, 9.0, 9.0]]))
    fm = top_features_by_mad(m, top_n=10, restrict_to={"tp53", "kras"})
    assert set(fm.feature_ids) == {"TP53:ENSG00000141510.19", "KRAS:ENSG00000133703.11"}


def test_heatmap_row_labels_trim_name_colon_id_but_not_junction_or_pathway_rows():
    """The heatmap endpoint composes FeatureMatrix.label() with
    gene_name_of_label() (see api/routes.py's /heatmap handler) so a
    count_novel_sjs()-style "<name>:<gene_id>" row displays as just the gene
    name, while junction (chr:start-end:strand) and pathway-style rows
    (whatever follows their colon isn't a real Ensembl gene id) pass through
    unchanged."""
    gene_fm = FeatureMatrix(
        kind="gene",
        feature_ids=["TP53:ENSG00000141510.19", "RF00019:ENSG00000200344.1"],
        samples=["s1"], values=np.zeros((2, 1)), n_junctions=0,
    )
    assert [gene_name_of_label(gene_fm.label(f)) for f in gene_fm.feature_ids] == ["TP53", "RF00019"]

    junction_fm = FeatureMatrix(
        kind="junction", feature_ids=["chr1:100-200:+"], samples=["s1"],
        values=np.zeros((1, 1)), n_junctions=1,
    )
    assert gene_name_of_label(junction_fm.label("chr1:100-200:+")) == "chr1:100-200:+"

    pathway_fm = FeatureMatrix(
        kind="gene", feature_ids=["xCell:aDC%HPCA%1.txt"], samples=["s1"],
        values=np.zeros((1, 1)), n_junctions=0,
    )
    assert gene_name_of_label(pathway_fm.label("xCell:aDC%HPCA%1.txt")) == "xCell:aDC%HPCA%1.txt"
