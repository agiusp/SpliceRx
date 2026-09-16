"""services.junction_type — the junction-level heatmap's alternate
"<gene name>:<novel-splicing-event type>" row labels."""
import pytest

from sjvc.services.junction_metadata import JunctionMetadataError
from sjvc.services.junction_type import load_junction_type_index


def _csv(tmp_path):
    """One gene (GENEA) with a single annotated intron (100-900), plus rows
    exercising every classification branch against it, and one row with no
    gene annotation at all."""
    csv = tmp_path / "junction_type.csv"
    csv.write_text(
        "seqnames,start,end,strand,gencode_gene_id,gencode_gene_name,"
        "annotated,left_annotated,right_annotated\n"
        "chr1,100,900,+,ENSG_A,GENEA,1,,\n"                    # annotated
        "chr1,150,850,+,ENSG_A,GENEA,0,1,1\n"                  # isoform_switch
        "chr1,100,750,+,ENSG_A,GENEA,0,1,0\n"                  # alt_3p (donor known)
        "chr1,300,900,+,ENSG_A,GENEA,0,0,1\n"                  # alt_5p (acceptor known)
        "chr1,200,800,+,ENSG_A,GENEA,0,0,0\n"                  # novel_exon (nested in row 1)
        "chr1,2000,2100,+,ENSG_A,GENEA,0,0,0\n"                # novel
        "chr2,5000,6000,+,NA,NA,0,0,0\n"                       # no gene -> unlabelled
    )
    return csv


ROWS = [
    "chr1:100-900:+",
    "chr1:150-850:+",
    "chr1:100-750:+",
    "chr1:300-900:+",
    "chr1:200-800:+",
    "chr1:2000-2100:+",
    "chr2:5000-6000:+",
]


def test_classification_covers_every_category(tmp_path):
    idx = load_junction_type_index(_csv(tmp_path), tmp_path)
    assert idx.has_annotation_detail
    labels = idx.labels_for(ROWS)
    assert labels == [
        "GENEA:annotated",
        "GENEA:isoform_switch",
        "GENEA:alt_3p",
        "GENEA:alt_5p",
        "GENEA:novel_exon",
        "GENEA:novel",
        None,  # no gene annotation for this row
    ]


def test_labels_for_is_order_preserving_and_handles_unknown_rownames(tmp_path):
    idx = load_junction_type_index(_csv(tmp_path), tmp_path)
    labels = idx.labels_for(["chr1:2000-2100:+", "chr9:1-2:+", "chr1:100-900:+"])
    assert labels == ["GENEA:novel", None, "GENEA:annotated"]


def test_without_annotation_detail_columns_falls_back_to_unclassified(tmp_path):
    csv = tmp_path / "no_detail.csv"
    csv.write_text(
        "seqnames,start,end,strand,gencode_gene_id,gencode_gene_name\n"
        "chr1,100,900,+,ENSG_A,GENEA\n"
    )
    idx = load_junction_type_index(csv, tmp_path)
    assert not idx.has_annotation_detail
    assert idx.labels_for(["chr1:100-900:+"]) == ["GENEA:unclassified"]


def test_requires_coordinate_columns(tmp_path):
    csv = tmp_path / "bad.csv"
    csv.write_text("chrom,start,end\n1,2,3\n")
    with pytest.raises(JunctionMetadataError, match="missing column"):
        load_junction_type_index(csv, tmp_path)


def test_requires_gene_name_column(tmp_path):
    csv = tmp_path / "no_gene_name.csv"
    csv.write_text("seqnames,start,end,strand\nchr1,100,900,+\n")
    with pytest.raises(JunctionMetadataError, match="gene-name"):
        load_junction_type_index(csv, tmp_path)
