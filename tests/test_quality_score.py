from bioscouter_core.models.unified import DataSource, OmicsType, UnifiedDataset
from bioscouter_core.services.quality_score import calculate_and_set_quality


def test_quality_score_is_normalized_fraction():
    dataset = UnifiedDataset(
        id="pride:PXD000001",
        accession="PXD000001",
        source=DataSource.PRIDE,
        source_url="https://www.ebi.ac.uk/pride/archive/projects/PXD000001",
        omics_type=OmicsType.PROTEOMICS,
        title="Example proteomics dataset",
        description="A sufficiently detailed proteomics dataset description for scoring.",
        organism=["Homo sapiens"],
        sample_count=12,
    )

    scored = calculate_and_set_quality(dataset)

    assert scored.quality_score is not None
    assert 0 <= scored.quality_score <= 1
    assert "metadata_readiness_breakdown" in scored.extensions

