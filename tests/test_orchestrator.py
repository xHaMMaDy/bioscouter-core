import pytest

from bioscouter_core.models.unified import DataSource, OmicsType, UnifiedDataset
from bioscouter_core.orchestrator import BioScouterCoreSearch, SearchConfig


class DummyAdapter:
    def __init__(self, source, datasets):
        self.source = source
        self.datasets = datasets

    async def search(self, query, max_results=50, organism=None, **kwargs):
        return self.datasets[:max_results]

    async def close(self):
        return None


def dataset(source, accession, secondary=None, score=0.5):
    return UnifiedDataset(
        id=f"{source.value}:{accession}",
        accession=accession,
        source=source,
        source_url=f"https://example.org/{accession}",
        secondary_accession=secondary or [],
        omics_type=OmicsType.PROTEOMICS,
        title=f"Dataset {accession}",
        description="A normalized public dataset record.",
        organism=["Homo sapiens"],
        sample_count=20,
        relevance_score=score,
    )


@pytest.mark.asyncio
async def test_search_merges_cross_source_secondary_accessions():
    pride = dataset(DataSource.PRIDE, "PXD000001", secondary=["MSV000001"], score=0.8)
    massive = dataset(DataSource.MASSIVE, "MSV000001", secondary=["PXD000001"], score=0.9)
    search = BioScouterCoreSearch(
        config=SearchConfig(max_results=10),
        adapters={
            DataSource.PRIDE: DummyAdapter(DataSource.PRIDE, [pride]),
            DataSource.MASSIVE: DummyAdapter(DataSource.MASSIVE, [massive]),
        },
    )

    response = await search.search("TMT proteomics", sources=[DataSource.PRIDE, DataSource.MASSIVE])

    assert response.total_results == 1
    assert response.datasets[0].source == DataSource.PRIDE
    assert set(response.datasets[0].merged_sources) == {DataSource.PRIDE, DataSource.MASSIVE}
    assert response.datasets[0].relevance_score == 0.9


@pytest.mark.asyncio
async def test_search_keeps_adapter_failures_isolated():
    class FailingAdapter:
        async def search(self, query, max_results=50, organism=None, **kwargs):
            raise RuntimeError("source down")

    good = dataset(DataSource.PRIDE, "PXD000002", score=0.7)
    search = BioScouterCoreSearch(
        config=SearchConfig(max_results=10),
        adapters={
            DataSource.PRIDE: DummyAdapter(DataSource.PRIDE, [good]),
            DataSource.MASSIVE: FailingAdapter(),
        },
    )

    response = await search.search("proteomics", sources=[DataSource.PRIDE, DataSource.MASSIVE])

    assert response.total_results == 1
    assert response.results_by_source == {"pride": 1}

