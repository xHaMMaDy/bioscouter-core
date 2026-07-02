from bioscouter_core.models.unified import DataSource, OmicsType, detect_omics_types, get_sources_for_omics


def test_detects_proteomics_query():
    assert OmicsType.PROTEOMICS in detect_omics_types("TMT proteomics breast cancer")


def test_proteomics_routes_to_public_proteomics_sources():
    sources = get_sources_for_omics([OmicsType.PROTEOMICS])

    assert DataSource.PRIDE in sources
    assert DataSource.CPTAC in sources
    assert DataSource.MASSIVE in sources

