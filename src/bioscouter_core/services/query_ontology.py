"""Deterministic concept normalization for search/evaluation queries.

This module is intentionally conservative. It expands only curated biomedical
terms whose identifiers are stable enough for reviewer-facing evaluation and
keeps the original query text as the lead term.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Literal


LEXICON_VERSION = "bioscouter-query-ontology-2026-07-02"
ConceptCategory = Literal["disease", "tissue", "organism", "assay"]


@dataclass(frozen=True)
class Concept:
    category: ConceptCategory
    canonical: str
    identifiers: tuple[str, ...]
    synonyms: tuple[str, ...]


@dataclass(frozen=True)
class QueryExpansion:
    original_query: str
    expanded_query: str
    matched_concepts: tuple[Concept, ...]
    version: str = LEXICON_VERSION

    @property
    def changed(self) -> bool:
        return self.expanded_query != self.original_query


CONCEPTS: tuple[Concept, ...] = (
    Concept(
        category="disease",
        canonical="Alzheimer disease",
        identifiers=("MONDO:0004975",),
        synonyms=("alzheimer disease", "alzheimer's disease", "alzheimers disease", "ad"),
    ),
    Concept(
        category="disease",
        canonical="breast cancer",
        identifiers=("MONDO:0007254", "DOID:1612"),
        synonyms=("breast cancer", "breast carcinoma", "mammary cancer"),
    ),
    Concept(
        category="disease",
        canonical="lung cancer",
        identifiers=("MONDO:0008903", "DOID:1324"),
        synonyms=("lung cancer", "lung carcinoma", "pulmonary cancer"),
    ),
    Concept(
        category="tissue",
        canonical="brain",
        identifiers=("UBERON:0000955",),
        synonyms=("brain", "cerebral", "cortex", "hippocampus"),
    ),
    Concept(
        category="tissue",
        canonical="liver",
        identifiers=("UBERON:0002107",),
        synonyms=("liver", "hepatic"),
    ),
    Concept(
        category="tissue",
        canonical="blood",
        identifiers=("UBERON:0000178",),
        synonyms=("blood", "whole blood", "peripheral blood", "pbmc", "pbmcs"),
    ),
    Concept(
        category="organism",
        canonical="Homo sapiens",
        identifiers=("NCBITaxon:9606",),
        synonyms=("human", "humans", "homo sapiens", "h. sapiens"),
    ),
    Concept(
        category="organism",
        canonical="Mus musculus",
        identifiers=("NCBITaxon:10090",),
        synonyms=("mouse", "mice", "murine", "mus musculus", "m. musculus"),
    ),
    Concept(
        category="assay",
        canonical="RNA-seq",
        identifiers=("EFO:0008896",),
        synonyms=("rna-seq", "rnaseq", "rna seq", "transcriptome sequencing"),
    ),
    Concept(
        category="assay",
        canonical="single-cell RNA-seq",
        identifiers=("EFO:0008913",),
        synonyms=("scrna-seq", "single-cell rna-seq", "single cell rna seq", "single-cell transcriptomics"),
    ),
    Concept(
        category="assay",
        canonical="tandem mass tag proteomics",
        identifiers=("LOCAL_ASSAY:TMT",),
        synonyms=("tmt", "tmt proteomics", "tandem mass tag"),
    ),
    Concept(
        category="assay",
        canonical="liquid chromatography-mass spectrometry",
        identifiers=("LOCAL_ASSAY:LCMS",),
        synonyms=("lc-ms", "lc ms", "liquid chromatography mass spectrometry"),
    ),
)


def _tokens(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.casefold()).strip()


def _contains_phrase(query_tokens: str, phrase: str) -> bool:
    phrase_tokens = _tokens(phrase)
    if not phrase_tokens:
        return False
    return f" {phrase_tokens} " in f" {query_tokens} "


def match_concepts(query: str, concepts: Iterable[Concept] = CONCEPTS) -> tuple[Concept, ...]:
    """Return curated concepts matched in a query."""
    query_tokens = _tokens(query)
    matches: list[Concept] = []
    for concept in concepts:
        if any(_contains_phrase(query_tokens, synonym) for synonym in concept.synonyms):
            matches.append(concept)
    return tuple(matches)


def expand_query(query: str, *, max_terms_per_concept: int = 4) -> QueryExpansion:
    """Expand a search query with curated canonical terms and identifiers."""
    original = query.strip()
    matched = match_concepts(original)
    if not matched:
        return QueryExpansion(original_query=original, expanded_query=original, matched_concepts=())

    original_tokens = _tokens(original)
    seen = {original_tokens}
    additions: list[str] = []
    for concept in matched:
        terms = [concept.canonical, *concept.identifiers, *concept.synonyms[:max_terms_per_concept]]
        for term in terms:
            term_key = _tokens(term)
            if term_key and term_key not in seen and not _contains_phrase(original_tokens, term):
                seen.add(term_key)
                additions.append(term)

    expanded = " ".join([original, *additions]).strip()
    return QueryExpansion(
        original_query=original,
        expanded_query=expanded,
        matched_concepts=matched,
    )
