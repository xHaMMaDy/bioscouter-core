"""Metadata Readiness Calculator for unified datasets.

The score summarizes discoverability and metadata availability. It is not a
measurement of experimental validity, biological quality, or reproducibility.
The serialized ``quality_score`` and ``quality_breakdown`` names are retained
for API compatibility with existing clients.
"""

import structlog
from typing import Optional, Dict, Any
from bioscouter_core.models.unified import UnifiedDataset, CurationLevel

logger = structlog.get_logger(__name__)


class QualityScoreCalculator:
    """Calculate a metadata-readiness heuristic from repository metadata."""
    
    # Weight distribution for metadata-readiness components
    WEIGHTS = {
        "sample_count": 0.20,       # Availability and scale of sample metadata
        "metadata_completeness": 0.25,  # Complete metadata = better documented
        "publications": 0.20,       # Linked publications improve traceability
        "curation_level": 0.15,     # Source curation status
        "data_availability": 0.10,  # Download links available
        "description_quality": 0.10,  # Description detail
    }
    
    # Sample count thresholds for scoring
    SAMPLE_THRESHOLDS = {
        "excellent": 100,  # >= 100 samples = full score
        "good": 50,        # >= 50 samples = 0.8
        "moderate": 20,    # >= 20 samples = 0.6
        "low": 5,          # >= 5 samples = 0.4
    }
    
    # Curation level scores
    CURATION_SCORES = {
        CurationLevel.CURATED: 1.0,      # Expert-curated
        CurationLevel.COMMUNITY: 0.65,   # Community-submitted with validation
        CurationLevel.AUTO: 0.35,        # Auto-indexed
    }
    
    def calculate_quality_score(self, dataset: UnifiedDataset) -> Dict[str, Any]:
        """
        Calculate a metadata-readiness score for a dataset.
        
        Returns:
            Dict with:
                - score: Overall metadata-readiness score (0.0 - 1.0)
                - breakdown: Individual component scores
                - is_estimated: True (always, since calculated from metadata)
        """
        breakdown = {}
        
        # 1. Sample count score
        breakdown["sample_count"] = self._score_sample_count(dataset.sample_count)
        
        # 2. Metadata completeness
        breakdown["metadata_completeness"] = self._score_metadata_completeness(dataset)
        
        # 3. Publications
        breakdown["publications"] = self._score_publications(dataset)
        
        # 4. Curation level
        breakdown["curation_level"] = self._score_curation_level(dataset.curation_level)
        
        # 5. Data availability
        breakdown["data_availability"] = self._score_data_availability(dataset)
        
        # 6. Description quality
        breakdown["description_quality"] = self._score_description_quality(dataset.description)
        
        # Calculate weighted average
        total_score = sum(
            breakdown[key] * self.WEIGHTS[key]
            for key in self.WEIGHTS
        )
        
        return {
            "score": round(total_score, 3),
            "breakdown": breakdown,
            "is_estimated": True,
        }
    
    def _score_sample_count(self, sample_count: int) -> float:
        """Score based on number of samples."""
        if sample_count >= self.SAMPLE_THRESHOLDS["excellent"]:
            return 1.0
        elif sample_count >= self.SAMPLE_THRESHOLDS["good"]:
            return 0.8
        elif sample_count >= self.SAMPLE_THRESHOLDS["moderate"]:
            return 0.6
        elif sample_count >= self.SAMPLE_THRESHOLDS["low"]:
            return 0.4
        elif sample_count > 0:
            return 0.2
        return 0.0
    
    def _score_metadata_completeness(self, dataset: UnifiedDataset) -> float:
        """Score based on how complete the metadata is."""
        fields_to_check = [
            bool(dataset.title),
            bool(dataset.description),
            bool(dataset.organism),
            bool(dataset.omics_type),
            bool(dataset.submission_date),
            bool(dataset.contributors),
            bool(dataset.tissue),
            bool(dataset.disease),
        ]
        
        filled_count = sum(fields_to_check)
        return filled_count / len(fields_to_check)
    
    def _score_publications(self, dataset: UnifiedDataset) -> float:
        """Score based on publication links."""
        has_pubmed = len(dataset.pubmed_ids) > 0
        has_doi = bool(dataset.doi)
        has_citation = bool(dataset.citation)
        
        if has_pubmed and has_doi:
            return 1.0
        elif has_pubmed or has_doi:
            return 0.7
        elif has_citation:
            return 0.4
        return 0.0
    
    def _score_curation_level(self, curation_level: Optional[CurationLevel]) -> float:
        """Score based on curation status from source."""
        if curation_level is None:
            return 0.5  # Unknown = neutral score
        return self.CURATION_SCORES.get(curation_level, 0.5)
    
    def _score_data_availability(self, dataset: UnifiedDataset) -> float:
        """Score based on download link availability."""
        if not dataset.download_links:
            return 0.3  # Some data is usually available even without explicit links
        
        link_count = len(dataset.download_links)
        if link_count >= 5:
            return 1.0
        elif link_count >= 2:
            return 0.8
        return 0.6
    
    def _score_description_quality(self, description: Optional[str]) -> float:
        """Score based on description length/detail."""
        if not description:
            return 0.0
        
        length = len(description)
        if length >= 500:
            return 1.0
        elif length >= 200:
            return 0.8
        elif length >= 100:
            return 0.6
        elif length >= 50:
            return 0.4
        return 0.2


# Singleton instance
_calculator_instance: Optional[QualityScoreCalculator] = None


def get_quality_calculator() -> QualityScoreCalculator:
    """Get the singleton metadata-readiness calculator instance."""
    global _calculator_instance
    if _calculator_instance is None:
        _calculator_instance = QualityScoreCalculator()
    return _calculator_instance


def calculate_and_set_quality(dataset: UnifiedDataset) -> UnifiedDataset:
    """
    Calculate metadata readiness and store it in compatibility fields.
    
    This is a convenience function to be called after creating a UnifiedDataset.
    Only sets quality_score if it's not already set.
    
    Returns:
        The dataset with metadata-readiness fields populated
    """
    if dataset.quality_score is not None:
        return dataset
    
    calculator = get_quality_calculator()
    result = calculator.calculate_quality_score(dataset)
    dataset.quality_score = result["score"]
    
    # Store breakdown in extensions for potential UI use
    if "quality_breakdown" not in dataset.extensions:
        dataset.extensions["quality_breakdown"] = result["breakdown"]
        dataset.extensions["quality_is_estimated"] = True
    dataset.extensions.setdefault(
        "metadata_readiness_breakdown",
        result["breakdown"],
    )
    dataset.extensions.setdefault("metadata_readiness_is_heuristic", True)
    
    return dataset
