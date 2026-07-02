# Code Availability

Suggested manuscript text:

> The scientific core of BioScouter is available at `https://github.com/xHaMMaDy/bioscouter-core` under the MIT license. The repository includes the public data-source adapters, unified dataset schemas, deterministic query expansion, metadata-readiness scoring, faceting, search orchestration, benchmark queries, frozen benchmark outputs, relevance-labeling protocol, returned independent-label summaries, and scoring scripts. The deployed web platform at `https://bioscouter.com` is a managed service; production-only components such as authentication, user accounts, credits, payment handling, admin tools, deployment configuration, private observability, and production secrets are not included.

This URL should be mirrored in the manuscript Code Availability statement.

## Reviewer Access

Reviewers can evaluate the scientific code without using the production platform:

```bash
git clone https://github.com/xHaMMaDy/bioscouter-core.git
cd bioscouter-core
python -m pip install -e ".[dev]"
python -m pytest
bioscouter-core-search "TMT proteomics breast cancer" --source pride --max-results 10
```

If reviewers need to compare against the hosted system, provide the reviewer coupon/token separately in the manuscript submission portal or confidential reviewer notes, not in the public repository.
