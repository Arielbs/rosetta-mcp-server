# Rosetta MCP Server: Training Data Analysis & Optimization Plan

**Author:** Ariel J. Ben-Sasson (molCore)
**Date:** April 2026
**Purpose:** Document all Q&A data sources for Rosetta/PyRosetta/Biotite, quantify their relevance, and design an ML pipeline to optimize the MCP server's response quality.

---

## 1. Research Methodology

### 1.1 Search Strategy

We searched for Q&A content across 15+ platforms using three keyword filters:
- **Primary:** `rosetta` (protein modeling software, not E. coli strains or Rosetta Stone)
- **Secondary:** `pyrosetta` (Python interface)
- **Tertiary:** `biotite` (Python structural biology library)

For each platform, we measured:
- Total content volume (questions, answers, posts)
- Rosetta/PyRosetta/Biotite-specific content (filtered)
- Quality signals available (accepted answers, votes, resolution status)
- Scrapability (API, HTML, blocked)

### 1.2 Platforms Investigated

| # | Platform | Type | Method Used | Date Checked |
|---|----------|------|-------------|--------------|
| 1 | RosettaCommons Forum | Domain forum | WebFetch (HTML) | Apr 2026 |
| 2 | GitHub RosettaCommons/rosetta | Issues + Discussions | GitHub API | Apr 2026 |
| 3 | GitHub PyRosetta.notebooks | Issues | GitHub API | Apr 2026 |
| 4 | GitHub (all repos) | Issues search | GitHub Search API | Apr 2026 |
| 5 | PyRosetta.org FAQ | Official FAQ | WebFetch | Apr 2026 |
| 6 | RosettaCommons Docs | Documentation | WebFetch | Apr 2026 |
| 7 | Rosetta Workshops (Meiler Lab) | Tutorials | WebSearch | Apr 2026 |
| 8 | BioStars | Bioinformatics Q&A | site: search | Apr 2026 |
| 9 | ResearchGate | Academic Q&A | Topic pages | Apr 2026 |
| 10 | StackOverflow | Dev Q&A | API search | Apr 2026 |
| 11 | Bioinformatics StackExchange | Domain Q&A | API search | Apr 2026 |
| 12 | Reddit (r/bioinformatics, r/computationalbiology) | Social Q&A | API search | Apr 2026 |
| 13 | Quora | General Q&A | site: search | Apr 2026 |
| 14 | CCL (Computational Chemistry List) | Mailing list archive | WebSearch | Apr 2026 |
| 15 | Bioinformatics Discord/Slack | Chat | WebSearch | Apr 2026 |

---

## 2. Source-by-Source Analysis

### 2.1 RosettaCommons Forum (forum.rosettacommons.org)

**Status:** Read-only archive (closed mid-2024)
**URL:** https://forum.rosettacommons.org/forum/

| Category | Topics | Posts |
|----------|--------|-------|
| Rosetta 3 | 2,070 | 12,446 |
| PyRosetta | 628 | 2,473 |
| ROSIE | 207 | 523 |
| Rosetta++ | 206 | 585 |
| **Total** | **3,111** | **16,027** |

**Relevance to MCP:** 100% -- every question is about Rosetta/PyRosetta
**Quality signals:** Reply count, view count, thread resolution (last poster = original asker thanking)
**Scrapability:** HTML, read-only, stable URLs (node/{id} format)
**Key insight:** The PyRosetta category (628 topics) is the most directly relevant -- these are Python API questions that map directly to `pyrosetta_introspect` and `get_rosetta_help` tool gaps.

**Sample question types:**
- "How do I set up a MoveMap for FastRelax?" → tests `get_rosetta_help`
- "What's the difference between ref2015 and talaris2014?" → tests `get_rosetta_info`
- "How to convert XML to PyRosetta?" → tests `xml_to_pyrosetta`
- "InterfaceAnalyzer returning 0" → tests error guidance

### 2.2 GitHub RosettaCommons/rosetta

**URL:** https://github.com/RosettaCommons/rosetta

| Type | Count | Status |
|------|-------|--------|
| Discussions | 374 | Active (new Q&A go here since forum closed) |
| Issues (open) | 38 | Active |
| Issues (closed) | 34 | Resolved |
| **Total** | **446** | |

**Relevance to MCP:** 100%
**Quality signals:** Upvotes, marked answers (discussions), closed status + labels (issues), author (core dev vs user)
**Scrapability:** GitHub API (paginated, 100/page, 4+ pages of discussions)
**Key insight:** Discussions are the current replacement for the forum. Most recent and active source. Error messages in issues are valuable for improving `validate_xml` and error handling.

### 2.3 GitHub PyRosetta.notebooks Issues

**URL:** https://github.com/RosettaCommons/PyRosetta.notebooks/issues

| Count | 34 issues |
|-------|-----------|

**Relevance:** 100% -- PyRosetta-specific Python questions
**Quality signals:** Closed = resolved
**Scrapability:** GitHub API
**Key insight:** Small but high-quality -- these are notebook/tutorial-related issues, often about API changes or missing documentation.

### 2.4 GitHub Global Search (all repos)

Searched across ALL GitHub repos for issues mentioning these keywords:

| Keyword | Issues Found |
|---------|-------------|
| `rosetta protein` | 323 |
| `pyrosetta` | 448 |
| `biotite protein` | 131 |
| **Total (with overlap)** | **~700 unique** |

**Relevance:** ~60% -- many are from downstream projects using Rosetta, not Rosetta itself
**Quality signals:** Closed status, labels
**Scrapability:** GitHub Search API (rate limited)
**Key insight:** The 448 `pyrosetta` issues across all repos represent real-world usage patterns. Many are in research lab repos where users hit API issues.

### 2.5 PyRosetta FAQ

**URL:** https://www.pyrosetta.org/faq

| Count | ~42 entries |
|-------|-------------|

**Relevance:** 100% -- curated by PyRosetta maintainers
**Quality signals:** Highest quality -- official answers
**Scrapability:** HTML
**Key insight:** These should be directly embedded into `get_rosetta_help` extended responses. Gold standard answers.

### 2.6 RosettaCommons Documentation

**URL:** https://docs.rosettacommons.org/docs/latest/

| Type | Estimated Pages |
|------|----------------|
| Application docs | ~100 |
| Scripting docs | ~80 |
| Rosetta basics | ~60 |
| Developer docs | ~50 |
| Getting started/FAQ | ~20 |
| **Total** | **~300+ pages** |

**Relevance:** 100% -- authoritative reference
**Quality signals:** Official, versioned
**Scrapability:** HTML
**Key insight:** Already partially used by `get_rosetta_help` (fetches live docs). Could index for keyword mapping to improve alias resolution.

### 2.7 BioStars

**URL:** https://www.biostars.org

| Search | Estimated Matches |
|--------|-------------------|
| `rosetta` tag/keyword | ~50-100 questions |
| `pyrosetta` keyword | ~20-30 questions |
| `biotite` keyword | ~10-20 questions |
| **Total** | **~80-150** |

**Relevance:** ~70% -- some "rosetta" results are about E. coli Rosetta strains
**Quality signals:** Votes, accepted answers, view count (3.1K views on top question)
**Scrapability:** API blocked (403), HTML scrapable with care
**Key insight:** Cross-platform questions are valuable -- users asking "should I use Rosetta or AlphaFold?" reveal what concepts need comparison answers.

### 2.8 ResearchGate

**URL:** https://www.researchgate.net/topic/

| Topic | Total Q&As | Rosetta-specific |
|-------|-----------|-----------------|
| Protein Design | 35 | ~5-8 |
| Protein Engineering | 332 | ~10-15 |
| Protein Structure | 925 | ~15-20 |
| Protein Modeling | 280 | ~20-30 |
| Protein Folding | 285 | ~10-15 |
| Protein-Protein Interaction | 817 | ~5-10 |
| Proteins (general) | 1,708 | ~5-10 |
| **Total** | **4,382** | **~70-108** |

**Relevance:** ~2% of total (most are general biology, not Rosetta-specific)
**Quality signals:** Expert answers (with credentials), follower count
**Scrapability:** HTML (no API for Q&A)
**Key insight:** Low volume but high quality. Questions like "How to directly get the score of a protein structure using Rosetta energy function?" with expert answers are perfect training data. NOT worth bulk scraping -- cherry-pick the ~70-100 Rosetta-specific ones.

### 2.9 StackOverflow

| Search | Results |
|--------|---------|
| `rosetta protein` | 0 |
| `pyrosetta` | 0 |

**Relevance:** None
**Key insight:** Rosetta community doesn't use StackOverflow. All Q&A happens on domain-specific platforms.

### 2.10 Bioinformatics StackExchange

| Search | Results |
|--------|---------|
| `rosetta protein` | 0 |
| `pyrosetta` | 1 |

**Relevance:** Negligible

### 2.11 Reddit

| Subreddit | Rosetta-specific posts |
|-----------|----------------------|
| r/bioinformatics | ~10-20 |
| r/computationalbiology | ~5-10 |
| r/proteinengineering | ~5-10 |
| **Total** | **~20-40** |

**Relevance:** Low -- mostly "which tool should I use?" meta-questions
**Quality signals:** Upvotes, comment depth
**Scrapability:** Reddit API
**Key insight:** Useful for understanding how beginners phrase Rosetta questions (alias mining) but answers are often shallow.

### 2.12 Quora

| Rosetta-specific questions | ~3-5 |
|---------------------------|------|

**Relevance:** Negligible

### 2.13 CCL (Computational Chemistry List)

**URL:** https://ccl.net/
**Status:** Shutting down mid-2025, 30-year archive

| Total messages | ~100,000+ |
|---------------|-----------|
| Rosetta-specific | ~50-100 |

**Relevance:** ~0.1% of total
**Key insight:** Mostly GROMACS/AMBER/CHARMM discussions. Not worth mining for Rosetta specifically.

### 2.14 Discord/Slack

| Platform | Rosetta-specific messages |
|----------|-------------------------|
| Bioinformatics Discord | ~30-50 (estimated) |
| r/bioinformatics Slack | ~20-30 (estimated) |

**Relevance:** Low, ephemeral
**Key insight:** Not worth the access effort for the volume.

---

## 3. Consolidated Data Inventory

### 3.1 Final Source Priority

| Priority | Source | Rosetta-specific Q&As | Quality | Effort to Scrape |
|----------|--------|----------------------|---------|-----------------|
| **P0** | RosettaCommons Forum | 3,111 | High | Medium (HTML) |
| **P0** | GitHub Discussions | 374 | High | Low (API) |
| **P1** | GitHub Issues (RosettaCommons) | 72 | High | Low (API) |
| **P1** | PyRosetta FAQ | 42 | Highest | Low (HTML) |
| **P1** | GitHub Issues (global pyrosetta) | ~448 | Medium | Medium (API, rate limits) |
| **P2** | BioStars | ~100 | Medium | Medium (HTML scrape) |
| **P2** | ResearchGate (cherry-pick) | ~80 | High | Medium (HTML) |
| **P2** | RosettaCommons Docs (index) | ~300 pages | Authoritative | Low (HTML) |
| **P3** | Reddit | ~30 | Low | Low (API) |
| **P3** | GitHub Issues (global biotite) | ~131 | Medium | Medium (API) |
| **Skip** | StackOverflow, Quora, CCL, Discord/Slack | <10 each | Low | Not worth |

### 3.2 Total Actionable Data

| Category | Questions | Answers/Posts | Unique Q&A Pairs |
|----------|-----------|---------------|-----------------|
| P0 (must scrape) | 3,485 | ~17,000 | ~3,485 |
| P1 (should scrape) | ~562 | ~1,500 | ~562 |
| P2 (nice to have) | ~510 | ~2,000 | ~510 |
| P3 (low priority) | ~161 | ~400 | ~161 |
| **Total** | **~4,718** | **~20,900** | **~4,718** |

---

## 4. Data Schema

### 4.1 Structured Q&A Record

```json
{
  "id": "rosetta_forum_10538",
  "source": "rosettacommons_forum",
  "source_url": "https://forum.rosettacommons.org/node/10538",
  "category": "rosetta3",
  "title": "InterfaceAnalyzer total score = 0.000",
  "question_text": "I'm running InterfaceAnalyzerMover on my complex and all scores return 0...",
  "question_author": "user123",
  "question_date": "2020-03-15",
  "answers": [
    {
      "text": "You need to specify the interface string correctly. Use 'A_B' format...",
      "author": "jadolfbr",
      "author_is_core_dev": true,
      "date": "2020-03-16",
      "is_accepted": true,
      "votes": 3,
      "has_code_snippet": true,
      "code_snippets": ["iam = InterfaceAnalyzerMover('A_B')"]
    }
  ],
  "tags": ["InterfaceAnalyzerMover", "scoring", "interface", "mover"],
  "view_count": 2341,
  "reply_count": 4,
  "resolved": true,
  "resolution_signal": "accepted_answer",

  "mcp_relevance": {
    "relevant_tools": ["get_rosetta_help", "pyrosetta_introspect"],
    "current_tool_handles": false,
    "gap_type": "missing_alias",
    "suggested_alias": "InterfaceAnalyzer -> InterfaceAnalyzerMover",
    "suggested_help_text": "Ensure interface string uses chain separator format 'A_B'"
  },

  "quality_score": 0.92,
  "quality_signals": {
    "has_accepted_answer": true,
    "answer_by_core_dev": true,
    "has_working_code": true,
    "high_view_count": true,
    "resolved": true
  }
}
```

### 4.2 Quality Score Calculation

```
quality_score = weighted_sum(
  has_accepted_answer     * 0.25,   # Explicit resolution
  answer_by_expert        * 0.20,   # Author credibility
  has_working_code        * 0.20,   # Actionable answer
  resolved                * 0.15,   # Thread resolved
  high_view_count (>500)  * 0.10,   # Community interest
  multiple_confirming     * 0.10,   # Consensus
)
```

Expert detection: match author against known RosettaCommons core developers (from GitHub contributor list).

---

## 5. Train/Test/Validation Split Strategy

### 5.1 Split Ratios

| Set | Ratio | Purpose | Count (est.) |
|-----|-------|---------|-------------|
| **Training** | 70% | Mine aliases, expand help, find gaps | ~3,300 |
| **Validation** | 15% | Tune quality thresholds, test alias coverage | ~710 |
| **Test** | 15% | Final evaluation -- never seen during development | ~710 |
| **Total** | 100% | | ~4,720 |

### 5.2 Split Methodology

**NOT random.** Split by source and time to avoid data leakage:

| Set | Sources | Rationale |
|-----|---------|-----------|
| **Training** | Forum (Rosetta 3, ROSIE, Rosetta++), GitHub Discussions (pre-2025), BioStars, ResearchGate | Largest, oldest, most comprehensive |
| **Validation** | Forum (PyRosetta category), PyRosetta FAQ, GitHub Issues | PyRosetta-specific -- tests tool coverage for Python API |
| **Test** | GitHub Discussions (2025-2026), GitHub global issues (pyrosetta keyword) | Most recent, never-seen questions from current users |

This ensures:
- Training on historical data
- Validating on curated/PyRosetta-specific data
- Testing on recent, real-world questions the MCP will actually face

### 5.3 Cross-validation for Alias Mining

For alias mining specifically, use 5-fold cross-validation within the training set:
- Fold 1-4: Extract candidate aliases from question titles
- Fold 5: Test if the alias resolves correctly
- Rotate folds

---

## 6. Training Pipeline

### 6.1 Overview

```
Scrape → Structure → Analyze → Extract → Optimize → Evaluate
  P0-P2     Schema     Gap       Aliases    Code       Test
 sources    4.2        detect    + help     changes    suite
```

### 6.2 Phase 1: Scrape & Structure

**Input:** Raw HTML/JSON from sources
**Output:** Structured Q&A records (schema 4.1)

Scripts to build:
```
scripts/scrape/
  forum_scraper.py        # RosettaCommons forum (3,111 topics)
  github_discussions.py   # GitHub Discussions API
  github_issues.py        # GitHub Issues API (RosettaCommons + global)
  pyrosetta_faq.py        # PyRosetta FAQ page
  biostars_scraper.py     # BioStars HTML scrape
  researchgate_picker.py  # Cherry-pick Rosetta-specific Q&As
  docs_indexer.py         # Index RosettaCommons docs pages for keywords
```

**Output format:** `data/raw/{source}/{id}.json`
**Deduplicated:** `data/structured/all_qa.jsonl` (one record per line)

### 6.3 Phase 2: Gap Analysis

**Input:** Structured Q&A records
**Output:** Gap report showing where the MCP fails

For each question:
1. Extract the core topic/concept from the title
2. Simulate MCP tool calls:
   - `get_rosetta_help(topic)` -- does it return useful content?
   - `rosetta_to_biotite(topic)` -- does it find a mapping?
   - `xml_to_pyrosetta` -- if XML is present, can it translate?
3. Compare MCP response to the known good answer
4. Classify gap type:

| Gap Type | Example | Fix |
|----------|---------|-----|
| `missing_alias` | "relax" doesn't resolve to FastRelax | Add to TOPIC_ALIASES |
| `missing_biotite_mapping` | "contact order" has no Biotite entry | Add to BIOTITE_ROSETTA_MAPPINGS |
| `missing_xml_element` | `<LoopModeler>` not in translator | Add to XML_TO_PYROSETTA_MAPPINGS |
| `shallow_help` | Static help is 1 sentence | Expand EXTENDED_HELP |
| `missing_url_pattern` | Docs URL probe fails | Add URL pattern to searchRosettaWebDocs |
| `wrong_answer` | MCP gives incorrect API signature | Fix mapping data |
| `no_gap` | MCP already handles this correctly | No action |

**Output:** `data/analysis/gap_report.json`

```json
{
  "total_questions": 4718,
  "no_gap": 4200,
  "gaps": {
    "missing_alias": 180,
    "missing_biotite_mapping": 45,
    "missing_xml_element": 30,
    "shallow_help": 120,
    "missing_url_pattern": 60,
    "wrong_answer": 5
  }
}
```

### 6.4 Phase 3: Extract Improvements

**Input:** Gap report + Q&A records
**Output:** Candidate improvements to MCP code

#### 3a. Alias Mining

From question titles, extract natural language → Rosetta concept mappings:

```python
# Input: "How do I relax my protein structure?"
# Output: {"relax": "FastRelax", "relax protein": "FastRelax"}

# Input: "Calculating binding energy with Rosetta"
# Output: {"binding energy": "Ddg", "calculate binding": "InterfaceAnalyzerMover"}
```

Method:
1. Tokenize question titles
2. Match tokens against known Rosetta class names
3. Extract the surrounding natural language as alias candidates
4. Rank by frequency (more questions using the same phrasing = better alias)
5. Manually review top candidates

#### 3b. Help Content Expansion

For questions where `get_rosetta_help` fails but a high-quality forum answer exists:
1. Extract the best answer (highest quality_score)
2. Summarize into a help snippet (200-500 chars)
3. Add to EXTENDED_HELP dictionary

#### 3c. Biotite Mapping Expansion

For questions comparing Rosetta to other tools:
1. Identify the Rosetta function and the alternative
2. If the alternative is Biotite or has a Biotite equivalent, create a new mapping entry
3. Include code examples from the answers

#### 3d. XML Element Expansion

For XML-related questions where `xml_to_pyrosetta` can't translate:
1. Identify the XML element name
2. Look up the PyRosetta class name from the answer
3. Add to XML_TO_PYROSETTA_MAPPINGS with attribute setters

### 6.5 Phase 4: Optimize Server Code

Apply extracted improvements to `rosetta_mcp_wrapper.js`:

| Target | What changes | File location |
|--------|-------------|---------------|
| TOPIC_ALIASES | Add ~100-200 new aliases | `getRosettaHelp()` |
| SEARCH_ALIASES | Add ~50-100 new Biotite search terms | `rosettaToBiotite()` |
| EXTENDED_HELP | Expand 5 topics to 20+ | `getRosettaHelp()` |
| BIOTITE_ROSETTA_MAPPINGS | Add ~5-10 new entries | Top-level const |
| XML_TO_PYROSETTA_MAPPINGS | Add ~10-15 new elements | Top-level const |
| URL patterns | Add ~5-10 new doc URL patterns | `searchRosettaWebDocs()` |

### 6.6 Phase 5: Evaluate

Run the test set (710 never-seen questions) through the optimized MCP:

**Metrics:**

| Metric | Definition | Current | Target |
|--------|-----------|---------|--------|
| **Hit rate** | % of questions where MCP returns useful content | ~97% | 99%+ |
| **Alias coverage** | % of natural language phrases that resolve | ~60% | 90%+ |
| **Help quality** | Average response length for valid queries | ~2000 chars | ~3000 chars |
| **Biotite coverage** | % of Biotite-related queries that find mappings | ~85% | 95%+ |
| **XML coverage** | % of XML elements the translator handles | ~37/100+ | 60+ |
| **False positive rate** | % of responses that give wrong information | <1% | <0.5% |

**Evaluation method:**
1. For each test question, call the relevant MCP tool
2. Compare response to the known-good forum answer
3. Score: 1 (correct and helpful), 0.5 (partially correct), 0 (wrong or missing)
4. Compute weighted average by quality_score of the test question

---

## 7. Timeline

| Week | Phase | Deliverable |
|------|-------|-------------|
| 1 | Scrape P0 sources | 3,485 structured Q&A records |
| 1 | Scrape P1 sources | 562 additional records |
| 2 | Gap analysis | Gap report with counts per category |
| 2 | Alias mining | 200+ candidate aliases, manually reviewed |
| 3 | Code optimization | Updated aliases, help, mappings |
| 3 | Evaluation on test set | Hit rate and coverage metrics |
| 4 | Scrape P2 sources | Additional 510 records |
| 4 | Second optimization round | Refinements based on P2 data |
| Ongoing | Monitor GitHub Discussions | Monthly alias/mapping updates |

---

## 8. Infrastructure

### 8.1 Storage

```
data/
  raw/                          # Raw scraped HTML/JSON per source
    rosettacommons_forum/       # 3,111 files
    github_discussions/         # 374 files
    github_issues/              # 72 + 448 files
    pyrosetta_faq/              # 42 files
    biostars/                   # ~100 files
    researchgate/               # ~80 files
  structured/
    all_qa.jsonl                # Deduplicated, structured records
    train.jsonl                 # 70% split
    validation.jsonl            # 15% split
    test.jsonl                  # 15% split
  analysis/
    gap_report.json             # Gap analysis results
    alias_candidates.json       # Mined aliases
    coverage_report.json        # Evaluation metrics
```

### 8.2 Estimated Size

| Dataset | Records | Est. Size |
|---------|---------|-----------|
| Raw data | ~4,718 Q&As + 20,900 answers | ~50-100 MB |
| Structured | ~4,718 records | ~20 MB |
| Analysis outputs | Reports + candidates | ~5 MB |

### 8.3 Tools Required

- Python 3.10+ with `requests`, `beautifulsoup4`, `lxml`
- GitHub API token (for rate limits)
- No GPU needed -- this is data processing, not neural network training
- The "training" is data-driven expansion of lookup tables, not model fine-tuning

---

## 9. Key Insight: This is NOT Neural Network Training

The optimization approach is **lookup table expansion**, not model fine-tuning:

| Approach | What we do | What we don't do |
|----------|-----------|-----------------|
| Mine aliases from question titles | YES | Train an NLP model |
| Expand help text with forum answers | YES | Fine-tune an LLM |
| Add Biotite mappings from comparisons | YES | Train a seq2seq translator |
| Add XML elements from user questions | YES | Train a code generator |

The MCP server is a **deterministic tool** -- it returns data from lookup tables. Making it better means making the tables more complete, not training a neural network. The "ML" is in the data mining (clustering, pattern extraction), not in the server itself.

This is actually a strength: the responses are **reproducible, verifiable, and debuggable** -- unlike an LLM's probabilistic output.
