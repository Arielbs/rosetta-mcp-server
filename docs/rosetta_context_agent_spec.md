# Rosetta Context Agent: Project Specification

**Codename:** `rosetta-context-agent`
**Author:** Ariel J. Ben-Sasson (molCore)
**Date:** April 2026
**Status:** Specification (pre-development)
**Relationship:** Separate product from `rosetta-mcp-server` (which becomes one tool this agent can call)

---

## 1. Vision

A RAG-powered expert agent that answers Rosetta/PyRosetta/Biotite questions with the depth of a senior computational biophysicist -- grounded in 4,718 real community Q&As, 20,900 answers, and continuously improving from user feedback.

**Current MCP** = lookup tables → deterministic answers → fast, simple, limited
**This agent** = knowledge retrieval + LLM reasoning → expert-quality answers → deeper, adaptive, growing

---

## 2. Architecture

```
User Question
    │
    ▼
┌─────────────────────────────┐
│  Rosetta Context Agent      │
│  (MCP Server / HTTP)        │
│                             │
│  ┌───────────────────────┐  │
│  │ Tool Router            │  │
│  │ Decides: RAG? Lookup?  │  │
│  │ Both? Direct answer?   │  │
│  └───────┬───────────────┘  │
│          │                  │
│  ┌───────┼───────────────┐  │
│  │       │               │  │
│  ▼       ▼               ▼  │
│ ┌─────┐ ┌──────────┐ ┌────┐│
│ │RAG  │ │rosetta-  │ │Web ││
│ │Search│ │mcp-server│ │Docs││
│ │(Vec.)│ │(existing)│ │Fetch│
│ └──┬──┘ └────┬─────┘ └─┬──┘│
│    │         │          │   │
│    ▼         ▼          ▼   │
│  ┌───────────────────────┐  │
│  │ Context Assembler     │  │
│  │ Combines: RAG results │  │
│  │ + tool outputs + docs │  │
│  └───────────┬───────────┘  │
│              │              │
│              ▼              │
│  ┌───────────────────────┐  │
│  │ LLM (Claude API)      │  │
│  │ Synthesizes expert     │  │
│  │ answer from context    │  │
│  └───────────┬───────────┘  │
│              │              │
│              ▼              │
│  ┌───────────────────────┐  │
│  │ Feedback Collector    │  │
│  │ Logs Q, A, user       │  │
│  │ rating (if provided)  │  │
│  └───────────────────────┘  │
└─────────────────────────────┘
    │
    ▼
Answer + Sources
```

---

## 3. Separation from rosetta-mcp-server

| Aspect | rosetta-mcp-server | rosetta-context-agent |
|--------|-------------------|----------------------|
| Repo | `Arielbs/rosetta-mcp-server` | `Arielbs/rosetta-context-agent` (new) |
| Purpose | Deterministic tool lookups | Expert reasoning + RAG |
| Dependency | None | Calls rosetta-mcp-server tools + Claude API |
| Compute | Lightweight (Node.js) | Heavier (vector DB + LLM API calls) |
| Storage | Stateless | Vector DB + feedback DB |
| Cost | $5/mo | $30-80/mo (see Section 7) |
| Users | Anonymous, no account | Optional account for feedback |
| Data | No user data stored | Stores questions, answers, feedback |
| Privacy | Minimal | Needs proper data policy |

**The current MCP becomes a dependency**, not a competitor. The agent calls it for fast lookups (aliases, Biotite mappings, XML translation) and uses RAG for deeper questions.

---

## 4. Components

### 4.1 Knowledge Base (Vector DB)

**Contents:**
| Source | Records | Est. Tokens | Embedding Cost (one-time) |
|--------|---------|-------------|--------------------------|
| RosettaCommons Forum | 3,111 Q&As | ~3M tokens | $0.06 |
| GitHub Discussions | 374 | ~400K | $0.01 |
| GitHub Issues | ~520 | ~500K | $0.01 |
| PyRosetta FAQ | 42 | ~50K | <$0.01 |
| BioStars | ~100 | ~100K | <$0.01 |
| ResearchGate (cherry-picked) | ~80 | ~100K | <$0.01 |
| RosettaCommons Docs (indexed) | ~300 pages | ~1.5M | $0.03 |
| **Total** | **~4,527** | **~5.6M tokens** | **~$0.12** |

**Embedding model:** OpenAI `text-embedding-3-small` ($0.02/1M tokens)
**Vector dimensions:** 1536
**Total vectors:** ~4,527 (tiny -- fits in memory)

**Vector DB options:**
| Option | Cost | Pros | Cons |
|--------|------|------|------|
| **SQLite + sqlite-vss** | $0 (embedded) | Zero infra, ships with app | Limited features |
| **Qdrant (self-hosted on Railway)** | $5-10/mo | Full-featured, fast | Extra service |
| **Pinecone Serverless** | $0 free tier (up to 100K vectors) | Zero ops, scale to zero | Vendor lock-in |
| **ChromaDB (embedded)** | $0 (embedded) | Python-native, simple | Less mature |

**Recommendation:** Start with **ChromaDB embedded** or **sqlite-vss** (zero cost, ~5K vectors easily fits in memory). Move to Qdrant/Pinecone only if scale demands it.

### 4.2 LLM Integration (Claude API)

**When the agent needs Claude API:**
- Synthesizing answers from multiple RAG results
- Reasoning about biophysics problems
- Explaining error messages in context

**When it does NOT need Claude API:**
- Direct lookup answers (use rosetta-mcp-server tools)
- FAQ matches (return cached answer directly)
- Exact match in knowledge base (return forum answer)

**Estimated Claude API usage:**
| Scenario | Input tokens | Output tokens | Cost per call |
|----------|-------------|---------------|--------------|
| RAG synthesis (3 context docs + question) | ~3,000 | ~500 | ~$0.012 (Haiku) |
| Complex reasoning | ~5,000 | ~1,000 | ~$0.030 (Sonnet) |
| Simple lookup (no API needed) | 0 | 0 | $0.00 |

**Model selection strategy:**
- **Haiku 4.5** ($1/$5 per 1M tokens): For straightforward RAG synthesis
- **Sonnet 4.6** ($3/$15 per 1M tokens): For complex biophysics reasoning
- **No LLM**: For direct lookup matches (>50% of queries)

### 4.3 Feedback Collection

**What we collect per interaction:**
```json
{
  "id": "uuid",
  "timestamp": "2026-04-04T10:30:00Z",
  "question": "How do I relax my protein?",
  "answer": "FastRelax performs all-atom relaxation...",
  "sources_used": ["forum_10538", "pyrosetta_faq_12"],
  "tools_called": ["search_rosetta_knowledge", "get_rosetta_help"],
  "response_time_ms": 1200,
  "user_feedback": null,  // null = no feedback, 1-5 = rating, "text" = comment
  "user_id": "anonymous",  // or hashed session ID
  "model_used": "haiku-4.5"
}
```

**Feedback types:**
| Type | How collected | Used for |
|------|--------------|----------|
| Implicit | Response time, tool calls | Monitor system health |
| Explicit (rating) | Thumbs up/down or 1-5 stars | Track answer quality |
| Explicit (correction) | "Actually, the correct API is..." | Update knowledge base |
| Implicit (follow-up) | User asks a follow-up = first answer was incomplete | Identify gaps |

**Storage:** SQLite database alongside the vector DB. Lightweight, no separate service.

### 4.4 Self-Improvement Loop

```
Weekly batch job:
  1. Collect all feedback from past week
  2. Identify low-rated answers (rating < 3)
  3. Identify unanswered questions (no good RAG match)
  4. For corrections: update knowledge base entry
  5. For gaps: flag for manual review
  6. Re-embed updated entries
  7. Generate quality report
```

**NOT automatic retraining.** Human reviews flagged items before knowledge base changes. The system suggests improvements; a human approves them.

---

## 5. MCP Tools (the new agent exposes)

| Tool | Description | Uses LLM? |
|------|-------------|-----------|
| `ask_rosetta_expert` | Main tool: ask any Rosetta/PyRosetta/Biotite question, get expert answer with sources | Yes (RAG + synthesis) |
| `search_rosetta_knowledge` | Semantic search over 4,500+ Q&As, return top matches | No (vector search only) |
| `submit_feedback` | Rate an answer or provide correction | No |
| `get_similar_questions` | Find related questions from the knowledge base | No (vector search) |
| All existing `rosetta-mcp-server` tools | Proxied through for deterministic lookups | No |

---

## 6. Tech Stack

| Component | Technology | Why |
|-----------|-----------|-----|
| Runtime | Node.js or Python | Python preferred (ChromaDB, embeddings are Python-native) |
| Vector DB | ChromaDB (embedded) | Zero cost, ~5K vectors, Python-native |
| Embeddings | OpenAI `text-embedding-3-small` | Cheapest, good quality, 1536 dims |
| LLM | Claude API (Haiku for simple, Sonnet for complex) | Best reasoning, your existing relationship |
| Feedback DB | SQLite | Embedded, zero ops |
| Hosting | Railway ($5/mo hobby) | Same platform as current MCP |
| Existing MCP | Call via HTTP to `mcp.molcore.bio/mcp` | Reuse existing tools |

---

## 7. Cost Estimation

### 7.1 One-Time Setup Costs

| Item | Cost |
|------|------|
| Scraping + structuring 4,718 Q&As | Developer time (your time or agent) |
| Embedding 5.6M tokens | $0.12 |
| ChromaDB setup | $0 (embedded) |
| **Total one-time** | **~$0.12 + dev time** |

### 7.2 Monthly Operating Costs

| Item | Low usage (100 queries/mo) | Medium (1,000 queries/mo) | High (10,000 queries/mo) |
|------|--------------------------|--------------------------|--------------------------|
| Railway hosting | $5 | $5 | $10 |
| Claude API (Haiku, 50% of queries need LLM) | $0.60 | $6 | $60 |
| Claude API (Sonnet, 10% complex queries) | $0.30 | $3 | $30 |
| OpenAI Embeddings (new questions) | $0.01 | $0.02 | $0.20 |
| Storage (SQLite + ChromaDB) | $0 (included in Railway) | $0 | $0 |
| **Monthly total** | **~$6** | **~$14** | **~$100** |

### 7.3 Scaling Cost Drivers

| Scale factor | Impact |
|-------------|--------|
| More queries | Linear LLM API cost increase |
| Larger knowledge base (10K+ docs) | Negligible (embedding is cheap, ChromaDB handles it) |
| More feedback data | Negligible (SQLite handles millions of rows) |
| Upgrading to Sonnet for all queries | 5x LLM cost increase |

### 7.4 Comparison to Current MCP

| | Current MCP | New Agent (medium usage) |
|---|---|---|
| Monthly cost | $5 | ~$14 |
| LLM API calls | 0 | ~500/mo |
| Storage | 0 | ~100MB |
| Maintenance | None | Weekly feedback review (~1 hr) |

---

## 8. Self-Improvement Cost

### 8.1 Weekly Improvement Cycle

| Task | Automated? | Cost |
|------|-----------|------|
| Collect feedback | Yes | $0 |
| Re-embed updated entries | Yes | <$0.01/week |
| Generate quality report | Yes | <$0.01/week |
| Manual review of flagged items | No (human) | ~1 hour/week |
| **Weekly total** | | **<$0.05 + 1 hr human time** |

### 8.2 Quarterly Knowledge Base Refresh

| Task | Cost |
|------|------|
| Scrape new GitHub Discussions (since last scrape) | Dev time |
| Embed new entries (~50-100 per quarter) | <$0.01 |
| Review quality metrics, tune retrieval | 2-3 hours |
| **Quarterly total** | **~$0.01 + 3 hrs** |

---

## 9. Railway Deployment Estimate

### 9.1 Resource Requirements

| Component | CPU | Memory | Disk |
|-----------|-----|--------|------|
| Python app (FastAPI) | 0.5 vCPU | 512MB | - |
| ChromaDB (in-process) | 0.5 vCPU | 256MB (5K vectors) | 50MB |
| SQLite (feedback) | negligible | negligible | 10MB |
| **Total** | **1 vCPU** | **~768MB** | **~60MB** |

### 9.2 Railway Plan Fit

| Railway Plan | Cost | CPU | Memory | Fits? |
|-------------|------|-----|--------|-------|
| Hobby ($5/mo, $5 credit) | $5 | 8 shared vCPU | 8GB | Yes, easily |
| Pro ($20/mo) | $20 | Priority compute | 32GB | Overkill |

**The $5/mo Hobby plan handles this comfortably.** The app uses <1GB memory and <1 vCPU. The main cost is the Claude API, not the hosting.

---

## 10. Privacy & Data Considerations

| Data Type | Collected | Stored | Shared | Retention |
|-----------|-----------|--------|--------|-----------|
| User questions | Yes | Yes (SQLite) | No | 1 year |
| Agent answers | Yes | Yes (SQLite) | No | 1 year |
| User feedback (ratings) | Optional | Yes | Aggregated only | Indefinite |
| User corrections | Optional | Yes | May update KB | Indefinite |
| User identity | No (anonymous by default) | Session hash only | No | Session |
| PDB files / sequences | No (not uploaded) | No | No | N/A |

**Privacy policy needed:** Separate from current MCP. Must explain data collection for improvement purposes.

---

## 11. Open Source Considerations

| Aspect | Approach |
|--------|----------|
| Code | Open source (MIT) -- same as current MCP |
| Knowledge base | Open (scraped from public sources) |
| Feedback data | Private (user interactions, not published) |
| Embeddings | Reproducible (script to regenerate from public data) |
| Claude API key | User brings their own, or use molCore's for hosted version |

**Self-hosting option:** Anyone can clone, scrape, embed, and run with their own Claude API key. The hosted version at molCore adds convenience + shared feedback loop.

---

## 12. Implementation Phases

| Phase | Deliverable | Timeline | Cost |
|-------|------------|----------|------|
| 0. Scrape & embed knowledge base | 4,527 embedded Q&A vectors | 1 week | $0.12 |
| 1. Basic RAG tool (`search_rosetta_knowledge`) | Semantic search, no LLM | 1 week | $0 |
| 2. Expert answer tool (`ask_rosetta_expert`) | RAG + Claude synthesis | 1 week | ~$5/mo API |
| 3. Feedback collection | Rating + corrections | 3 days | $0 |
| 4. Self-improvement loop | Weekly batch, quality reports | 3 days | ~$0.05/week |
| 5. Deploy to Railway | `agent.molcore.bio` | 1 day | $5/mo |
| 6. Claude Skill companion | Packaged knowledge for claude.ai | 1 week | $0 |
| **Total** | | **~4-5 weeks** | **~$10-15/mo ongoing** |

---

## 13. Success Metrics

| Metric | Baseline (current MCP) | Target (agent) |
|--------|----------------------|----------------|
| Question hit rate | 97% | 99.5%+ |
| Answer quality (user rating) | N/A (no feedback) | 4.0+/5.0 avg |
| Average response depth | ~200 chars (lookup) | ~1000 chars (synthesized) |
| Sources cited per answer | 0-1 | 2-3 |
| Self-improvement rate | 0 (static) | ~5 KB entries updated/month |
| Time to answer new topic | Never (until manual update) | <1 week (auto-flagged) |
