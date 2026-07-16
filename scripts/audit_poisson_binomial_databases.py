#!/usr/bin/env python3
"""Reproducible bibliographic search for the finite first-descent theorem.

This script records discovery metadata and a compact primary-text screening
ledger.  Detailed mathematical exclusions remain in the companion audit note;
neither layer may infer exclusions from database snippets alone.
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any


USER_AGENT = "erdos-993-literature-audit/1.0 (reproducible scholarly search)"
CUTOFF = "2026-07-16"
CUTOFF_DATE = date.fromisoformat(CUTOFF)
CUTOFF_ARXIV = CUTOFF.replace("-", "") + "2359"


ZBMATH_QUERIES = [
    "Poisson binomial mode variance",
    "Poisson binomial adjacent ratios",
    "Poisson binomial Turan inequality",
    "Poisson binomial discrete curvature",
    "Poisson binomial strong log concavity",
    "Poisson binomial mode",
    "Poisson binomial log concavity",
    "Bernoulli sum mode variance",
    "three consecutive probabilities variance",
    "adjacent probabilities mode variance",
    "Turan determinant Bernoulli",
    "real-rooted coefficient ratio variance",
]

ZBMATH_FREE_QUERIES = [
    '"Poisson binomial" & mode & variance',
    'ab:"Poisson binomial" & (mode | modal)',
    'ab:"Poisson binomial" & (ratio | curvature | Turan)',
    'ut:"Poisson binomial" & (ratio | curvature | Turan)',
    '(ab:"Bernoulli sum" | ab:"Bernoulli sums") & '
    '(mode | modal | curvature | "adjacent ratio")',
    '(ab:"three consecutive probabilities" | ab:"adjacent probabilities") & '
    '(variance | mode)',
    '(ab:"real-rooted" | ab:"real zeros") & "coefficient ratio" & variance',
    '(ab:"ultra log-concave" | ab:"strong log-concavity") & '
    '(mode | curvature) & variance',
    'rft:"Poisson binomial" & (mode | ratio | curvature)',
]

ARXIV_QUERIES = [
    'all:"poisson binomial"',
    'all:"Poisson binomial" AND (all:mode OR all:modal)',
    'all:"Poisson binomial" AND all:"log-concavity"',
    'all:"Poisson binomial" AND (all:ratio OR all:curvature OR all:Turan)',
    'all:"Bernoulli sums" AND (all:mode OR all:modal OR all:curvature)',
]

OPENALEX_QUERIES = [
    (
        "poisson_binomial_mode",
        'title_and_abstract.search:"Poisson binomial",title_and_abstract.search:mode',
    ),
    (
        "poisson_binomial_mode_variance",
        'title_and_abstract.search:"Poisson binomial",title_and_abstract.search:mode,'
        "title_and_abstract.search:variance",
    ),
    (
        "poisson_binomial_adjacent_ratio",
        'title_and_abstract.search:"Poisson binomial",title_and_abstract.search:"adjacent ratio"',
    ),
    (
        "poisson_binomial_turan",
        'title_and_abstract.search:"Poisson binomial",title_and_abstract.search:Turan',
    ),
    (
        "bernoulli_sum_mode_variance",
        'title_and_abstract.search:"Bernoulli sum",title_and_abstract.search:mode,'
        "title_and_abstract.search:variance",
    ),
]

SEEDS = [
    {
        "label": "Hillion--Johnson",
        "doi": "10.1214/14-AOP973",
        "openalex_id": "W1768331708",
    },
    {
        "label": "Bobkov--Marsiglietti--Melbourne",
        "doi": "10.1017/S096354832100016X",
        "openalex_id": "W3045400294",
    },
    {
        "label": "Darroch",
        "doi": "10.1214/AOMS/1177703287",
        "openalex_id": "W1999942773",
    },
    {
        "label": "Pitman",
        "doi": "10.1006/JCTA.1997.2747",
        "openalex_id": "W1969474871",
    },
    {
        "label": "Duembgen--Wellner",
        "doi": "10.1016/J.SPL.2020.108862",
        "openalex_id": "W3037956058",
    },
    {
        "label": "Baillon--Cominetti--Vaisman",
        "doi": "10.1017/S0963548315000127",
        "openalex_id": "W1570731788",
    },
]

SCREENING_LEDGER = [
    {
        "label": "Hillion--Johnson 2016",
        "identifier": "doi:10.1214/14-AOP973",
        "primary_text_locator": "Theorem A.2 (78), Corollary A.3 (79)",
        "decision": "known_input",
        "reason": "Supplies the two cubic propagation inequalities, but no first-descent mass windows or variance synthesis.",
    },
    {
        "label": "Bobkov--Marsiglietti--Melbourne 2022",
        "identifier": "doi:10.1017/S096354832100016X",
        "primary_text_locator": "Theorem 1.1, Corollary 3.2",
        "decision": "known_input",
        "reason": "Supplies the lattice max-atom/variance inequality, but contains no adjacent-mass curvature bound.",
    },
    {
        "label": "Darroch 1964",
        "identifier": "doi:10.1214/AOMS/1177703287",
        "primary_text_locator": "main mode theorem",
        "decision": "exclude_no_implication",
        "reason": "Locates a mode relative to the mean without quantifying a three-mass descent.",
    },
    {
        "label": "Pitman 1997",
        "identifier": "doi:10.1006/JCTA.1997.2747",
        "primary_text_locator": "equations (20)--(21)",
        "decision": "exclude_no_implication",
        "reason": "Bounds ratios through tilted laws, not neighbouring-ratio curvature normalized by the original variance.",
    },
    {
        "label": "Duembgen--Wellner 2020",
        "identifier": "doi:10.1016/J.SPL.2020.108862",
        "primary_text_locator": "Proposition 2",
        "decision": "exclude_no_implication",
        "reason": "Gives qualitative decrease of an index-weighted density ratio, with no finite variance-scale constant.",
    },
    {
        "label": "Johnson 2017",
        "identifier": "arxiv:1507.06268",
        "primary_text_locator": "Lemma 5.1, Example 5.2",
        "decision": "exclude_no_implication",
        "reason": "The c-log-concavity constant is in a different normalization and can be arbitrarily smaller than 1/V.",
    },
    {
        "label": "Baillon--Cominetti--Vaisman 2016",
        "identifier": "doi:10.1017/S0963548315000127",
        "primary_text_locator": "Theorem 1",
        "decision": "exclude_no_implication",
        "reason": "Sharp one-atom anti-concentration contains no adjacent masses.",
    },
    {
        "label": "Gnedin 2024",
        "identifier": "arxiv:2408.06477",
        "primary_text_locator": "mode and cross-modality results",
        "decision": "exclude_no_implication",
        "reason": "Refines modal location and stability, not local descent magnitude.",
    },
    {
        "label": "Tang--Tang survey",
        "identifier": "doi:10.1214/22-STS852",
        "primary_text_locator": "shape and approximation survey",
        "decision": "exclude_no_overlap",
        "reason": "No finite first-descent curvature theorem or implication was found.",
    },
    {
        "label": "Marsiglietti--Melbourne ULC(n) estimate",
        "identifier": "doi:10.1093/imrn/rnag023",
        "primary_text_locator": "after Definition 1.3",
        "decision": "closest_nonimplication",
        "reason": "The same-three-mass inequality is normalized by summand count; an explicit variance-one family makes its bound tend to zero.",
    },
    {
        "label": "Fontana--Semeraro 2025",
        "identifier": "doi:10.1214/25-ECP741; arxiv:2410.13920",
        "primary_text_locator": "full paper",
        "decision": "exclude_population_mismatch",
        "reason": "Studies fibers of generally dependent Bernoulli representations and does not characterize the independent subclass.",
    },
    {
        "label": "Brazitikos--Pafis 2026",
        "identifier": "arxiv:2601.15444",
        "primary_text_locator": "full paper",
        "decision": "exclude_qualitative_only",
        "reason": "Uses ordinary discrete log-concavity for random-polytope thresholds, with no quantitative modal curvature.",
    },
    {
        "label": "Gaxiola--Melbourne--Pigno--Pollard 2025",
        "identifier": "arxiv:2505.05793",
        "primary_text_locator": "Theorems 1.3/4.4 and 1.4/4.5",
        "decision": "exclude_one_atom_only",
        "reason": "Bounds variance through one selected or maximal atom; no adjacent masses occur.",
    },
    {
        "label": "Jakimiuk--Murawski--Nayar--Slobodzianiuk 2024",
        "identifier": "doi:10.1016/j.disc.2024.114020",
        "primary_text_locator": "Theorems 1--2",
        "decision": "exclude_one_atom_only",
        "reason": "Optimizes atom/variance quantities for log-concave and ULC laws, not three-mass curvature.",
    },
    {
        "label": "Aravinda 2024",
        "identifier": "doi:10.1016/j.disc.2023.113683",
        "primary_text_locator": "Theorem 1.4",
        "decision": "exclude_one_atom_only",
        "reason": "Gives a maximum-atom inequality without adjacent ratios.",
    },
    {
        "label": "Marsiglietti--Melbourne Littlewood--Offord note",
        "identifier": "arxiv:2510.25869",
        "primary_text_locator": "Theorem 1.1, Proposition 3.4",
        "decision": "exclude_one_atom_only",
        "reason": "Gives maximum-atom bounds for weighted and Poisson--binomial sums, not local curvature.",
    },
    {
        "label": "Kontorovich 2026",
        "identifier": "arxiv:2601.04079",
        "primary_text_locator": "full paper",
        "decision": "exclude_no_overlap",
        "reason": "Treats variance, unimodality, maximal atoms, and total variation without neighbouring-ratio estimates.",
    },
    {
        "label": "Broadie--Petkova 2026",
        "identifier": "arxiv:2603.09019",
        "primary_text_locator": "full paper",
        "decision": "exclude_no_overlap",
        "reason": "Treats Poisson-trinomial decomposition, real-rootedness, and mode location without descent magnitude.",
    },
    {
        "label": "Kurauskas 2026",
        "identifier": "arxiv:2603.11043",
        "primary_text_locator": "full paper",
        "decision": "exclude_asymptotic_only",
        "reason": "Gives asymptotic concentration comparisons, not a finite local-curvature inequality.",
    },
    {
        "label": "Kovacevic 2026",
        "identifier": "arxiv:2605.11831",
        "primary_text_locator": "full paper",
        "decision": "exclude_population_or_quantity_mismatch",
        "reason": "Studies entropy of ternary sums and parity-conditioned ULC, not the target three-mass bound.",
    },
]


def request(
    base_url: str,
    params: dict[str, Any] | None = None,
    *,
    attempts: int = 4,
) -> tuple[bytes, dict[str, Any]]:
    url = base_url
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(req, timeout=45) as response:
                body = response.read()
                return body, {
                    "requested_url": url,
                    "final_url": response.geturl(),
                    "http_status": response.status,
                    "content_type": response.headers.get("Content-Type"),
                    "bytes": len(body),
                }
        except urllib.error.HTTPError as exc:
            if exc.code not in (429, 500, 502, 503, 504) or attempt + 1 == attempts:
                raise
            delay = 2 ** (attempt + 1)
            retry_after = exc.headers.get("Retry-After")
            if retry_after and retry_after.isdigit():
                delay = max(delay, int(retry_after))
            time.sleep(delay)
        except urllib.error.URLError:
            if attempt + 1 == attempts:
                raise
            time.sleep(2 ** (attempt + 1))
    raise RuntimeError("unreachable")


def request_json(
    base_url: str,
    params: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    body, transport = request(base_url, params)
    return json.loads(body.decode("utf-8")), transport


def error_record(exc: Exception) -> dict[str, str]:
    return {"error_type": type(exc).__name__, "error": str(exc)}


def title_text(value: Any) -> str | None:
    if isinstance(value, str):
        return value
    if not isinstance(value, dict):
        return None
    parts = [value.get("title"), value.get("subtitle"), value.get("addition")]
    return ": ".join(str(part) for part in parts if part)


def year_cutoff_status(year: Any) -> dict[str, Any]:
    """Mark coarse-year records without pretending they have day-level precision."""
    try:
        numeric_year = int(year)
    except (TypeError, ValueError):
        return {
            "cutoff_status": "unknown_date",
            "eligible_through_cutoff": None,
        }
    if numeric_year < CUTOFF_DATE.year:
        return {
            "cutoff_status": "confirmed_before_cutoff_year",
            "eligible_through_cutoff": True,
        }
    if numeric_year > CUTOFF_DATE.year:
        return {
            "cutoff_status": "confirmed_after_cutoff_year",
            "eligible_through_cutoff": False,
        }
    return {
        "cutoff_status": "cutoff_year_day_unavailable",
        "eligible_through_cutoff": None,
    }


def exact_date_cutoff_status(raw_date: Any) -> dict[str, Any]:
    if not isinstance(raw_date, str) or not raw_date:
        return {"cutoff_status": "unknown_date", "eligible_through_cutoff": None}
    try:
        parsed = date.fromisoformat(raw_date[:10])
    except ValueError:
        return {"cutoff_status": "unparsed_date", "eligible_through_cutoff": None}
    eligible = parsed <= CUTOFF_DATE
    return {
        "cutoff_status": "confirmed_on_or_before_cutoff" if eligible else "confirmed_after_cutoff",
        "eligible_through_cutoff": eligible,
    }


def zbmath_zero_hit(exc: urllib.error.HTTPError) -> dict[str, Any] | None:
    """Accept a zbMATH 404 as zero hits only when its JSON status says so."""
    if exc.code != 404:
        return None
    raw = exc.read()
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return None
    status = payload.get("status") or {}
    internal = str(status.get("internal_code", ""))
    if status.get("execution_bool") is not False or not any(
        marker in internal.lower() for marker in ("zero results", "no results")
    ):
        return None
    return {
        "http_status": 404,
        "total": 0,
        "returned": 0,
        "results": [],
        "zbmath_status": status,
        "interpretation": "Verified zbMATH JSON status reports a successful zero-result query.",
    }


def zbmath_searches() -> list[dict[str, Any]]:
    searches = []
    endpoint = "https://api.zbmath.org/v1/document/_structured_search"
    for query in ZBMATH_QUERIES:
        params = {"Anywhere": query, "page": 0, "results_per_page": 100}
        record: dict[str, Any] = {"field": "Anywhere", "query": query, "params": params}
        try:
            payload, transport = request_json(endpoint, params)
            status = payload.get("status", {})
            record.update(
                {
                    "transport": transport,
                    "total": status.get("nr_total_results"),
                    "returned": status.get("nr_request_results"),
                    "results": [
                        {
                            "id": item.get("id"),
                            "title": title_text(item.get("title")),
                            "year": item.get("year"),
                            "document_type": item.get("document_type"),
                            "zbmath_url": item.get("zbmath_url"),
                            "identifier": item.get("identifier"),
                            **year_cutoff_status(item.get("year")),
                        }
                        for item in payload.get("result", [])
                    ],
                }
            )
        except urllib.error.HTTPError as exc:
            zero_hit = zbmath_zero_hit(exc)
            if zero_hit is not None:
                record.update(zero_hit)
            else:
                record.update(error_record(exc))
        except Exception as exc:  # preserve partial audit if one service fails
            record.update(error_record(exc))
        searches.append(record)
        time.sleep(0.15)

    endpoint = "https://api.zbmath.org/v1/document/_search"
    for query in ZBMATH_FREE_QUERIES:
        params = {"search_string": query, "page": 0, "results_per_page": 100}
        record = {"field": "free syntax", "query": query, "params": params}
        try:
            payload, transport = request_json(endpoint, params)
            status = payload.get("status", {})
            record.update(
                {
                    "transport": transport,
                    "total": status.get("nr_total_results"),
                    "returned": status.get("nr_request_results"),
                    "results": [
                        {
                            "id": item.get("id"),
                            "title": title_text(item.get("title")),
                            "year": item.get("year"),
                            "document_type": item.get("document_type"),
                            "zbmath_url": item.get("zbmath_url"),
                            "identifier": item.get("identifier"),
                            **year_cutoff_status(item.get("year")),
                        }
                        for item in payload.get("result", [])
                    ],
                }
            )
        except urllib.error.HTTPError as exc:
            zero_hit = zbmath_zero_hit(exc)
            if zero_hit is not None:
                record.update(zero_hit)
            else:
                record.update(error_record(exc))
        except Exception as exc:
            record.update(error_record(exc))
        searches.append(record)
        time.sleep(0.15)
    return searches


def arxiv_searches() -> list[dict[str, Any]]:
    searches = []
    endpoint = "https://export.arxiv.org/api/query"
    atom = {"a": "http://www.w3.org/2005/Atom", "o": "http://a9.com/-/spec/opensearch/1.1/"}
    for query in ARXIV_QUERIES:
        cutoff_query = (
            f"({query}) AND submittedDate:[199101010000 TO {CUTOFF_ARXIV}]"
        )
        params = {
            "search_query": cutoff_query,
            "start": 0,
            "max_results": 100,
            "sortBy": "submittedDate",
            "sortOrder": "descending",
        }
        record: dict[str, Any] = {
            "query": query,
            "cutoff_query": cutoff_query,
            "params": params,
        }
        try:
            body, transport = request(endpoint, params)
            root = ET.fromstring(body)
            total_node = root.find("o:totalResults", atom)
            entries = []
            for entry in root.findall("a:entry", atom):
                entries.append(
                    {
                        "id": (entry.findtext("a:id", default="", namespaces=atom)),
                        "title": " ".join(
                            entry.findtext("a:title", default="", namespaces=atom).split()
                        ),
                        "published": entry.findtext("a:published", default="", namespaces=atom),
                        "updated": entry.findtext("a:updated", default="", namespaces=atom),
                        **exact_date_cutoff_status(
                            entry.findtext("a:published", default="", namespaces=atom)
                        ),
                    }
                )
            record.update(
                {
                    "transport": transport,
                    "total": int(total_node.text) if total_node is not None else None,
                    "returned": len(entries),
                    "results": entries,
                }
            )
        except Exception as exc:
            record.update(error_record(exc))
        searches.append(record)
        time.sleep(3.1)
    return searches


def openalex_summary(item: dict[str, Any]) -> dict[str, Any]:
    publication_date = item.get("publication_date")
    return {
        "id": item.get("id"),
        "doi": item.get("doi"),
        "title": item.get("display_name") or item.get("title"),
        "year": item.get("publication_year"),
        "publication_date": publication_date,
        "type": item.get("type"),
        "cited_by_count": item.get("cited_by_count"),
        **exact_date_cutoff_status(publication_date),
    }


def openalex_searches() -> dict[str, Any]:
    endpoint = "https://api.openalex.org/works"
    discovery = []
    for label, query_filter in OPENALEX_QUERIES:
        cutoff_filter = f"{query_filter},to_publication_date:{CUTOFF}"
        params = {
            "filter": cutoff_filter,
            "per-page": 100,
            "select": "id,doi,display_name,publication_year,publication_date,type,cited_by_count",
        }
        record: dict[str, Any] = {
            "label": label,
            "filter": query_filter,
            "cutoff_filter": cutoff_filter,
            "params": params,
        }
        try:
            payload, transport = request_json(endpoint, params)
            record.update(
                {
                    "transport": transport,
                    "total": payload.get("meta", {}).get("count"),
                    "returned": len(payload.get("results", [])),
                    "results": [openalex_summary(item) for item in payload.get("results", [])],
                }
            )
        except Exception as exc:
            record.update(error_record(exc))
        discovery.append(record)
        time.sleep(0.15)

    citations = []
    for seed in SEEDS:
        params = {
            "filter": f"cites:{seed['openalex_id']},to_publication_date:{CUTOFF}",
            "per-page": 200,
            "select": "id,doi,display_name,publication_year,publication_date,type,cited_by_count",
        }
        record = {"seed": seed, "params": params}
        try:
            payload, transport = request_json(endpoint, params)
            record.update(
                {
                    "transport": transport,
                    "total": payload.get("meta", {}).get("count"),
                    "returned": len(payload.get("results", [])),
                    "results": [openalex_summary(item) for item in payload.get("results", [])],
                }
            )
        except Exception as exc:
            record.update(error_record(exc))
        citations.append(record)
        time.sleep(0.15)
    return {"discovery_searches": discovery, "forward_citations": citations}


def semantic_scholar_citations() -> list[dict[str, Any]]:
    records = []
    fields = "title,year,externalIds,citationCount"
    for seed in SEEDS:
        encoded_doi = urllib.parse.quote("DOI:" + seed["doi"], safe=":")
        endpoint = f"https://api.semanticscholar.org/graph/v1/paper/{encoded_doi}/citations"
        record: dict[str, Any] = {"seed": seed, "page_size": 100}
        citing = []
        transports = []
        next_offset: int | None = 0
        pages = 0
        while next_offset is not None and pages < 10:
            params = {"limit": 100, "offset": next_offset, "fields": fields}
            try:
                payload, transport = request_json(endpoint, params)
            except Exception as exc:
                record.update(error_record(exc))
                break
            transports.append(transport)
            pages += 1
            for edge in payload.get("data", []):
                paper = edge.get("citingPaper") or {}
                citing.append(
                    {
                        "paperId": paper.get("paperId"),
                        "title": paper.get("title"),
                        "year": paper.get("year"),
                        "externalIds": paper.get("externalIds"),
                        "citationCount": paper.get("citationCount"),
                        **year_cutoff_status(paper.get("year")),
                    }
                )
            next_offset = payload.get("next")
            time.sleep(1.1)
        record.update(
            {
                "transports": transports,
                "pages": pages,
                "next": next_offset,
                "complete": next_offset is None,
                "returned": len(citing),
                "results": citing,
            }
        )
        records.append(record)
    return records


def crossref_seeds() -> list[dict[str, Any]]:
    records = []
    for seed in SEEDS:
        endpoint = "https://api.crossref.org/works/" + urllib.parse.quote(seed["doi"], safe="")
        record: dict[str, Any] = {"seed": seed}
        try:
            payload, transport = request_json(endpoint)
            item = payload.get("message", {})
            record.update(
                {
                    "transport": transport,
                    "title": (item.get("title") or [None])[0],
                    "type": item.get("type"),
                    "published": item.get("published"),
                    "is_referenced_by_count": item.get("is-referenced-by-count"),
                    "reference_count": item.get("reference-count"),
                    "URL": item.get("URL"),
                }
            )
        except Exception as exc:
            record.update(error_record(exc))
        records.append(record)
        time.sleep(0.15)
    return records


def mathscinet_probe() -> dict[str, Any]:
    endpoint = "https://mathscinet.ams.org/mathscinet/api/publications/search"
    params = {
        "query": 'ti:"Poisson binomial"',
        "currentPage": 1,
        "pageSize": 20,
        "sort": "relevance",
    }
    record: dict[str, Any] = {"params": params}
    try:
        body, transport = request(endpoint, params)
        prefix = body[:1000].decode("utf-8", errors="replace")
        title = None
        marker = "<title>"
        if marker in prefix and "</title>" in prefix:
            title = prefix.split(marker, 1)[1].split("</title>", 1)[0].strip()
        record.update(
            {
                "transport": transport,
                "page_title": title,
                "result_payload_available": transport.get("content_type", "").startswith(
                    "application/json"
                ),
                "interpretation": (
                    "Search redirected to institutional WAYF; no MathSciNet result list or review text was accessible."
                    if "connect.liblynx.com/wayf/" in transport.get("final_url", "")
                    else "Inspect transport and content type manually."
                ),
            }
        )
    except Exception as exc:
        record.update(error_record(exc))
    return record


def build_audit() -> dict[str, Any]:
    generated_at = datetime.now(timezone.utc)
    return {
        "schema_version": 1,
        "generated_at": generated_at.isoformat(),
        "literature_cutoff": CUTOFF,
        "cutoff_policy": {
            "arxiv": "server-side submittedDate upper bound through 2026-07-16 23:59",
            "openalex": "server-side to_publication_date:2026-07-16 on discovery and citations",
            "zbmath_open": (
                "result metadata exposes year only; pre/post-cutoff years are marked, and "
                "cutoff-year records are explicitly day-ambiguous"
            ),
            "semantic_scholar": (
                "citation metadata exposes year only; pre/post-cutoff years are marked, and "
                "cutoff-year records are explicitly day-ambiguous"
            ),
            "crossref": "fixed seed metadata only; no discovery result set",
            "same_day_snapshot": generated_at.date() == CUTOFF_DATE,
        },
        "target": {
            "population": "finite nondegenerate Poisson-binomial laws",
            "condition": "V >= 1 and first strict descent D is nonterminal",
            "inequality": "V * (1 - f[D-1] * f[D+1] / f[D]^2) >= 1/4",
        },
        "method_note": (
            "These are discovery results and citation-graph metadata. Primary full text, not snippets, "
            "controls theorem-overlap decisions. Hit counts are database-specific and may change."
        ),
        "screening_ledger": {
            "scope": (
                "Compact ledger for every mathematically plausible theorem overlap promoted to "
                "primary-full-text review; bulk irrelevant title screening remains represented by "
                "the query result lists."
            ),
            "screened_on": CUTOFF,
            "decisions": SCREENING_LEDGER,
        },
        "zbmath_open": {"searches": zbmath_searches()},
        "arxiv": {"searches": arxiv_searches()},
        "openalex": openalex_searches(),
        "semantic_scholar": {"forward_citations": semantic_scholar_citations()},
        "crossref": {"seed_metadata": crossref_seeds()},
        "mathscinet": mathscinet_probe(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/poisson_binomial_novelty_queries_20260716.json"),
    )
    args = parser.parse_args()
    audit = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
