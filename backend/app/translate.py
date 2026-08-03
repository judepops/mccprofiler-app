"""Turn a question into a query. The only part of the app that uses a model.

Retrieval lives in `query.py` and is plain pandas, so the provider here is
swappable without touching anything that produces genes. That was the point of
the split: a different model, or no model at all, changes how a query is
authored and never changes what a query returns.

Two providers, selected by `LLM_PROVIDER` in `.env`:

    gemini   (default) free tier, needs GEMINI_API_KEY
    claude   needs ANTHROPIC_API_KEY

The schema below is deliberately flat rather than a discriminated union. A
union expressed with `anyOf` is fine for Claude but sits outside the OpenAPI
subset Gemini's structured output accepts, so one flat filter shape with a
`type` discriminator and optional fields works on both. Validation of the
combination happens in `query.py`, where it belongs anyway.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

DEFAULT_PROVIDER = "gemini"

# Overridable because model ids move faster than this file does. Checked
# against `client.models.list()` on 2026-08-03: the 2.0 family still lists but
# its free tier is now zero-quota, which surfaces as a 429 reading
# `limit: 0` rather than as a 404, so a working key looks broken. If the
# default starts failing, list the models the key actually serves rather than
# guessing an id.
DEFAULT_MODELS = {
    "gemini": os.environ.get("GEMINI_MODEL", "gemini-3.5-flash-lite"),
    "claude": os.environ.get("ANTHROPIC_MODEL", "claude-opus-5"),
}


class TranslationUnavailable(RuntimeError):
    """No usable credentials or SDK. The caller degrades to manual queries."""


class TranslationFailed(RuntimeError):
    """The provider was reachable but did not return a usable query."""


# ---------------------------------------------------------------------------
# schema
# ---------------------------------------------------------------------------


def flat_schema(axes: list[str], cohorts: list[str], groups: list[str],
                features: list[str]) -> dict:
    """One filter shape, portable across providers.

    `type` says which fields matter; the rest are optional and ignored when
    irrelevant. Less elegant than a union, but it survives the trip through
    two different structured-output implementations.

    Pass `cohorts=[]` to remove external reference sets from the query language
    entirely, which is how the ask endpoint calls it. Selection then runs on
    contact architecture alone and external membership survives only as
    annotation on the results. Dropping the field from the schema rather than
    filtering afterwards is deliberate: a model cannot emit a filter it has no
    vocabulary for.
    """
    types = ["axis", "archetype", "feature"] + (["cohort"] if cohorts else [])
    props: dict[str, Any] = {
        "type": {"type": "string", "enum": types},
        "axis": {"type": "string", "enum": axes},
        "group": {"type": "string", "enum": groups},
    }
    if cohorts:
        props["cohort"] = {"type": "string", "enum": cohorts}
    return {
        "type": "object",
        "properties": {
            "filters": {
                "type": "array",
                "description": "Conditions combined with AND. Use as few as possible.",
                "items": {
                    "type": "object",
                    "properties": {
                        **props,
                        "feature": {"type": "string", "enum": features},
                        "direction": {"type": "string", "enum": ["high", "low"]},
                        "percentile": {
                            "type": "number",
                            "description": "Cutoff, default 25.",
                        },
                    },
                    "required": ["type"],
                },
            },
            "interpretation": {
                "type": "string",
                "description": "One plain sentence restating the query, naming each "
                               "direction explicitly so the user can spot a misreading.",
            },
            "unsupported": {
                "type": "string",
                "description": "Any part of the question this dataset cannot answer, "
                               "such as TADs, compartments, insulation, expression "
                               "level, or membership of an external gene set. Empty "
                               "if the question is fully answerable. Filters should "
                               "still cover the answerable part.",
            },
            "reasoning": {
                "type": "array",
                "description": "Your working, one short step per element. Say which "
                               "term you mapped, what you mapped it to, and why that "
                               "pole or feature. Mention anything you considered and "
                               "rejected. Three to six steps.",
                "items": {"type": "string"},
            },
        },
        "required": ["filters", "interpretation", "reasoning"],
    }


# `CONTEXT.md` at the repo root. Loaded once and cached, because it is static
# and the prompt is rebuilt on every question.
_CONTEXT_PATH = Path(__file__).resolve().parents[2] / "CONTEXT.md"
_context_cache: str | None = None


def load_context() -> str:
    """The literature glossary, or an empty string if it is missing.

    Missing is survivable and must not be fatal: without it the model still
    translates, it just no longer knows that "housekeeping" is a cohort rather
    than the similarly named group. Absence is reported by `status()` so the
    degradation is visible rather than silent.
    """
    global _context_cache
    if _context_cache is None:
        try:
            _context_cache = _CONTEXT_PATH.read_text(encoding="utf-8")
        except OSError:
            _context_cache = ""
    return _context_cache


def build_prompt(pole_lines: str, feature_reference: str = "") -> str:
    """The system instruction. Identical across providers on purpose."""
    context = load_context()
    glossary = (
        "\n\nREFERENCE GLOSSARY. This tells you what a user's term means in the "
        "literature and what it maps to in this dataset. Two things to take from "
        "it. Where it says a term maps to nothing, do NOT invent a filter for it, "
        "set `unsupported` instead. Where it warns that an obvious-looking mapping "
        "is wrong, such as housekeeping being a cohort and not the similarly named "
        "group, follow the warning over the resemblance.\n\n"
        f"{context}\n"
    ) if context else ""
    return (
        "You translate a biologist's question about gene contact architecture into a "
        "structured query. You do NOT answer the question or name genes. You only "
        "build the query; a deterministic filter runs it.\n\n"
        f"Axes and which end is which:\n{pole_lines}\n\n"
        "Rules:\n"
        "- Use the fewest filters that capture the question. Extra filters silently "
        "shrink the result to a handful.\n"
        "- Map direction words to the correct POLE, not to 'high'. 'Long-range' is the "
        "LOW end of PC2, not the high end.\n"
        "- PC1 is peak COUNT and reach, not total signal. It correlates with "
        "total_mcc at only +0.03. If a question is about how much signal a gene has, "
        "use the total_mcc feature directly; PC1 will not answer it.\n"
        "- Every axis name describes only part of its axis, between 22 and 35 percent "
        "of the loading mass. Treat a name as the dominant contrast, not a "
        "definition, and prefer a named FEATURE over an axis when the question names "
        "something a feature measures directly.\n"
        "- Selection runs on contact architecture ONLY. You have no filter for "
        "external gene sets, and there is deliberately no vocabulary for one. If a "
        "question asks for a biological category that is not architectural, such as "
        "immune, housekeeping, essential, conserved or expressed, that part is "
        "`unsupported`. Do NOT approximate it with an axis or a group: the point of "
        "this tool is that architecture is measured independently of borrowed "
        "labels, and selecting on those labels would destroy that independence. "
        "Membership is reported alongside the results afterwards instead.\n"
        "- Each filter object only needs the fields its type uses: axis filters need "
        "`axis` and `direction`; cohort filters need `cohort`; archetype filters need "
        "`group`; feature filters need `feature` and `direction`.\n"
        "- Always fill `interpretation` with one plain sentence naming each direction "
        "explicitly, so the user can spot a misreading.\n"
        "- If part of the question asks for something this dataset does not hold, "
        "put that in `unsupported`, in plain words, and build filters only for the "
        "part that IS answerable. An approximation offered silently is worse than "
        "a stated gap.\n"
        + ("\n\n" + feature_reference if feature_reference else "")
        + glossary +
        "\nReturn only the JSON object."
    )


# ---------------------------------------------------------------------------
# providers
# ---------------------------------------------------------------------------


def _translate_gemini(question: str, schema: dict, system: str, model: str) -> dict:
    key = os.environ.get("GEMINI_API_KEY") or os.environ.get("GOOGLE_API_KEY")
    if not key:
        raise TranslationUnavailable(
            "No Gemini credentials. Copy .env.example to .env, add GEMINI_API_KEY "
            "(free from aistudio.google.com/apikey), and restart the backend."
        )
    try:
        from google import genai
        from google.genai import types
    except ImportError:
        raise TranslationUnavailable(
            "The google-genai SDK is not installed. Run `pip install google-genai` "
            "in the mccapp env."
        ) from None

    try:
        client = genai.Client(api_key=key)
        resp = client.models.generate_content(
            model=model,
            contents=question,
            config=types.GenerateContentConfig(
                system_instruction=system,
                response_mime_type="application/json",
                response_schema=schema,
                temperature=0,  # translation should be reproducible
            ),
        )
        return json.loads(resp.text)
    except Exception as e:  # noqa: BLE001
        msg = str(e)
        # A quota refusal arrives as several hundred characters of nested JSON
        # and reads like a credentials problem, which it is not. Say which of
        # the two it actually is.
        if "RESOURCE_EXHAUSTED" in msg or "429" in msg:
            kind = ("has no free-tier quota for it" if "limit: 0" in msg
                    else "is rate limited on it")
            raise TranslationFailed(
                f"model {model!r} refused the request: this key {kind}. "
                "Set GEMINI_MODEL in .env to a model the key serves, or wait "
                "and retry."
            ) from None
        raise TranslationFailed(msg) from None


def _translate_claude(question: str, schema: dict, system: str, model: str) -> dict:
    if not (os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN")):
        raise TranslationUnavailable(
            "No Anthropic credentials. Copy .env.example to .env, add "
            "ANTHROPIC_API_KEY, and restart the backend."
        )
    try:
        import anthropic
    except ImportError:
        raise TranslationUnavailable(
            "The anthropic SDK is not installed. Run `pip install anthropic`."
        ) from None

    try:
        client = anthropic.Anthropic()
        resp = client.messages.create(
            model=model,
            max_tokens=2000,
            system=system,
            output_config={
                "effort": "low",  # translation, not reasoning
                "format": {"type": "json_schema", "schema": schema},
            },
            messages=[{"role": "user", "content": question}],
        )
        if resp.stop_reason == "refusal":
            raise TranslationFailed("the question was declined by safety classifiers")
        text = next((b.text for b in resp.content if b.type == "text"), None)
        if not text:
            raise TranslationFailed("no query returned")
        return json.loads(text)
    except TranslationFailed:
        raise
    except Exception as e:  # noqa: BLE001
        raise TranslationFailed(str(e)) from None


_PROVIDERS = {"gemini": _translate_gemini, "claude": _translate_claude}


def active_provider() -> str:
    return os.environ.get("LLM_PROVIDER", DEFAULT_PROVIDER).strip().lower()


def status() -> dict[str, Any]:
    """Whether translation is usable, for the UI to show before a question."""
    provider = active_provider()
    keys = {
        "gemini": bool(os.environ.get("GEMINI_API_KEY") or os.environ.get("GOOGLE_API_KEY")),
        "claude": bool(os.environ.get("ANTHROPIC_API_KEY")
                       or os.environ.get("ANTHROPIC_AUTH_TOKEN")),
    }
    ctx = load_context()
    return {
        "provider": provider,
        "model": DEFAULT_MODELS.get(provider),
        "available": keys.get(provider, False),
        "credentials_present": keys,
        # Reported so a missing glossary is visible. Without it the model still
        # answers, it just loses the mappings and the traps, which is a quiet
        # quality drop rather than an error.
        "context_loaded": bool(ctx),
        "context_chars": len(ctx),
        "note": "Only the plain-English box needs a model. Building and running "
                "queries by hand works regardless.",
    }


def translate(question: str, schema: dict, pole_lines: str,
              feature_reference: str = "") -> dict:
    provider = active_provider()
    fn = _PROVIDERS.get(provider)
    if fn is None:
        raise TranslationUnavailable(
            f"unknown LLM_PROVIDER {provider!r}; expected one of "
            f"{', '.join(_PROVIDERS)}"
        )
    return fn(question, schema, build_prompt(pole_lines, feature_reference),
              DEFAULT_MODELS[provider])
