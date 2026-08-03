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
from typing import Any

DEFAULT_PROVIDER = "gemini"

# Overridable because model ids move faster than this file does.
DEFAULT_MODELS = {
    "gemini": os.environ.get("GEMINI_MODEL", "gemini-2.0-flash"),
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
    """
    return {
        "type": "object",
        "properties": {
            "filters": {
                "type": "array",
                "description": "Conditions combined with AND. Use as few as possible.",
                "items": {
                    "type": "object",
                    "properties": {
                        "type": {
                            "type": "string",
                            "enum": ["axis", "cohort", "archetype", "feature"],
                        },
                        "axis": {"type": "string", "enum": axes},
                        "cohort": {"type": "string", "enum": cohorts},
                        "group": {"type": "string", "enum": groups},
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
        },
        "required": ["filters", "interpretation"],
    }


def build_prompt(pole_lines: str) -> str:
    """The system instruction. Identical across providers on purpose."""
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
        "- PC1 is contact amount, not biology. Only use it if the question is about "
        "how much signal a gene has.\n"
        "- Prefer a cohort filter for biological categories such as immune, "
        "housekeeping or essential, rather than expressing them as axis positions.\n"
        "- Each filter object only needs the fields its type uses: axis filters need "
        "`axis` and `direction`; cohort filters need `cohort`; archetype filters need "
        "`group`; feature filters need `feature` and `direction`.\n"
        "- Always fill `interpretation` with one plain sentence naming each direction "
        "explicitly, so the user can spot a misreading.\n"
        "Return only the JSON object."
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
        raise TranslationFailed(str(e)) from None


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
    return {
        "provider": provider,
        "model": DEFAULT_MODELS.get(provider),
        "available": keys.get(provider, False),
        "credentials_present": keys,
        "note": "Only the plain-English box needs a model. Building and running "
                "queries by hand works regardless.",
    }


def translate(question: str, schema: dict, pole_lines: str) -> dict:
    provider = active_provider()
    fn = _PROVIDERS.get(provider)
    if fn is None:
        raise TranslationUnavailable(
            f"unknown LLM_PROVIDER {provider!r}; expected one of "
            f"{', '.join(_PROVIDERS)}"
        )
    return fn(question, schema, build_prompt(pole_lines), DEFAULT_MODELS[provider])
