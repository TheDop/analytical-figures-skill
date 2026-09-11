# VISION.md — the ideal analytical-chemistry agent (product ideas, NOT skill features)

A planning doc, deliberately separate from `ROADMAP.md`. **ROADMAP = features buildable
in THIS skill now** (reference docs, dataclasses, gates, generators). **VISION = the ideal
analytical-chemistry *agent/product*** — ideas that need platform infrastructure (a memory
substrate, connectors, a persistent problem-state store, agent orchestration) we are *not*
building here. Kept so the thinking isn't lost, and so we don't fool ourselves that these
are a weekend's work in the skill.

Source: the "ideal analytical-chemistry agentic product" thought exercise (2026-07-03),
prompted by Anthropic's Claude Science launch. The exercise walked turn-by-turn through a
method-development project (an IND: quantify API + excipients + identify all impurities) and
kept asking *what capability does this turn demand?* and *is it buildable in the skill?*

## The turn-one thesis

A great analytical agent's *first* move is not "design a method" — it's to convert a wish
("identify **all** impurities") into a **risk-managed, phase-appropriate spec** by
decomposing the goal, characterising the analyte, gating against the lab's real
constraints, and asking only the questions that fork the plan. It earns the right to
propose by first proving it understands why the goal as written can't be met. The value of
a method is knowing its limits — so the agent's job is to state them.

## The four primitives (+ one cross-cutting property)

Three were named readily; the fourth is the load-bearing one that a naive product flattens.

1. **Problem model / decomposition** — one goal → the coupled sub-methods, *with the
   couplings* (the impurity method's diluent must keep a labile API alive in the autosampler,
   or you manufacture your own degradants). Dependency-graph reasoning; "flexible" = the
   ability to *backtrack* when a downstream stage invalidates an upstream choice.
2. **Constraint model as an admissibility gate** — instruments/budget/phase are *typed,
   first-class* and actively **reject** proposals, and can declare the **goal itself
   infeasible** ("identify all impurities to 0.05% — no MS here → not reachable; control at
   source"). Instrument constraints are `{capability, throughput, lead-time/scheduling,
   cost-to-access}`, not booleans. Lab/researcher profile is cross-project; goal state is
   project-scoped.
3. **Knowledge connection** — pull molecular/impurity/regulatory facts from databases +
   literature, each carrying **provenance + confidence** (experimental vs predicted; monograph
   vs vendor PDF). Prices/availability via *live search*, never recall.
4. **★ The epistemic layer** (the missing primitive) — the agent's model of *its own
   knowledge state*. Three kinds of not-knowing, routed differently: *knowable by lookup* →
   fetch; *knowable only from the user* → ask; *unknowable/unprovable* (completeness of "all
   impurities") → design around as managed risk. Plus **value-of-information questioning**
   (ask only the unknowns that fork the plan) and a **scoping stop-criterion** (done when no
   remaining unknown would change the approach beyond tolerance).

**Cross-cutting: a legible, user-correctable, traceable problem-state artifact.** The
researcher owns the regulatory accountability, so must *see and edit* the agent's model of the
problem; every fact carries provenance + confidence; everything traces to the goal so a
changed fact (mesylate not HCl; starch not lactose) **propagates and flags what's now
invalidated**. This is what makes the pipeline flexible rather than a brittle recipe.

## Why these are parked (not in ROADMAP)

Each needs infrastructure outside a figure/analysis skill: a persistent memory store (2, the
profile/goal split), live DB/literature connectors (3), a dependency graph with invalidation
(1, the cross-cutting property), and an uncertainty-tracking state object (4). That's an
*agent platform*, not a reference doc. We explicitly are not building it — but the exercise
was generative, and if a real platform ever gets built, this is the spec.

## What DID land in the skill (the achievable slivers)

The exercise's payoff was a handful of **form** artifacts shaped like what the skill already
does — logged/built in `ROADMAP.md`:

- **`references/cocrystal_id.md`** (BUILT 2026-07-03) — primitive #2's admissibility gate +
  #1's decision layer, specialised to cocrystal-vs-salt-vs-mixture ID. The
  method-design sibling of `figure_selection.md`.
- **Provenance stamp + number-fidelity critic** (BUILT, earlier 2026-07-03) — primitive #3's
  "provenance + confidence" and the cross-cutting traceability, at the figure/script level.
- **Web-search-then-cite discipline** (SKILL.md) — primitive #3's "facts looked up, not
  recalled."
- **Parked as skill-features** (see ROADMAP): a general capability→requirement matrix, an
  excipient→degradation-liability map, a forced-degradation design generator (sibling of
  `doe.py`), a cited procurement/budget ledger.

The guarded-**exploratory-data-analysis** marquee (actor-critic EDA that proposes findings
which the gates then dispose of) lives in `ROADMAP.md` — it straddles vision and skill: the
orchestration is platform, but the *gates that make exploration safe* are already ours.
