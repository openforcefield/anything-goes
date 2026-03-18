# Copilot Instructions for `anything-goes`

This repository collects community examples, demos, and scripts built on the [Open Force Field](https://openforcefield.org/) software stack. Contributions are low-friction by design, but a reviewer (human or AI) should be able to understand **what** a contribution does and **why** it exists. Use the checklist below to assess any incoming contribution.

---

## Contribution Review Checklist

When evaluating a pull request or new directory, check that the contribution satisfies the criteria below. Flag items that are missing or unclear.

### 1. Purpose is stated

> *Can someone reading the directory name and its README (or notebook title/intro cell) understand the goal within 30 seconds?*

- [ ] There is a brief description somewhere (README, notebook intro cell, or inline comments) of **what problem this solves** or **what concept it demonstrates**.
- [ ] The description explains *why* someone would want to do this, not just *how* the code works.

### 2. Reproducibility is achievable

> *Could someone re-run this with reasonable effort?*

- [ ] An `env.yaml`, `conda-lock.yml` **or** the required packages are obvious from the notebook/script imports.
- [ ] Input files referenced in the code are present in the directory (or fetched from a URL that is documented).
- [ ] No hard-coded absolute paths that only work on the author's machine.

### 3. Scope is self-contained

> *Does the contribution stand alone within its directory?*

- [ ] All files belonging to the contribution are inside a single subdirectory (or are clearly labelled root-level notebooks).
- [ ] The directory does not silently depend on files in sibling directories without documentation.

### 4. README or table entry is present or updated

> *Is the summary table in [README.md](../README.md) up to date?*

- [ ] A row for the new directory exists in the root `README.md` summary table.
- [ ] The description is ≤ 30 words and mentions the primary tool(s) and the outcome.

---

## What does NOT need to be perfect

This is an `anything-goes` repo. The following are explicitly **not required**:

- Unit tests or CI integration.
- Support for multiple operating systems or Python versions.
- Long-form documentation.
- Guaranteed reproducibility forever (dependencies drift; that is accepted).
- Code review for style or best practices (unless the author requests it).

---

## Example of a well-described contribution

> **`optimize-dimer/`** — Loads a bromobenzene–water dimer from two SDF files, parameterizes with Sage 2.0 via Interchange, energy-minimizes with OpenMM, and overlays initial vs. final coordinates. Includes `README.md`, input SDF files, and a self-contained notebook.

This example passes all checklist items: purpose is clear, packages named, inputs present, self-contained directory, and a table row exists.

---

## Flagging a contribution

If a contribution is missing critical context, leave a PR comment such as:

```
This contribution is missing:
- [ ] A brief description of what it demonstrates or solves.
- [ ] An environment file or package version list.
Please add these before merging so others can benefit from your work.
```

You do **not** need to block the PR — the goal is to nudge contributors toward adding just enough information for others to understand and reuse their work.
