# User preferences

Recorded from work with Jane Rigby, 2026-11-24. Read before editing this project.

- Ask the user when uncertain; do not guess at their intent or consequential scientific choices.
- Do not make anything up. Distinguish verified facts from assumptions, hypotheses, and unknowns.
- Be succinct. Avoid excessive verbosity, generic AI filler, and canned language.
- End-users are scientists. State the code's purpose near the top of each documentation file, then explain its scientific use and relevant limitations.
- **Check with the user before changing file formats that may affect downstream work.** Preserve existing header markers, separators, keyword names, column names/order, units, and naming conventions unless the user explicitly authorizes the specific change. If an instruction is ambiguous, describe the existing format and ask what should change. For example, a request to remove verbose CSV header text does not authorize changing established `##` metadata lines to `#`.
- Add the requested measurements without changing existing extraction calculations or spectral values unless asked.
- Keep scientific documentation concise. Explain why a feature is useful, as well as what it does. State provisional criteria and their need for further testing.
- Avoid verbose explanatory text and unrequested status/reason keywords in data-file headers. Put diagnostic details in reports or documentation.
- Keep plot labels concise and readable. Avoid redundant labels and annotations that interfere with legends or data.
- Ask before adding unrequested comparison data or overlays, such as L2 exposures alongside L3 profiles.
- Preserve user-edited selection files. Create starter files only if absent; never overwrite previous choices.
- When asked to “just tell me,” answer without modifying files.

These preferences do not authorize restoring or changing any existing file format. Such a change still requires the user's approval.
