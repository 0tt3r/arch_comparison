# Earlier notes

The current complete narrative is [../main.tex](../main.tex), built as
[../main.pdf](../main.pdf). These three earlier notes retain their original
TeX, bibliography dependencies, and PDFs. Their historical assumptions and
claims have not been rewritten to match later results.

- `analytical_forms`: initial analytical laws and proxy curves.
- `nisq_refinement`: native compilation, exact cones, and factory reassessment.
- `refinement_v3`: parallel compilation, dynamics calibration, cultivation, and STAR retries.
- `chronological_2026-09-23`: the complete manuscript before its restructuring into a setup-first
  narrative (newest results first, earlier stages following). Its source uses the main
  `references.bib` and `figures/`; the PDF is the build of commit b4b9231.

From the repository root, `python python/build_manuscript.py --archive refinement_v3`
rebuilds a historical note in `build/tex/archive/`. It supplies the old relative
figure paths without replacing the archived source or PDF.
