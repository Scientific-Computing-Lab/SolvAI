# Journal manuscript source

`main.tex` builds the complete reading copy, including the Supplementary Information. For a journal portal that requests a separately paginated supplement, compile the main document first and then the submission wrapper:

```bash
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary_submission.tex
```

The second build reads `main.aux` for two references to main figures. Its content is a submission-only copy of the Supplementary Information with those references made explicit. The scientific text, tables and figures are otherwise the same as in `main.tex`.

The public experimental records and Supplementary Data 1–6 are linked from the root README. The retrospective source-surrogate challenge is disclosed in the manuscript as an aggregate-only release until its fit-level records and code are deposited.
