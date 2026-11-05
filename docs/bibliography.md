# Bibliography conventions

## Setup
- **System:** natbib (author-year, round brackets) with BibTeX and the `plainnat` style. Both are loaded in `structecon.sty`.
- **Master bibliography:** one consolidated `references.bib` at the repo root, printed in the back matter.
- **Building:** `make book` runs BibTeX. Single-chapter builds (`make chNN`) print citations as `[?]`; that is expected.

## Citation commands
| Text form | Command |
|---|---|
| "Golub and Jackson (2010) show …" | `\citet{golub2010naive}` |
| "… as shown before (Golub and Jackson 2010)" | `\citep{golub2010naive}` |
| With a page | `\citet[p.~93]{foucault1978history}`, `\citep[ch.~4]{furtado1961desenvolvimento}` |
| Several works | `\citep{a2001x, b2005y}` |
| Year only | `\citeyear{key}` |
| Author only | `\citeauthor{key}` |

**Further Reading** sections stay prose, with names replaced by `\citet{...}`.

**Epigraph attributions** keep their typeset form (author, *title*, year). Add `\nocite{key}` after the epigraph so the work appears in the bibliography.

## Keys
Format: **first author's surname + year + first significant word of the title**, all lowercase ASCII.
- Drop accents and particles: de Paula → `depaula`; Hernández → `hernandez`.
- Skip a/an/the/on/of/la/el/o/a/de in the title.
- Institutional authors use a short name: `cepal1990transformacion`, `worldbank2022poverty`.

Examples:
- `acemoglu2012network`
- `golub2010naive`
- `prebisch1949desarrollo`
- `fernandez2009women`
- `cardoso1969dependencia`

**Translations and editions.** Cite the edition you quote. If you mention both the original and a translation, use two entries, e.g. `cardoso1969dependencia` and `cardoso1979dependency`.

## Entries
- **Types:** `@article`, `@book`, `@incollection`, `@techreport` (working papers), `@misc` (datasets, websites).
- **Required:**
  - author (`and`-separated, "Surname, Given");
  - title;
  - year;
  - journal / publisher / booktitle as appropriate.
- **Include when known:** volume, number, pages, doi.
- **Encoding:** non-ASCII as LaTeX escapes, e.g. `{\'e}`, `{\~n}`, `{\c{c}}`. Protect capitals in titles with braces, e.g. `{L}atin {A}merica`.
- **Never invent a field.** If you are not sure of a volume or page, omit it. If you are not sure the work exists as described, or of its venue or year, add `x-verify = {what to check}`; BibTeX ignores unknown fields. Never create an entry for a work you cannot identify.
