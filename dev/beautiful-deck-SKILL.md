---
name: beautiful-deck
description: Create, restructure or audit a Beamer slide deck written in R Markdown (.Rmd) with APA citations from the project's .bib, held to the Rhetoric of Decks standard and a zero-warning compile.
---

# Beautiful Deck for R Markdown Beamer

A deck is not a document dump. It is rhetoric under time pressure. The standard is not coverage; it is whether a listener who hears each slide once leaves believing the claim.

This skill adapts Scott Cunningham's `beautiful_deck` skill and the essays it sits on (*The Rhetoric of Decks*, *Aristotle for Academic Decks*, in the MixtapeTools repository, github.com/scunning1975/MixtapeTools, folder `presentations/`) to William Garcia's workflow: decks are written in `.Rmd`, rendered with `rmarkdown::beamer_presentation` and xelatex on a fixed house template, read saved pipeline results, and cite a reference-manager `.bib` in APA style. If the repository is on disk or the network allows, read `presentations/rhetoric_of_decks.md` and `presentations/aristotle.md` once before the first deck in a project. The rules below are the operating standard either way.

## The standard

1. **Titles are assertions, not labels.** Each title is a sentence that states what the slide proves. Read in sequence, the titles alone carry the argument.
2. **One idea per slide.** One figure, one table, one equation, one number, or one contrast. If the slide needs "also", it needs a second slide.
3. **Ethos, pathos, logos, in that order, every time.** Earn trust, then make it matter, then prove it. The order fixes the opening of the deck and the order of the audit. Only the proportions change with the audience.
4. **Beauty is function.** Every element earns its place or is deleted. White space reads as confidence. Never shrink type to fit; cut or split.
5. **Compile clean. Zero warnings. Always.** No overfull or underfull boxes, no LaTeX, biber, pandoc or R warnings, no missing glyphs. The gate script decides, not an impression of the log.

Two rules come from this workflow and rank with the five:

6. **Nothing is typed from memory.** A citation is a key that exists in the project's `.bib`, and the author and year on the slide are printed by LaTeX from that entry. A statistic is inline R that reads a saved result, or a literature value checked in the source and recorded. What cannot be verified is a visible placeholder, never a guess.
7. **Existing files are not modified.** Pipeline scripts, the `.bib`, `beamer_preamble.tex` and existing decks stay as they are. A new or restructured deck is a new file. Edit an existing deck in place only when asked to.

## Step 0: Look before writing

Read the project before producing anything.

- Find the deck folder, any existing `.Rmd` decks, `beamer_preamble.tex`, the `.bib` named in their YAML, `CLAUDE.md`, and the saved result objects the decks load (in the dissertation: `.rds` files reached through `need()`, `out_path()` and `cfg`).
- Read an existing deck's setup chunk for loaders and palettes, and reuse them. In the dissertation, sex is always Female `#CB181D` and Male `#2171B5`, every other data colour comes from the fixed PuOr vector, and UW purple `#4B2E83` is reserved for titles and typography.
- Read the source content the deck will argue from (report, notebook, paper PDFs in the project). Every claim on a slide must trace to a file actually read in this session.
- Check the tools: `Rscript`, `rmarkdown`, `xelatex`, `biber`. With TinyTeX, missing pieces install with `tinytex::tlmgr_install(c("biber", "biblatex", "biblatex-apa", "appendixnumberbeamer"))`. If R or LaTeX is not available where this is running, say so and deliver the `.Rmd` unrendered and labelled as unverified. Never report a clean compile that was not run.

Then install the toolchain next to the deck, once per deck folder, from the blocks at the end of this skill:

| File | Rule |
| --- | --- |
| `beamer_preamble.tex` | The house look. If the folder already has one, use it untouched. Write the embedded copy only when none exists. |
| `beamer_addons.tex` | Write it if absent. If present, keep it: it may hold repairs confirmed for this project's `.bib`. |
| `deck_gate.R` | Write it if absent, exactly as given. |

Do not design a new theme. The same template across every deck of a project is part of the speaker's credibility. Design effort goes into structure and layout inside it.

## Step 1: Triage

Settle four things before outlining. Ask only for what the request and the project do not answer. When working unattended, choose, and state the choice at the top of the outline.

1. **Source.** What the deck argues from, and whether it is a new deck or a restructuring of an existing one.
2. **Audience.** Pick one row and commit to it.

| Audience | Ethos | Pathos | Logos | What that means on the slides |
| --- | --- | --- | --- | --- |
| Committee or co-author working deck | 30% | 10% | 60% | Decisions and open questions documented, uncertainty kept visible, a little more text allowed |
| Research seminar or defense | 40% | 10% | 50% | Design before results, one estimate at a time, a skeptic slide before the close |
| Conference talk, 15 to 20 minutes | 35% | 15% | 50% | Main claim on the title slide, no literature review, about one content slide per minute |
| Teaching lecture | 20% | 35% | 45% | Clarity over compression, worked steps, recap slides allowed |
| Policy or public audience | 25% | 45% | 30% | Human scale first, minimal jargon, one number per slide |

3. **Minutes.** Plan about one content slide per minute and pass `--minutes=N` to the gate. Everything else is backup.
4. **The one sentence.** What the audience should remember tomorrow. It goes on the closing slide and every main slide must serve it. If it cannot be written yet, the deck is not ready to outline.

## Step 2: Outline, then stop for approval

Write `<deck>_outline.md` beside the deck: the triage answers, then every slide as its title, its single element (figure, table, number, equation, contrast), and the saved object or source that element comes from. Show it and wait for approval before writing slides, unless working unattended. Fixing the arc here costs minutes; fixing it after forty slides costs the afternoon.

**Act I, three slides at most after the title, in this order.**

- *Title slide.* The talk title is itself the main claim, so the conclusion is never buried.
- *Ethos: why trust this.* The question, the evidence base, the method in plain words, and the honest boundary of the work. In the dissertation the boundary is that post-2019 quantities are target paths and counterfactuals, not forecasts.
- *Pathos: why care.* One number or one contrast at a scale this audience feels. For a seminar the feeling is intellectual: a fact that should not be true, or a gap between what a target requires and what has ever been observed.

Never open with an agenda, "Today I will talk about", a definition, or a literature review.

**Act II, the argument (logos), most of the deck.** Each slide makes one sub-claim and shows its evidence. Put the design before the results. When a method or a new quantity is introduced, move from intuition to formalism and never the reverse: plain-language story, then a concrete example, then a picture, then the equation. Code appears only in teaching and working decks, between picture and equation. Alternate dense slides with light ones so the deck breathes.

**Act III, resolution.** For seminar, defense and external audiences: a skeptic slide that names the strongest objection, answers it, and says what remains uncertain. This replaces a bulleted "Limitations" slide. Then implications if the audience needs them, then the closing slide with the one sentence. Never "Questions?", never "Thank you", never a bulleted summary.

**Backup.** Everything that is coverage rather than argument: full state tables, extra equations, robustness, sensitivity. Backup slides come after the closing slide and are numbered separately. A deck being restructured usually loses half its main slides to backup; that is the point.

Sections (`#` headings) make divider slides. Use three or four for a long deck and none for a short one.

## Step 3: Sources

### Citations

1. Cite only keys that exist in the `.bib` named in the YAML. Find a key by searching the `.bib` for the author or title words. Never construct a key from a pattern.
2. Never type an author name or a year in slide text. Write `\textcite{key}` for "Bennett et al. (2020)" and `\parencite{key}` for "(Bennett et al., 2020)". `\citet`, `\citep`, `@key` and `[@key]` also work. For "the Karlsson framework", write "the priority-conditions framework of `\textcite{Karlsson2025}`".
3. Never add or change an entry in the `.bib`. If a claim needs a source that is not there, put `[CITATION NEEDED: the claim, and what source is believed to exist]` on the slide and list it in the final report. The user adds it through their reference manager.
4. A citation supports a specific claim. Before attributing a number, finding or method to a source, open that source (the PDF in the project, or the `abstract` field of the entry) and confirm it, then record where on the slide: `<!-- src: Karlsson2025, abstract -->` or `<!-- src: Jamison2024, p. N -->`. If the source cannot be opened, cite it only for what its stored title and abstract establish, or mark `[VERIFY: ...]`.
5. Read the citation ledger the gate prints and confirm each row is the work intended. Suspected errors in the `.bib` (a family and given name swapped, a group parsed as a person, a missing year, a title stored in Title Case) are reported to the user with the evidence. After the user confirms, add one `\step` to the deck-specific block at the bottom of `beamer_addons.tex`. Never correct bibliographic data silently.
6. The reference list is generated: `citation_package: biblatex` appends "References" frames after the last slide, containing exactly the cited entries in APA 7 form. Do not write a references slide by hand and never use `\nocite{*}`.

What `beamer_addons.tex` does at render time, without touching the `.bib`: drops abstracts, keywords, months and ISSNs; escapes stray `%`, `&`, `#`; trims DOIs that carry an asset path or tracking suffix; keeps the first of several URLs; collapses page ranges like `e123-e123`; moves `issue` to `number`. Titles print exactly as stored, because automatic sentence-casing destroys acronyms and proper nouns. The gate lists each repair so the user can fix the library at the source.

### Numbers

- Every statistic from the analysis is inline R on a saved object, in titles as well as bodies: `` `r fmt(x)` ``. Never retype a number from a table, an earlier deck, the outline or memory.
- The deck loads saved results and never runs the pipeline. A missing file stops the render with the name of the script that produces it. No placeholder numbers, ever.
- Every title that states a result gets a guard in the setup chunk, `claim(condition, "title")`, so that if the data change, the render stops instead of the title becoming false.
- A value from the literature is typed once, with its `<!-- src: ... -->` comment on the slide.

## Step 4: Write the slides

Start from the skeleton at the end of this skill. For each slide: confirm its role in the outline, write the title, choose the one element, add only the text that labels or sets it up.

**Titles.** A full sentence with a verb, at most two lines (about 78 characters). "Results", "Methods", "Data sources", "Key findings", "Summary" and "Limitations" are labels and fail the gate. Write the finding ("Child mortality fell twice as fast as adult mortality"), the design ("Each state is compared with the North Atlantic, cause by cause"), or the boundary ("These are target paths, not forecasts").

**Text.** Text labels, sets up, or concludes; it does not explain, because the speaker does. Aim for under 25 words a slide; the gate flags more than 40. Bullets are for genuinely parallel items, four at most. Usually a structure is hiding in a list: a sequence, a contrast, a hierarchy, a causal chain. Show it with layout (two columns, a number, a table) instead.

**Callout boxes.** At most one per slide, used for its meaning: `block` grey for a definition, `alertblock` red for an open question, `exampleblock` green for a takeaway, `uncertaintyblock` amber for a caveat. A takeaway box that repeats the title is noise; the title is the takeaway.

**Figures.**

- One message per figure. The slide title states it, so the plot has no title of its own.
- Label directly instead of using a legend wherever possible. Colour the series the slide is about; grey the context. Reuse the project's palettes; the same variable keeps the same colour in every deck.
- Draw at slide size: `fig.width = 6`, `fig.height = 2.7` with a source line (2.9 without, 2.4 under a two-line title), `base_size = 11`, so figure text matches slide text. Nothing in a main-slide figure may be smaller than size 7. A plot with 51 state labels fails this: label the handful that carry the point, grey the rest, and put the full version in backup.
- Order categories by meaning (set factor levels), never by the alphabet.
- `dev = "cairo_pdf"`, so en dashes and accents in labels survive.
- Put `\sourceline{...}` under every figure and table: unit, population, years, data source.

**Tables.** `kbl(format = "latex", booktabs = TRUE) |> kable_styling(position = "center")`. At most about six rows and four columns on a main slide, with the cell or column that matters emphasised (`column_spec(..., bold = TRUE)`). No captions (the slide title is the caption), no `hold_position`, no `scale_down`, no reduced `font_size`. A table that does not fit is a backup slide or a figure.

**Equations.** One displayed equation per slide, after the picture that motivates it. Name every symbol in a line beneath. Mark the term that carries the insight with `\textcolor{uwpurple}{...}`.

**Never on a main slide:** `\tiny`, `\scriptsize`, `{.shrink}`, reduced table fonts, two figures, a figure beside a bullet list, decoration.

**R Markdown specifics.**

- Markdown is not parsed inside raw LaTeX environments such as `\begin{block}`. Use `\textbf{}` there. Citation commands and inline R work everywhere.
- Two columns: pandoc's `:::: columns` and `::: {.column width="48%"}` blocks.
- Closing slide: `## {.standout}` followed by the sentence. It is full-bleed, unnumbered, and ends the main numbering.
- A big number: `\bignumber{value}{what it counts}`.
- Keep `warning = TRUE` in the chunk options. A warning printed on a slide is a defect to fix at its cause, not to hide.
- A raw LaTeX command placed between two slides lands inside the earlier frame. Slides written with `results = 'asis'` are invisible to the gate's source checks; do not generate titles that way.
- A symbol the font lacks is silently dropped from the PDF; the gate reports it as a missing character. Use math mode for it.
- Diagrams in TikZ: add `\usepackage{tikz}` through `header-includes`. Give every node explicit `minimum width` and `minimum height`, put a direction (`above`, `below`, `left`, `right`) on every edge label, never use `scale` on a finished diagram, and define styles that take `#1` in the header, never inside a frame.

## Step 5: The gate

```
Rscript deck_gate.R path/to/deck.Rmd --check          # source checks, seconds
Rscript deck_gate.R path/to/deck.Rmd --minutes=20     # full gate
```

The full gate checks the source, deletes stale intermediates, renders, and reads the logs. It passes only at zero in every row: LaTeX errors, overfull and underfull boxes, LaTeX warnings, missing characters, biber warnings, R warnings printed in slides, R and pandoc warnings. It also fails on a cited key that is not in the `.bib`, a hand-typed attribution, a label title, a missing closing slide, hidden R warnings, and open `[VERIFY]` or `[CITATION NEEDED]` markers (`--draft` lets those through for a work-in-progress render and nothing else).

Run `--check` while writing and the full gate after every round of edits, not only at the end. A deck is not described as finished, clean or ready unless the last command run on it was a passing full gate.

`REVIEW` items do not fail the gate, but each one is either fixed or answered in the final report with a reason. `BIB NOTES` go to the user as a list of things to fix in the reference manager.

| The gate reports | Fix |
| --- | --- |
| Overfull `\vbox` | The slide holds too much. Cut text, lower `fig.height` by 0.2, or split. Not a smaller font. |
| Overfull `\hbox` | A token or table is too wide. Rephrase, shorten column heads, drop a column. |
| Underfull `\hbox` | A forced line break or a narrow column. Rephrase; remove the `\\`. |
| Missing character | The font lacks the glyph. Use math mode or another character. |
| Citation undefined, biber WARN | The key is not in the `.bib`. Citation rule 3. |
| R warning in slides | Fix the cause in the chunk: handle `NA`s explicitly, set factor levels. |
| Font cannot be found | Calibri is not installed. Do not edit the template. Tell the user; for a test render only, use a copy of the preamble with Carlito, which has the same metrics. |

**Circuit breaker.** If three different fixes for the same error have failed, stop editing. Quote the log line, list the three attempts and why each failed, and ask. Every further blind fix makes the file harder to repair.

## Step 6: Audits, in order

A passing gate means the deck compiles, not that it persuades. Rasterise the PDF (`pdftoppm -r 90 -png deck.pdf pages`, or `pdftools::pdf_convert`) and look at every page; do not audit from memory of the source. If the session can start a fresh reviewer that has not seen the drafting, give it the PDF, the `.Rmd` and these lists, and have it report without editing. Otherwise do the passes yourself as separate read-throughs.

**Ethos: would this audience trust the speaker?**

- Is the evidence base stated before the first result? Is the boundary of the work stated plainly?
- Does any title claim more than its slide shows? Does any conclusion outrun the design?
- Is the strongest objection named and answered, with what remains uncertain?
- Does every attributed claim carry a citation, and every statistic a source?

**Pathos: does this audience have a reason to care, and can they stay with it?**

- Does the opening raise a problem this audience recognises, at a scale they feel?
- Is there one clear point of tension the deck resolves?
- Where would a listener get lost or anxious, and is there a slide that lets them back in?

**Logos: is the argument visible?**

- Read the storyboard the gate prints. Do the titles alone make the argument, in order, with no gaps?
- Does each slide play exactly one inferential role, and each figure or table prove exactly one claim?
- Is the closing sentence earned by the slides before it?

**Attention.** Rate each main slide for what it gives the listener against what it costs them. Split the overloaded ones, merge or cut the empty ones. One deliberately dense slide for effect is allowed; more is drift.

**Graphics.** For every figure and table, from the rendered page:

- Do the numbers shown match the saved objects? Recompute two per figure in R.
- Are labels unclipped, unoverlapped and legible? Is anything below size 7?
- Are axes and category order sensible? Are colours consistent with the project palette across slides?
- Is the key value emphasised? Does every figure and table have its source line?
- On the References pages: does every entry show authors, year, title and source, with nothing printing as "n.d." or as a person's name that is really a group?

Apply every fix, then run the full gate again.

## Step 7: Deliver

Files: `<deck>.Rmd`, `<deck>.pdf`, `<deck>_outline.md`, plus the three toolchain files if they were new. Report briefly:

- the PDF, the count of content slides and backup slides, and the gate's last line;
- each `REVIEW` item and how it was resolved;
- `BIB NOTES` and suspected `.bib` errors, for the user to fix at the source;
- every `[CITATION NEEDED]` and `[VERIFY]` that remains, and any claim cut because it could not be sourced.

## Converting an existing deck (natbib, label titles)

1. Copy it to a new file; leave the original alone.
2. YAML: `citation_package: biblatex`; add `beamer_addons.tex` after the preamble in `in_header`; add `biblio-title: References`; remove any `biblio-style`.
3. Delete the hand-made references slide. Existing `\citep` and `\citet` keep working.
4. Set `warning = TRUE` and `dev = "cairo_pdf"` in the chunk options.
5. Run `--check`. It lists every hand-typed attribution, label title, overloaded slide and typed statistic. That list is the restructuring plan: write the one sentence, rebuild the arc as in Step 2, move coverage to backup.

## Toolchain files

### `beamer_addons.tex`

```latex
% ============================================================
%  beamer_addons.tex -- load AFTER beamer_preamble.tex
%  APA 7 citations and reference list (biblatex-apa + biber),
%  closing slide, big-number and source-line helpers, and
%  render-time repairs for reference-manager exports.
%  The .bib file itself is never edited.
% ============================================================
\providecolor{uwpurple}{HTML}{4B2E83}

% -- APA 7 citations and reference list ---------------------------------------
\usepackage[style=apa, backend=biber, natbib=true]{biblatex}
\DeclareFieldFormat{apacase}{#1}            % print titles exactly as stored
\setlength{\bibhang}{0pt}
\setlength{\bibitemsep}{0.6\baselineskip}
\renewcommand*{\bibfont}{\footnotesize\raggedright}
\setcounter{biburlbigbreakpenalty}{5000}    % break before a DOI, not inside it
\setcounter{biburlbreakpenalty}{5000}
\setbeamertemplate{frametitle continuation}{}

\makeatletter
% reference frames are top-aligned on every page
\AtBeginBibliography{%
  \global\beamer@frametopskipautobreak=.45cm plus .5\paperheight\relax
  \global\beamer@framebottomskipautobreak=0pt plus 1fill\relax}
% keep the space between "(year)." and the title under biblatex-apa
\AtBeginDocument{%
  \patchcmd{\abx@macro@title}{\bbx@tempa\newblock\unspace}{\bbx@tempa\newblock}{}{}}

% -- Closing slide:  ## {.standout}  ------------------------------------------
%    Full-bleed and unnumbered. It ends the main deck: every frame after it
%    (backup slides, references) is numbered separately, so the footer total
%    counts main slides only.
\usepackage{appendixnumberbeamer}
\renewcommand{\appendixname}{Backup}
\AtBeginPart{}                              % no automatic part page
\newif\ifdeck@standout
\define@key{beamerframe}{standout}[true]{%
  \deck@standouttrue
  \begingroup
  \setkeys{beamerframe}{c,noframenumbering}%
  \setbeamercolor{background canvas}{bg=uwpurple}%
  \setbeamercolor{normal text}{fg=white}%
  \setbeamertemplate{footline}{}%
  \usebeamercolor[fg]{normal text}%
  \centering\Large}
% must be prepended: pandoc's ignorenonframetext mode skips appended code
\pretocmd{\beamer@reseteecodes}{%
  \ifdeck@standout\endgroup\deck@standoutfalse
    \ifnum\currentgrouplevel>0 \aftergroup\appendix\else\appendix\fi
  \fi}{}{}
\makeatother

% -- One number that carries the slide:  \bignumber{80\%}{of the gap ...} ------
\newcommand{\bignumber}[2]{%
  \begin{center}
    {\fontsize{44}{48}\selectfont\bfseries\textcolor{uwpurple}{#1}}\par\vspace{6pt}
    {\large #2}\par
  \end{center}}

% -- Source line under a figure, table or number ------------------------------
\newcommand{\sourceline}[1]{\par\vspace{4pt}{\footnotesize\textcolor{black!60}{#1}\par}}

% -- Render-time repairs for reference-manager exports ------------------------
% Regex rules: spaces are stripped, so write \s or \x{20}; write % & # as
% \x{25} \x{26} \x{23}; never put < in a pattern.
\DeclareSourcemap{
  \maps[datatype=bibtex]{
    \map{% fields APA never prints, which often carry unescaped characters
      \step[fieldset=abstract, null]
      \step[fieldset=keywords, null]
      \step[fieldset=month, null]
      \step[fieldset=issn, null]
      \step[fieldset=pmid, null]
    }
    \map{% a note that only repeats the DOI
      \step[fieldsource=note, match=\regexp{^\s*doi:}, final]
      \step[fieldset=note, null]
    }
    \map{% escape % & # in printed fields
      \step[fieldsource=title, match=\regexp{\\([\x{25}\x{26}\x{23}])}, replace=\regexp{$1}]
      \step[fieldsource=title, match=\regexp{([\x{25}\x{26}\x{23}])}, replace=\regexp{\\$1}]
      \step[fieldsource=journal, match=\regexp{\\([\x{25}\x{26}\x{23}])}, replace=\regexp{$1}]
      \step[fieldsource=journal, match=\regexp{([\x{25}\x{26}\x{23}])}, replace=\regexp{\\$1}]
    }
    \map{% DOI with an asset path or tracking suffix; several URLs in one field
      \step[fieldsource=doi, match=\regexp{^(10\.\d{4,9}/[^\s;]+?)(?:/ASSET/.*|;.*)$}, replace=\regexp{$1}]
      \step[fieldsource=url, match=\regexp{^(\S+)\s.*$}, replace=\regexp{$1}]
    }
    \map{% page range "e123-e123"; volume and issue appended to the journal name
      \step[fieldsource=pages, match=\regexp{^(\w+)\s*-+\s*\1$}, replace=\regexp{$1}]
      \step[fieldsource=journal, match=\regexp{\s+\d{4}\s+\d+:\d+$}, replace=\regexp{}]
    }
    \map[overwrite=false]{% reference managers write the issue number to "issue"
      \step[fieldsource=issue, fieldtarget=number]
    }
    % -- Deck-specific repairs: one step per confirmed error in the .bib -------
    \map{
      \step[fieldsource=author, match=\regexp{and\s+NCD\s+Countdown}, replace=\regexp{and\x{20}{NCD\x{20}Countdown}}]
      \step[fieldsource=author, match=\regexp{Horiuchi\s+Shiro}, replace=\regexp{Shiro\x{20}Horiuchi}]
    }
  }
}
```

### `deck_gate.R`

```r
#!/usr/bin/env Rscript
# deck_gate.R -- source checks, clean render and zero-warning gate for an
# R Markdown Beamer deck. Base R + rmarkdown only. Edits nothing.
#   Rscript deck_gate.R path/to/deck.Rmd            full gate (checks + render + log)
#   Rscript deck_gate.R path/to/deck.Rmd --check    source checks only, no render
#   Rscript deck_gate.R path/to/deck.Rmd --draft    allow open [VERIFY ...] markers
#   add --minutes=20 to check the slide count against the talk length
args <- commandArgs(trailingOnly = TRUE)
rmd  <- args[!startsWith(args, "--")][1]
if (is.na(rmd) || !file.exists(rmd)) stop("usage: Rscript deck_gate.R deck.Rmd [--check] [--draft]")
check_only <- "--check" %in% args
draft      <- "--draft" %in% args
minutes    <- as.numeric(sub("^--minutes=", "", grep("^--minutes=", args, value = TRUE)[1]))
rmd  <- normalizePath(rmd, winslash = "/")
dir  <- dirname(rmd)
stem <- tools::file_path_sans_ext(basename(rmd))
fails <- character(); warns <- character(); notes <- character()
FAIL <- function(...) fails <<- c(fails, paste0(...))
WARN <- function(...) warns <<- c(warns, paste0(...))
NOTE <- function(...) notes <<- c(notes, paste0(...))
cut_to <- function(x, n) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "~"), x)
squish <- function(x) gsub(r"(\s+)", " ", x)

# ---- read the deck -----------------------------------------------------------
src   <- readLines(rmd, warn = FALSE, encoding = "UTF-8")
ydash <- which(trimws(src) == "---")
ystop <- if (length(ydash) >= 2) ydash[2] else 0
yaml  <- src[seq_len(ystop)]
fence <- grepl(paste0(r"(^\s*)", strrep("`", 3)), src)
k     <- cumsum(fence)
in_chunk   <- (k %% 2 == 1) | fence
chunk_open <- fence & (k %% 2 == 1) & grepl(r"(\{\s*r\b)", src) & !grepl(r"(include\s*=\s*F)", src)
text <- ifelse(in_chunk, "", src)
text[seq_len(ystop)] <- ""
text <- gsub("<!--.*?-->", "", text, perl = TRUE)

yval <- function(key) {
  m <- grep(paste0(r"(^\s*)", key, r"(\s*:)"), yaml, value = TRUE)
  if (!length(m)) return(NA_character_)
  gsub(r"(^["']|["']$)", "", trimws(sub("^[^:]*:", "", m[1])))
}
if (!isTRUE(yval("citation_package") == "biblatex"))
  FAIL("YAML: output must set citation_package: biblatex")
if (!any(grepl(r"(beamer_addons\.tex)", yaml)))
  FAIL("YAML: includes: in_header must list beamer_addons.tex after the preamble")
if (!isTRUE(yval("latex_engine") == "xelatex")) WARN("YAML: latex_engine is not xelatex")
if (any(grepl(r"(^\s*(biblio-style|natbiboptions|csl)\s*:)", yaml)))
  FAIL("YAML: remove biblio-style / natbiboptions / csl; beamer_addons.tex sets APA")
if (any(grepl(r"(warning\s*=\s*(FALSE|F)\b|suppressWarnings\s*\()", src)))
  FAIL("R warnings are hidden (warning = FALSE or suppressWarnings); show and fix them")

# ---- bibliography ------------------------------------------------------------
read_bib <- function(path) {
  x  <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
  ch <- strsplit(x, "")[[1]]
  depth <- cumsum(ch == "{") - cumsum(ch == "}")
  st <- gregexpr(r"(@[A-Za-z]+\s*\{)", x)[[1]]
  out <- list()
  if (st[1] < 0) return(out)
  for (i in seq_along(st)) {
    s <- st[i]; o <- s + attr(st, "match.length")[i] - 1
    if (s > 1 && depth[s - 1] != 0) next
    type <- tolower(gsub("[@{[:space:]]", "", substr(x, s, o)))
    if (type %in% c("comment", "string", "preamble")) next
    e <- o + which(depth[(o + 1):length(depth)] == 0)[1]
    if (is.na(e)) next
    idx <- (o + 1):(e - 1)
    cm  <- idx[ch[idx] == "," & depth[idx] == 1]
    if (!length(cm)) next
    key <- trimws(paste(ch[(o + 1):(cm[1] - 1)], collapse = ""))
    b <- c(cm, e)
    f <- list()
    for (j in seq_len(length(b) - 1)) {
      if (b[j + 1] - b[j] < 2) next
      p <- paste(ch[(b[j] + 1):(b[j + 1] - 1)], collapse = "")
      m <- regmatches(p, regexec(r"(^\s*([A-Za-z_-]+)\s*=\s*(.*)$)", p))[[1]]
      if (length(m) == 3) f[[tolower(m[2])]] <- squish(gsub(r"(^[{"]|[}"]$)", "", trimws(m[3])))
    }
    out[[length(out) + 1]] <- list(key = key, type = type, f = f)
  }
  out
}
authors <- function(a) if (is.null(a)) character() else trimws(strsplit(a, r"(\s+and\s+)")[[1]])
surname <- function(n) gsub("[{}]", "", if (grepl(",", n)) trimws(sub(",.*", "", n)) else
  if (grepl(r"(^\{.*\}$)", n)) n else sub(r"(.*\s)", "", n))
addons <- if (file.exists(file.path(dir, "beamer_addons.tex")))
  readLines(file.path(dir, "beamer_addons.tex"), warn = FALSE) else character()
bibrel <- yval("bibliography")
bibpath <- if (is.na(bibrel)) NA else file.path(dir, bibrel)
bib <- list()
if (is.na(bibpath) || !file.exists(bibpath)) {
  FAIL("YAML: bibliography file not found: ", bibrel)
} else bib <- read_bib(bibpath)
bibkeys <- vapply(bib, function(e) e$key, "")
if (anyDuplicated(bibkeys)) FAIL("duplicate keys in .bib: ",
                                 paste(unique(bibkeys[duplicated(bibkeys)]), collapse = ", "))

# ---- citations ---------------------------------------------------------------
txt  <- paste(text, collapse = "\n")
cite_re <- r"(\\[A-Za-z]*cite[A-Za-z]*\*?(\[[^]]*\])*\{[^}]*\})"
lat  <- regmatches(txt, gregexpr(cite_re, txt))[[1]]
lkey <- trimws(unlist(strsplit(gsub(r"(^.*\{|\}$)", "", lat), ",")))
pan  <- regmatches(txt, gregexpr(r"((?<![\w\\.])@[A-Za-z][\w:+/-]*)", txt, perl = TRUE))[[1]]
cited <- unique(c(lkey, sub("^@", "", pan)))
cited <- cited[nzchar(cited)]
if ("*" %in% cited) FAIL("nocite{*} prints the whole library; cite only what a slide uses")
cited <- setdiff(cited, "*")
for (m in setdiff(cited, bibkeys)) FAIL("cited key not in ", basename(bibpath), ": ", m)

cat("\n== CITATION LEDGER (read each row against the source) ==\n")
ent <- bib[match(intersect(cited, bibkeys), bibkeys)]
for (e in ent) {
  f <- e$f; au <- authors(f$author)
  yr <- if (!is.null(f$year)) f$year else if (!is.null(f$date)) substr(f$date, 1, 4) else "n.d."
  cat(sprintf("  %-28s %s%s (%s). %s | %s\n", e$key,
              if (length(au)) surname(au[1]) else "[no author]",
              if (length(au) > 1) sprintf(" +%d", length(au) - 1) else "", yr,
              cut_to(if (is.null(f$title)) "[no title]" else f$title, 58),
              if (is.null(f$journal)) e$type else f$journal))
  if (is.null(f$author) && is.null(f$editor)) WARN(e$key, ": no author in .bib")
  if (is.null(f$year) && is.null(f$date)) WARN(e$key, ": no year in .bib (prints as n.d.)")
  if (is.null(f$title)) WARN(e$key, ": no title in .bib")
  if (e$type == "article" && is.null(f$journal)) WARN(e$key, ": article without journal")
  if (!is.null(f$doi) && !grepl(r"(^10\.\d{4,9}/[^\s;]+$)", f$doi, perl = TRUE) ||
      isTRUE(grepl("/ASSET/", f$doi)))
    NOTE(e$key, ": malformed DOI in .bib (trimmed at render)")
  if (isTRUE(grepl(r"(\s)", f$url))) NOTE(e$key, ": several URLs in one field (first kept at render)")
  if (any(grepl(r"((?<!\\)[%&#])", c(f$title, f$journal), perl = TRUE)))
    NOTE(e$key, ": unescaped % & # in title/journal (escaped at render)")
  if (isTRUE(grepl(r"(^(\w+)\s*-+\s*\1$)", f$pages, perl = TRUE)))
    NOTE(e$key, ": page range repeats one value (collapsed at render)")
  if (isTRUE(grepl(r"(\s\d{4}\s+\d+:\d+$)", f$journal)))
    NOTE(e$key, ": volume/issue appended to journal name (trimmed at render)")
  if (!is.null(f$title)) {
    w <- regmatches(f$title, gregexpr("[A-Za-z][A-Za-z'-]{3,}", f$title))[[1]][-1]
    if (length(w) >= 3 && mean(grepl("^[A-Z]", w)) > 0.6)
      WARN(e$key, ": title stored in Title Case; APA prints sentence case (fix in the reference manager)")
  }
  grp <- au[!grepl(",", au) & !grepl(r"(^\{.*\}$)", au) &
            grepl(r"(^[A-Z]{2,}\b|Collaborat|Group|Consortium|Committee|Countdown|Network|Organization|Team$)", au)]
  for (g in grp) {
    if (any(grepl(gsub(r"(\s+)", r"(\\s+)", g), addons, fixed = TRUE)))
      NOTE(e$key, ": group author '", g, "' (braced at render)")
    else WARN(e$key, ": author '", g, "' looks like a group name; needs {braces} or a deck-specific repair")
  }
}
jn <- unlist(lapply(bib, function(e) e$f$journal))
jnorm <- tolower(sub(r"(^the\s+)", "", jn, ignore.case = TRUE))
for (j in unique(jnorm[duplicated(jnorm)])) {
  v <- unique(jn[jnorm == j])
  if (length(v) > 1) NOTE("journal spelled inconsistently in .bib: ", paste(v, collapse = " / "))
}

# hand-typed attributions: author and year must be printed by a cite command
plain <- gsub(cite_re, " ", txt)
plain <- gsub("`r [^`]*`", " ", plain)
nm <- r"([^\s\d,;:()&\\{}]*)"
att_re <- paste0(r"(\b([A-Z][a-z])", nm, r"()((?:\s+(?:and|&)\s+[A-Z])", nm,
                 r"()?)((?:\s+et\s+al\.?)?)\\?,?\s*\(?((?:19|20)\d{2})[a-z]?\))")
att <- unique(regmatches(plain, gregexpr(att_re, plain, perl = TRUE))[[1]])
first_sn <- vapply(bib, function(e) { a <- authors(e$f$author); if (length(a)) surname(a[1]) else "" }, "")
byear <- vapply(bib, function(e) if (is.null(e$f$year)) "" else e$f$year, "")
for (a in att) {
  sn <- sub(paste0("^([A-Z][a-z]", nm, ").*$"), r"(\1)", a, perl = TRUE)
  yr <- regmatches(a, regexpr(r"((19|20)\d{2})", a))
  hit <- which(first_sn == sn & byear == yr)
  if (length(hit)) FAIL("hand-typed attribution '", squish(a), r"(': write \textcite{)", bibkeys[hit[1]], "}")
  else if (grepl(r"(et\s+al)", a)) FAIL("attribution '", squish(a), "' matches no entry in the .bib")
  else WARN("possible uncited attribution: '", squish(a), "'")
}

# ---- slides: titles and density ---------------------------------------------
h1 <- grepl(r"(^#\s+\S)", text)
h2 <- grepl(r"(^##(\s|$))", text) & !grepl("^###", text)
attr_re <- r"(\s*\{[^}]*\}\s*$)"
starts <- which(h2)
slides <- lapply(seq_along(starts), function(i) {
  a <- starts[i]; b <- if (i < length(starts)) starts[i + 1] - 1 else length(text)
  nxt <- which(h1[a:b]); if (length(nxt)) b <- a + nxt[1] - 2
  head <- sub(r"(^##\s*)", "", text[a])
  attr <- regmatches(head, regexpr(attr_re, head))
  body <- if (b > a) text[(a + 1):b] else character()
  bt <- paste(body, collapse = "\n")
  w <- paste(body[!grepl(r"(^\s*:::)", body)], collapse = "\n")
  w <- gsub(r"(\$\$[\s\S]*?\$\$|\$[^$\n]*\$)", " ", w, perl = TRUE)
  w <- gsub(cite_re, " ", w); w <- gsub("`r [^`]*`", " N ", w)
  w <- gsub(r"(\\(begin|end)\{[^}]*\}|\\[A-Za-z]+\*?(\[[^]]*\])?)", " ", w)
  sec <- sub(r"(^#\s+)", "", sub(attr_re, "", text[max(c(0, which(h1 & seq_along(h1) < a)))][1]))
  list(line = a, title = trimws(sub(attr_re, "", head)), attr = if (length(attr)) attr else "",
       sec = if (is.na(sec)) "" else sec,
       words = length(regmatches(w, gregexpr(r"(\S*[A-Za-z]\S*)", w))[[1]]),
       bullets = sum(grepl(r"(^\s*([-*+]|\d+\.)\s|\\item\b)", body)),
       visuals = sum(chunk_open[a:b]) + sum(grepl(r"(\\includegraphics|!\[)", body)),
       eqs = lengths(regmatches(bt, gregexpr(r"(\$\$)", bt))) %/% 2 +
         sum(grepl(r"(\\begin\{(equation|align))", body)),
       boxes = sum(grepl(r"(\\begin\{(block|alertblock|exampleblock|uncertaintyblock)\})", body)),
       sourced = any(grepl(r"(<!--\s*src:)", src[a:b])),
       first = squish(paste(body[nzchar(trimws(body))], collapse = " ")),
       stats = regmatches(w, gregexpr(r"(\d+(\.\d+)?\s*\\?%|\b\d+\.\d+\b|\b\d{1,3}(,\d{3})+\b)", w))[[1]])
})
is_close <- vapply(slides, function(s) grepl("standout", s$attr), NA)
is_bk <- vapply(slides, function(s) grepl("^(appendix|backup)", s$sec, ignore.case = TRUE), NA)
n_main <- min(c(which(is_close), which(is_bk) - 1, length(slides)))
secs <- unique(vapply(slides[seq_len(n_main)], function(s) s$sec, ""))
if (length(secs) > 5) WARN(length(secs), " section dividers in the main deck; three or four acts are enough")
label_re <- paste0("^(introduction|background|motivation|overview|outline|agenda|roadmap|aims?|objectives?|",
  "research questions?|data|data (sources?|overview)|methods?|methodology|approach|results?|findings|",
  "key findings|main results|discussion|interpretation|limitations?|conclusions?|summary|next steps|",
  "future work|planned extensions|recommendations|implications|contributions?|what is new|",
  "literature( review)?|references|questions|thank you|thanks|appendix)$|",
  r"(^(results?|methods?|discussion|limitations?|appendix|summary)\s*[:-])")
if (!any(is_close)) FAIL("no closing slide: end the main deck with  ## {.standout}  and the one sentence")
cat("\n== STORYBOARD (the titles alone must carry the argument) ==\n")
for (i in seq_along(slides)) {
  s <- slides[[i]]; main <- i <= n_main
  if (nzchar(s$sec) && (i == 1 || s$sec != slides[[i - 1]]$sec)) cat("     --", s$sec, "--\n")
  tag <- if (is_close[i]) "CLOSE " else if (main) sprintf("%5d ", i) else "  bk  "
  cat(tag, if (nzchar(s$title)) s$title else cut_to(s$first, 110), "\n")
  if (!main || is_close[i]) next
  where <- sprintf("slide %d (line %d) '%s': ", i, s$line, cut_to(s$title, 40))
  bare <- tolower(trimws(gsub("[?!.]+$", "", s$title)))
  nw <- length(strsplit(bare, r"(\s+)")[[1]])
  if (!nzchar(bare)) FAIL(where, "no title")
  else if (grepl(label_re, bare) || nw <= 2) FAIL(where, "title is a label; state the claim")
  if (nchar(s$title) > 78) WARN(where, "title over 78 characters; it will wrap")
  if (s$words > 40) WARN(where, s$words, " words on the slide; cut or split")
  if (s$bullets > 4) WARN(where, s$bullets, " bullets; find the structure or split")
  if (s$visuals > 1) WARN(where, s$visuals, " figures/tables; one per slide")
  if (s$visuals >= 1 && s$bullets >= 3) WARN(where, "figure or table plus a bullet list")
  if (s$eqs > 1) WARN(where, s$eqs, " displayed equations; one per step")
  if (s$boxes > 1) WARN(where, s$boxes, " callout boxes; one at most")
  if (length(s$stats) && !s$sourced)
    WARN(where, "typed statistic(s) ", paste(unique(trimws(s$stats)), collapse = ", "),
         "; use inline R, or verify in the source and record  <!-- src: key, page -->")
}
mw <- vapply(slides[seq_len(n_main)], function(s) s$words, 0)
n_content <- n_main - any(is_close)
cat(sprintf("\n  content slides: %d (plus title%s)   words per slide: median %g, max %g\n", n_content,
            if (any(is_close)) " and closing" else "", stats::median(mw), max(mw, 0)))
if (!is.na(minutes) && n_content > minutes)
  WARN(n_content, " content slides for a ", minutes, "-minute talk; plan about one per minute")
open <- grep(r"(\[(VERIFY|CITATION NEEDED|TODO))", src)
if (length(open)) (if (draft) WARN else FAIL)(
  length(open), " open marker(s) [VERIFY/CITATION NEEDED/TODO] at line(s) ", paste(open, collapse = ", "))

report <- function() {
  if (length(notes)) cat("\n== BIB NOTES (handled at render; fix at the source when convenient) ==\n",
                         paste0("  - ", notes, "\n"), sep = "")
  if (length(warns)) cat("\n== REVIEW (resolve or justify each) ==\n", paste0("  - ", warns, "\n"), sep = "")
  if (length(fails)) cat("\n== FAIL ==\n", paste0("  x ", fails, "\n"), sep = "")
  cat(sprintf("\nGATE %s: %d failure(s), %d review item(s)\n",
              if (length(fails)) "FAILED" else "PASSED", length(fails), length(warns)))
  quit(status = if (length(fails)) 1 else 0, save = "no")
}
if (check_only || length(fails)) report()

# ---- clean render ------------------------------------------------------------
if (!nzchar(Sys.which("biber"))) {
  FAIL("biber not found; with TinyTeX run tinytex::tlmgr_install(c('biber','biblatex','biblatex-apa'))")
  report()
}
p   <- function(ext) file.path(dir, paste0(stem, ".", ext))
aux <- c("aux", "bbl", "bcf", "blg", "log", "nav", "snm", "toc", "out", "run.xml", "vrb", "tex")
unlink(p(aux))
options(tinytex.clean = FALSE)
rw <- character()
withCallingHandlers(
  tryCatch(rmarkdown::render(rmd, quiet = TRUE, envir = new.env(),
                             output_options = list(keep_tex = TRUE, pandoc_args = "--fail-if-warnings")),
           error = function(e) FAIL("render stopped: ", conditionMessage(e))),
  warning = function(w) { rw <<- c(rw, conditionMessage(w)); invokeRestart("muffleWarning") })
rd <- function(ext) if (file.exists(p(ext))) readLines(p(ext), warn = FALSE) else character()
g  <- function(re, x) grep(re, x, value = TRUE, useBytes = TRUE)
log <- rd("log"); blg <- rd("blg"); tex <- rd("tex"); bbl <- rd("bbl")
wre <- r"(^(LaTeX( Font)? Warning|(Package|Class|Module) \S+ Warning|pdfTeX warning))"
counts <- c(
  "LaTeX errors"         = length(g("^!", log)),
  "overfull/underfull"   = length(g(r"(^(Over|Under)full \\[hv]box)", log)),
  "LaTeX warnings"       = length(g(wre, log)),
  "missing characters"   = length(g("^Missing character", log)),
  "biber warnings"       = length(g("WARN|ERROR", blg)),
  "R warnings in slides" = length(g(r"(^(## |\\#\\# )(Warning|Error))", tex)),
  "R/pandoc warnings"    = length(rw))
cat("\n== COMPILE ==\n"); for (n in names(counts)) cat(sprintf("  %-22s %d\n", n, counts[[n]]))
if (!length(log)) FAIL("no LaTeX log produced")
for (n in names(counts)[counts > 0]) FAIL(counts[[n]], " ", n)
show <- c(g(paste0("^!|^(Over|Under)full|^Missing character|", wre), log), g("WARN|ERROR", blg), rw)
if (length(show)) cat(paste0("    ", utils::head(show, 25), "\n"), sep = "")
inbbl <- unique(sub(r"(^\s*\\entry\{([^}]*)\}.*$)", r"(\1)", g(r"(^\s*\\entry\{)", bbl)))
for (m in setdiff(cited, inbbl)) FAIL("cited but absent from the reference list: ", m)
for (m in setdiff(inbbl, cited)) WARN("in the reference list but not cited in slide text: ", m)
pg <- sub(r"(.*\((\d+) pages?.*)", r"(\1)", g("^Output written on", log))
cat(sprintf("  pages: %s   references: %d\n", if (length(pg)) pg[1] else "?", length(inbbl)))
if (!length(fails)) { unlink(p(aux)); unlink(file.path(dir, paste0(stem, "_files")), recursive = TRUE) }
report()
```

### `beamer_preamble.tex` (only when the deck folder has none)

```latex
% ============================================================
%  beamer_preamble.tex -- house look for R Markdown Beamer decks
%  Default Beamer theme, UW purple titles and bullets, four
%  colour-coded callout boxes, plain footer. Engine: xelatex.
% ============================================================
\usetheme{default}
\usecolortheme{default}
\usefonttheme{professionalfonts}

% -- Fonts (xelatex). Carlito is the metric twin if Calibri is not installed. -
\usepackage{fontspec}
\setmainfont{Calibri}[
  BoldFont       = Calibri Bold,
  ItalicFont     = Calibri Italic,
  BoldItalicFont = Calibri Bold Italic
]
\setsansfont{Calibri}
\setmonofont{Courier New}

% -- Packages -----------------------------------------------------------------
\usepackage{booktabs}
\usepackage{xcolor}
\usepackage{colortbl}
\usepackage{array}
\usepackage{graphicx}
\usepackage{ragged2e}
\usepackage{longtable}
\usepackage{microtype}

% -- Palette ------------------------------------------------------------------
\definecolor{uwpurple}{HTML}{4B2E83}
\definecolor{uwgold}{HTML}{85754D}
\definecolor{rtslgreen}{HTML}{74C476}
\definecolor{baugrey}{HTML}{999999}
\definecolor{uncertaintyyellow}{HTML}{FFF3CD}
\definecolor{uncertaintyborder}{HTML}{856404}

\setbeamercolor{structure}{fg=uwpurple}
\setbeamercolor{title}{fg=uwpurple}
\setbeamercolor{subtitle}{fg=uwpurple}
\setbeamercolor{author}{fg=black}
\setbeamercolor{date}{fg=black}
\setbeamercolor{institute}{fg=black}
\setbeamercolor{frametitle}{fg=uwpurple}
\setbeamercolor{normal text}{fg=black}

% -- Callout boxes: each colour has one meaning -------------------------------
% block = grey (definition) | alertblock = red (open question)
% exampleblock = green (takeaway) | uncertaintyblock = amber (caveat)
\setbeamercolor{block title}{bg=baugrey!15, fg=black}
\setbeamercolor{block body}{bg=baugrey!6}
\setbeamercolor{block title alerted}{bg=red!15, fg=red!70!black}
\setbeamercolor{block body alerted}{bg=red!6}
\setbeamercolor{block title example}{bg=rtslgreen!18, fg=rtslgreen!80!black}
\setbeamercolor{block body example}{bg=rtslgreen!6}
\newenvironment{uncertaintyblock}[1]{%
  \setbeamercolor{block title}{bg=uncertaintyyellow, fg=uncertaintyborder}%
  \setbeamercolor{block body}{bg=uncertaintyyellow!30}%
  \begin{block}{#1}%
}{\end{block}}

% -- Bullets ------------------------------------------------------------------
\setbeamertemplate{itemize item}{\textcolor{uwpurple}{\textbullet}}
\setbeamertemplate{itemize subitem}{\textcolor{uwpurple}{--}}
\setbeamertemplate{enumerate item}{\textcolor{uwpurple}{\insertenumlabel.}}

% -- Footer: author | project | frame number ----------------------------------
\newcommand{\deckfootcenter}{PhD Dissertation}   % edit per project
\setbeamertemplate{footline}{%
  \leavevmode
  \hbox{%
    \begin{beamercolorbox}[wd=0.45\paperwidth, ht=2.25ex, dp=1ex, left]{author in head/foot}%
      \hspace{1em}\textcolor{black!60}{\footnotesize\insertshortauthor}%
    \end{beamercolorbox}%
    \begin{beamercolorbox}[wd=0.45\paperwidth, ht=2.25ex, dp=1ex, center]{title in head/foot}%
      \textcolor{black!60}{\footnotesize\deckfootcenter}%
    \end{beamercolorbox}%
    \begin{beamercolorbox}[wd=0.10\paperwidth, ht=2.25ex, dp=1ex, right]{date in head/foot}%
      \textcolor{black!60}{\footnotesize\insertframenumber{}/\inserttotalframenumber}%
      \hspace{0.5em}%
    \end{beamercolorbox}%
  }%
  \vskip0pt%
}
\setbeamertemplate{navigation symbols}{}

% -- Section divider (one per # heading) --------------------------------------
\AtBeginSection[]{%
  \begin{frame}
    \vfill
    \centering
    \begin{beamercolorbox}[sep=8pt, center, shadow=false, rounded=true]{block title}
      \usebeamerfont{title}\insertsectionhead\par
    \end{beamercolorbox}
    \vfill
  \end{frame}
}
```

### Deck skeleton

````markdown
---
title: "<The main claim of the talk, as one sentence>"
subtitle: "<Occasion, e.g. Dissertation Aim 1, committee meeting>"
author: "William Garcia"
institute: "University of Washington"
date: "`r format(Sys.Date(), '%B %Y')`"
output:
  beamer_presentation:
    includes:
      in_header: [beamer_preamble.tex, beamer_addons.tex]
    slide_level: 2
    latex_engine: xelatex
    citation_package: biblatex
bibliography: ../library/references.bib
biblio-title: References
classoption: "aspectratio=169"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = FALSE, message = FALSE, warning = TRUE,   # warnings stay visible; the gate counts them
  dev = "cairo_pdf", fig.align = "center",
  fig.width = 6.0, fig.height = 2.7, out.width = "100%"
)
library(ggplot2)
library(kableExtra)

# Saved results only. Rendering never runs the analysis and never shows a
# placeholder number. Use the project's own loader if it has one.
need <- function(path) {
  if (!file.exists(path)) stop("Missing result file: ", path, " (run the pipeline first)", call. = FALSE)
  readRDS(path)
}
# res <- need(here::here("output", "processed", "<object>.rds"))

# Reuse the project's palettes when they exist; the same variable keeps the
# same colour in every figure of every deck.
deck_font <- "Calibri"
col_key   <- "#542788"
col_muted <- "grey70"
theme_deck <- theme_minimal(base_size = 11, base_family = deck_font) +
  theme(panel.grid.minor = element_blank(), axis.title = element_blank())

# A title that states a result is guarded: if the data stop supporting it,
# the render stops instead of the slide saying something false.
claim <- function(ok, title) {
  if (!isTRUE(ok)) stop("Data no longer support the title: ", title, call. = FALSE)
}
# claim(n_on_pace == 0, "No state is on pace to halve premature death by 2050")
```

## <Ethos: the question, the evidence base and its boundary, as a claim>

:::: columns
::: {.column width="48%"}
**<What was estimated>**

<Data, units, years: one line.>
:::
::: {.column width="48%"}
**<How>**

<Method in plain words, with \textcite{key} for the design it follows.>
:::
::::

\medskip

\begin{uncertaintyblock}{Boundary}
<What this analysis is not, in one sentence.>
\end{uncertaintyblock}

## <Pathos: why it matters to this audience, stated with the number>

\bignumber{`r "<inline R>"`}{<what the number counts, in plain words>}

\sourceline{<Source and year. \parencite{key}>}

## <Logos: the first finding, as a sentence with its number>

```{r fig-finding}
# One message. No ggplot title: the slide title is the figure title.
# Direct labels, logical category order, key series in col_key, context in col_muted.
```

\sourceline{<Unit, population, years, data source.>}

## <The comparison that matters, as a sentence>

```{r tab-finding}
# kbl(tab, format = "latex", booktabs = TRUE, col.names = c(...)) |>
#   kable_styling(position = "center") |>
#   column_spec(3, bold = TRUE, color = "#4B2E83")
```

\sourceline{<What the columns are; uncertainty shown how.>}

## <Ethos again: the strongest objection, named before the audience raises it>

:::: columns
::: {.column width="48%"}
**A skeptic would say**

<The objection, stated fairly.>
:::
::: {.column width="48%"}
**What the evidence supports**

<The response, and what remains uncertain.>
:::
::::

## {.standout}

<The one sentence the audience should remember tomorrow.>

# Backup

## <Backup slide title: a claim here too>

<Detail that answers a likely question: full tables, equations, robustness.>
````
