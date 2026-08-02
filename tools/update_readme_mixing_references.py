"""Update README mixing plots with citations and references.

This script turns the generated LaTeX captions and BibTeX entries for the
mixing plots into a Markdown block that GitHub can render.
"""

from __future__ import annotations

import re
from pathlib import Path

from pylatexenc.latex2text import LatexNodes2Text


ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
TEX_FILE = ROOT / "tex_files" / "mixing_plots.tex"
BIB_FILE = ROOT / "tex_files" / "mixing_plots.bib"

BEGIN = "<!-- BEGIN GENERATED MIXING REFERENCES -->"
END = "<!-- END GENERATED MIXING REFERENCES -->"


def read_balanced(text: str, start: int, opener: str, closer: str, *, respect_quotes: bool = False) -> tuple[str, int]:
    depth = 0
    in_quote = False
    value = []
    i = start
    while i < len(text):
        char = text[i]
        if char == "\\":
            value.append(text[i : i + 2])
            i += 2
            continue
        if respect_quotes and char == '"':
            in_quote = not in_quote
            value.append(char)
            i += 1
            continue
        if respect_quotes and in_quote:
            value.append(char)
            i += 1
            continue
        if char == opener:
            depth += 1
            if depth > 1:
                value.append(char)
        elif char == closer:
            depth -= 1
            if depth == 0:
                return "".join(value), i + 1
            value.append(char)
        else:
            value.append(char)
        i += 1
    raise ValueError("Unbalanced BibTeX field")


def read_quoted(text: str, start: int) -> tuple[str, int]:
    value = []
    i = start + 1
    while i < len(text):
        char = text[i]
        if char == "\\":
            value.append(text[i : i + 2])
            i += 2
            continue
        if char == '"':
            return "".join(value), i + 1
        value.append(char)
        i += 1
    raise ValueError("Unbalanced quoted BibTeX field")


def split_bib_entries(text: str) -> dict[str, dict[str, str]]:
    entries: dict[str, dict[str, str]] = {}
    pos = 0
    while True:
        match = re.search(r"@(\w+)\s*\{\s*([^,\s]+)\s*,", text[pos:])
        if not match:
            return entries

        entry_start = pos + match.start()
        entry_open = text.find("{", entry_start, pos + match.end())
        key = match.group(2)
        body, entry_end = read_balanced(text, entry_open, "{", "}", respect_quotes=True)
        entries[key] = parse_bib_fields(body[len(key) + 1 :])
        pos = entry_end


def parse_bib_fields(body: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    i = 0
    while i < len(body):
        while i < len(body) and body[i] in " \t\r\n,":
            i += 1
        name_start = i
        while i < len(body) and (body[i].isalnum() or body[i] in "_-"):
            i += 1
        if i == name_start:
            break
        name = body[name_start:i].lower()
        while i < len(body) and body[i].isspace():
            i += 1
        if i >= len(body) or body[i] != "=":
            break
        i += 1
        while i < len(body) and body[i].isspace():
            i += 1
        if body[i] == "{":
            value, i = read_balanced(body, i, "{", "}")
        elif body[i] == '"':
            value, i = read_quoted(body, i)
        else:
            value_start = i
            while i < len(body) and body[i] != ",":
                i += 1
            value = body[value_start:i].strip()
        fields[name] = value.strip()
    return fields


def tex_to_text(value: str) -> str:
    value = re.sub(r"\\ensuremath\s*\{([^{}]+)\}", r"\1", value)
    value = value.replace("---", "-").replace("--", "-")
    text = LatexNodes2Text(math_mode="text").latex_to_text(value)
    return " ".join(text.split())


def tex_to_markdown_text(value: str) -> str:
    parts = re.split(r"(\$[^$]*\$)", value)
    converted = []
    for part in parts:
        if part.startswith("$") and part.endswith("$"):
            converted.append(part)
        else:
            converted.append(LatexNodes2Text(math_mode="verbatim").latex_to_text(part))
    return "".join(converted)


def format_author(raw_author: str) -> str:
    authors = [part.strip() for part in raw_author.split(" and ")]
    if not authors:
        return ""
    if len(authors) > 1 and authors[1].lower() == "others":
        return f"{format_one_author(authors[0])} et al."
    formatted = [format_one_author(author) for author in authors]
    if len(formatted) == 1:
        return formatted[0]
    if len(formatted) == 2:
        return f"{formatted[0]} and {formatted[1]}"
    return f"{', '.join(formatted[:-1])}, and {formatted[-1]}"


def format_one_author(author: str) -> str:
    author = tex_to_text(author)
    if "," in author:
        last, first = [part.strip() for part in author.split(",", 1)]
        return f"{first} {last}".strip()
    return author


def format_reference(number: int, key: str, entry: dict[str, str]) -> str:
    author = format_author(entry.get("author", entry.get("collaboration", key)))
    title = tex_to_text(entry.get("title", ""))
    journal = tex_to_text(entry.get("journal", ""))
    volume = tex_to_text(entry.get("volume", ""))
    pages = tex_to_text(entry.get("pages", ""))
    year = tex_to_text(entry.get("year", ""))
    note = tex_to_text(entry.get("note", ""))
    doi = tex_to_text(entry.get("doi", ""))
    eprint = tex_to_text(entry.get("eprint", ""))

    parts = [author]
    if title:
        parts.append(f'"{title}"')

    venue = " ".join(part for part in [journal, volume] if part)
    if pages:
        venue = f"{venue}, {pages}" if venue else pages
    if year:
        venue = f"{venue} ({year})" if venue else f"({year})"
    if venue:
        parts.append(venue)
    if note:
        parts.append(note)

    links = []
    if doi:
        links.append(f"[doi:{doi}](https://doi.org/{doi})")
    if eprint:
        arxiv_id = eprint.replace("arXiv:", "")
        links.append(f"[arXiv:{arxiv_id}](https://arxiv.org/abs/{arxiv_id})")
    if links:
        parts.append("; ".join(links))

    return f'{number}. <a id="mixing-ref-{number}"></a>{"; ".join(part for part in parts if part)}.'


def extract_figures(text: str) -> list[dict[str, str]]:
    figures = []
    pattern = re.compile(
        r"\\includegraphics(?:\[[^\]]+\])?\{(?P<image>[^}]+)\}%?\s*"
        r"\\caption\{(?P<caption>.*?)\}%",
        re.S,
    )
    for match in pattern.finditer(text):
        image = Path(match.group("image")).with_suffix(".png")
        figures.append(
            {
                "image": str(Path("plots") / "mixing" / image.name),
                "caption": " ".join(match.group("caption").split()),
            }
        )
    return figures


def replace_citations(caption: str, numbers: dict[str, int]) -> str:
    def citation(match: re.Match[str]) -> str:
        keys = [key.strip() for key in match.group(1).split(",")]
        return " [" + ", ".join(f"[{numbers[key]}](#mixing-ref-{numbers[key]})" for key in keys) + "]"

    caption = re.sub(r"~?\\cite\{([^}]+)\}", citation, caption)
    caption = caption.replace(r"\ ", " ")
    caption = tex_to_markdown_text(caption)
    caption = re.sub(r"\s+", " ", caption)
    return caption.strip()


def collect_citation_numbers(figures: list[dict[str, str]]) -> dict[str, int]:
    numbers: dict[str, int] = {}
    for figure in figures:
        for group in re.findall(r"\\cite\{([^}]+)\}", figure["caption"]):
            for key in [key.strip() for key in group.split(",")]:
                if key not in numbers:
                    numbers[key] = len(numbers) + 1
    return numbers


def render_block() -> str:
    figures = extract_figures(TEX_FILE.read_text(encoding="utf-8"))
    numbers = collect_citation_numbers(figures)
    bib_entries = split_bib_entries(BIB_FILE.read_text(encoding="utf-8"))

    lines = [
        BEGIN,
        "",
        "### Mixing Plot Citations",
        "",
        "The captions and references below are generated from `tex_files/mixing_plots.tex` and `tex_files/mixing_plots.bib`.",
        "",
    ]

    for figure in figures:
        alt = Path(figure["image"]).stem
        caption = replace_citations(figure["caption"], numbers)
        lines.extend(
            [
                f"![{alt}]({figure['image']})",
                "",
                caption,
                "",
            ]
        )

    lines.extend(["#### References", ""])
    ordered_keys = sorted(numbers, key=numbers.get)
    for key in ordered_keys:
        if key not in bib_entries:
            raise KeyError(f"No BibTeX entry found for {key}")
        lines.append(format_reference(numbers[key], key, bib_entries[key]))
    lines.extend(["", END])
    return "\n".join(lines)


def update_readme() -> None:
    readme = README.read_text(encoding="utf-8")
    block = render_block()
    if BEGIN in readme and END in readme:
        readme = re.sub(
            rf"{re.escape(BEGIN)}.*?{re.escape(END)}",
            lambda _: block,
            readme,
            flags=re.S,
        )
    else:
        marker = re.search(r"\n---\s*\n## Limits on the dimension-six", readme)
        if not marker:
            raise ValueError("Could not find README insertion point")
        readme = readme[: marker.start()] + f"\n{block}" + readme[marker.start() :]
    README.write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    update_readme()
