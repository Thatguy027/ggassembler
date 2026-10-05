"""``ggasm`` console entry point."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .core import baseline
from .core.library import Library, PlasmidEntry
from .core.seqio import write_genbank
from .levels import level2_cassette


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ggasm", description="GG Assembler")
    sub = parser.add_subparsers(dest="command")

    scan = sub.add_parser("scan", help="index folders of GenBank files and print what it found")
    scan.add_argument("folder", type=Path, nargs="+", help="one or more folders to index as one library")
    scan.add_argument("--no-recursive", dest="recursive", action="store_false",
                      help="do not descend into subfolders")
    scan.add_argument("--rescan", action="store_true", help="ignore the cached index")
    scan.add_argument("--json", action="store_true", help="print the index as JSON")
    scan.add_argument("--type", dest="part_type", help="only show parts of this type")
    scan.add_argument("--unknown", action="store_true", help="only show unclassified plasmids")

    build = sub.add_parser(
        "build", help="assemble a cassette from named library parts and write the .gbk"
    )
    build.add_argument("folder", type=Path, nargs="+", help="library folder(s)")
    build.add_argument("-p", "--part", action="append", default=[], required=True,
                       metavar="NAME", help="a part, repeated once per slot, in circle order")
    build.add_argument("-o", "--out", type=Path, help="where to write the product")
    build.add_argument("--name", default="cassette", help="name for the assembled plasmid")
    build.add_argument("--split-3", action="store_true", help="use 3a + 3b instead of 3")
    build.add_argument("--split-4", action="store_true", help="use 4a + 4b instead of 4")
    build.add_argument("--split-8", action="store_true", help="use 8a + 8b (integration)")
    build.add_argument("--composite-left", action="store_true", help="one 2\u00b73\u00b74 slot")
    build.add_argument("--composite-right", action="store_true", help="one 6\u00b77\u00b78 slot")

    check = sub.add_parser(
        "check",
        help="fail if detection finds new conflicts or unrecognised files",
    )
    check.add_argument("folder", type=Path, nargs="+", help="library folder(s) to index")
    check.add_argument("--baseline", type=Path, default=Path(".ggasm-baseline.json"),
                       help="the committed baseline to compare against")
    check.add_argument("--update", action="store_true",
                       help="write the current state as the new baseline and exit 0")
    check.add_argument("--no-recursive", dest="recursive", action="store_false")

    serve = sub.add_parser("serve", help="start the app and open it in a browser")
    serve.add_argument("folder", type=Path, nargs="+", help="library folder(s) to index")
    serve.add_argument("--port", type=int, default=8737)
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--no-recursive", dest="recursive", action="store_false")
    serve.add_argument("--no-browser", dest="browser", action="store_false",
                       help="do not open a browser tab")
    serve.add_argument("--cache-dir", type=Path, default=None,
                       help="where to keep the index (default: .ggasm/ beside the library)")
    serve.add_argument("--data-dir", type=Path, default=None,
                       help="where to keep shared curation (default: ggasm/ beside the library)")
    serve.add_argument("--reload", action="store_true",
                       help="restart the server when a .py file in the package changes")

    return parser


def cmd_check(args: argparse.Namespace) -> int:
    """Compare the library against a committed baseline; fail on anything new.

    Meant for a pre-commit hook on the plasmid repository. It deliberately does
    not fail on problems that were already there - a library of several hundred
    files always has some - only on ones this commit introduces. That is the
    difference between a hook people keep and a hook people disable.
    """
    library = Library(args.folder, recursive=args.recursive)
    library.scan()
    current = baseline.state(library)

    if args.update:
        args.baseline.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n",
                                 encoding="utf-8")
        print(f"baseline written to {args.baseline} "
              f"({len(current['conflicts'])} conflicts, "
              f"{len(current['unrecognised'])} unrecognised)")
        return 0

    if not args.baseline.exists():
        print(f"no baseline at {args.baseline}. Run `ggasm check --update` once to "
              f"record the current state, then commit it.", file=sys.stderr)
        return 2

    before = json.loads(args.baseline.read_text(encoding="utf-8"))
    problems = 0
    for kind, label in (("conflicts", "conflict"), ("unrecognised", "unrecognised")):
        known = set(before.get(kind, {}))
        for name in sorted(set(current[kind]) - known):
            print(f"new {label}: {name}: {current[kind][name]}", file=sys.stderr)
            problems += 1

    fixed = sum(
        len(set(before.get(kind, {})) - set(current[kind]))
        for kind in ("conflicts", "unrecognised")
    )
    if problems:
        print(f"\n{problems} new problem(s). Fix them, or re-baseline with "
              f"`ggasm check --update`.", file=sys.stderr)
        return 1

    print(f"{current['total']} plasmids, nothing new"
          + (f" ({fixed} fewer than the baseline)" if fixed else ""))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "scan":
        return cmd_scan(args)
    if args.command == "build":
        return cmd_build(args)
    if args.command == "check":
        return cmd_check(args)
    if args.command == "serve":
        return cmd_serve(args)
    parser.print_help()
    return 0


def cmd_scan(args: argparse.Namespace) -> int:
    for folder in args.folder:
        if not folder.is_dir():
            print(f"not a folder: {folder}", file=sys.stderr)
            return 2

    library = Library(args.folder, recursive=args.recursive)
    entries = library.scan(force=args.rescan)

    if args.part_type:
        entries = [e for e in entries if e.call.part_type == args.part_type]
    if args.unknown:
        entries = [e for e in entries if e.call.part_type is None]

    if args.json:
        print(json.dumps([e.to_dict() for e in entries], indent=2))
        return 0

    print_table(entries, library)
    return 0


COLUMNS = ("plasmid", "type", "conf", "5' -> 3'", "bp", "roles", "internal", "marker")


def print_table(entries: list[PlasmidEntry], library: Library) -> None:
    rows = [_row(e) for e in entries]
    widths = [
        max(len(header), *(len(row[i]) for row in rows)) if rows else len(header)
        for i, header in enumerate(COLUMNS)
    ]
    line = "  ".join(h.ljust(w) for h, w in zip(COLUMNS, widths, strict=True))
    print(line)
    print("  ".join("-" * w for w in widths))
    for row in rows:
        print("  ".join(cell.ljust(w) for cell, w in zip(row, widths, strict=True)))

    print()
    print(_summary(entries, library))
    for entry in entries:
        if entry.call.conflict:
            print(f"  conflict  {entry.name}: {entry.call.conflict}")
    for path, error in sorted(library.errors.items()):
        print(f"  skipped   {path}: {error}")


def _row(entry: PlasmidEntry) -> tuple[str, ...]:
    call = entry.call
    if call.five_prime and call.three_prime:
        overhangs = f"{call.five_prime} -> {call.three_prime}"
    elif entry.cassette_overhangs:
        five, three = entry.cassette_overhangs
        overhangs = f"{five} -> {three} (BsmBI)"
    else:
        overhangs = ""
    internal = ", ".join(
        f"{n}x{label}"
        for n, label in (
            (entry.sites.part_enzyme_internal, "BsaI"),
            (entry.sites.multigene_enzyme_internal, "BsmBI"),
            (entry.sites.linearizer_total, "NotI"),
        )
        if n
    )
    conf = call.confidence if call.source == "detected" else "manual"
    return (
        entry.name[:40],
        call.display_type,
        conf + ("!" if call.conflict else ""),
        overhangs,
        str(entry.length),
        ",".join(r for r in entry.roles if r != "part"),
        internal,
        entry.ecoli_marker or "",
    )


def _summary(entries: list[PlasmidEntry], library: Library) -> str:
    classified = sum(1 for e in entries if e.call.part_type)
    by_digest = sum(1 for e in entries if e.call.confidence == "digest")
    cassettes = sum(1 for e in entries if e.is_cassette)
    vectors = sum(1 for e in entries if e.is_multigene_vector)
    return (
        f"{len(entries)} plasmids, {classified} typed ({by_digest} by digest), "
        f"{cassettes} cassettes, {vectors} multigene vectors, "
        f"{len(library.errors)} unreadable"
    )


def cmd_build(args: argparse.Namespace) -> int:
    library = Library(args.folder)
    library.scan()

    design = level2_cassette.CassetteDesign(
        split_3=args.split_3,
        split_4=args.split_4,
        split_8=args.split_8,
        composite_left=args.composite_left,
        composite_right=args.composite_right,
        name=args.name,
    )
    active = level2_cassette.slots(design, library.scheme)
    if len(args.part) != len(active):
        expected = ", ".join(s.label for s in active)
        print(f"this design has {len(active)} slots ({expected}); "
              f"you gave {len(args.part)} part(s)", file=sys.stderr)
        return 2
    design.selections = {slot.key: name for slot, name in zip(active, args.part, strict=True)}

    result = level2_cassette.build(library, design)

    for slot in active:
        print(f"  {slot.label:<4} {slot.five_prime} -> {slot.three_prime}  {design.selections[slot.key]}")
    print()
    for issue in result.issues:
        print(f"  {issue}")

    if not result.ok:
        print(f"\nno product: {len(result.errors)} error(s)", file=sys.stderr)
        return 1

    print(f"\n{args.name}: {result.length} bp, {len(result.parts)} parts, "
          f"{len(result.junctions)} junctions, assembled with {result.enzyme}")
    for part in result.parts:
        print(f"  {part.start:>6}-{part.end:<6} {str(part.part_type):<5} {part.source_name}")

    if args.out:
        write_genbank(result.product, args.out)
        print(f"\nwrote {args.out}")
    return 0


def cmd_serve(args: argparse.Namespace) -> int:
    """Index the library and serve the app, optionally reloading on edits.

    ``--reload`` only watches ``*.py``, which is uvicorn's default and also the
    right answer here: the app serves ``web/`` with ``no-cache``, so a changed
    stylesheet or script is already one browser refresh away. Restarting the
    process for those would throw away the index and reload nothing new.

    Python is the opposite case. The routers register at ``create_app`` time,
    so a running server keeps serving the routes it started with no matter what
    the files on disk say - which looks exactly like a change that did not
    work.
    """
    import os
    import threading
    import webbrowser

    import uvicorn

    from .api.main import (
        ENV_CACHE_DIR, ENV_DATA_DIR, ENV_FOLDERS, ENV_RECURSIVE, PACKAGE, create_app,
    )

    for folder in args.folder:
        if not folder.is_dir():
            print(f"not a folder: {folder}", file=sys.stderr)
            return 2

    print(f"indexing {', '.join(str(f) for f in args.folder)} ...")

    # Several roots put the default state directories in whatever ancestor the
    # roots happen to share - which for a public library and a private one is
    # the home directory, or worse. Layering wants these said out loud.
    dirs = {"cache_dir": args.cache_dir, "data_dir": args.data_dir}
    if len(args.folder) > 1 and not (args.cache_dir and args.data_dir):
        print(f"  note: several roots, so state goes in "
              f"{Path(os.path.commonpath([str(f.resolve()) for f in args.folder]))}"
              f"; pass --cache-dir/--data-dir to put it somewhere deliberate")

    app = None
    if args.reload:
        # The app is built in the subprocess, not here, so index only far
        # enough to report the count and to fail now rather than inside a
        # reloader that would retry the same bad folder forever. The child
        # scans again, but against the mtime cache this leaves warm.
        library = Library(args.folder, recursive=args.recursive, **dirs)
        library.scan()
    else:
        app = create_app(args.folder, recursive=args.recursive, **dirs)
        library = app.state.library
    print(f"  {len(library.entries)} plasmids, {sum(1 for e in library.sorted_entries() if e.is_part)} usable parts")

    url = f"http://{args.host}:{args.port}/"
    if args.browser:
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    print(f"\nGG Assembler on {url}  (ctrl-c to stop)")

    if args.reload:
        os.environ[ENV_FOLDERS] = json.dumps([str(f.resolve()) for f in args.folder])
        os.environ[ENV_RECURSIVE] = "1" if args.recursive else "0"
        for key, value in ((ENV_CACHE_DIR, args.cache_dir), (ENV_DATA_DIR, args.data_dir)):
            os.environ[key] = str(value.resolve()) if value else ""
        print(f"  reloading on changes to {PACKAGE}/**/*.py")
        uvicorn.run(
            "ggassembler.api.main:from_env",
            factory=True,
            reload=True,
            reload_dirs=[str(PACKAGE)],
            host=args.host,
            port=args.port,
            log_level="warning",
        )
    else:
        uvicorn.run(app, host=args.host, port=args.port, log_level="warning")
    return 0
