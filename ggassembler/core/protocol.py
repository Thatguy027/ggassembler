"""Turn a finished design into a reaction you can pipette.

A Golden Gate reaction wants its DNA **equimolar**, not equal by mass: a 700 bp
promoter and a 10 kb backbone at the same ng/µL differ fourteen-fold in copy
number, and the short piece wins the ligation. So every calculation here starts
from femtomoles and converts to a volume through the one number that cannot be
read out of a GenBank file - what the prep on your bench measured at.

That is why `Library.set_concentration` exists and why a component with no
concentration comes back with its volume left blank and an issue attached,
rather than with a number quietly invented for it. A protocol that looks
complete but silently assumed 100 ng/µL is worse than one that says what it
is missing.

Molar mass follows the usual double-stranded estimate, 617.96 g/mol per base
pair plus 36.04 for the ends, which is accurate enough that the 20 fmol target
is the dominant uncertainty.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .enzymes import Enzyme

#: Average molar mass of one base pair of double-stranded DNA, g/mol.
BP_MASS = 617.96
#: The two hydroxyls and one proton the polymer carries at its ends, g/mol.
END_MASS = 36.04

#: What each piece of DNA contributes, in femtomoles. The YTK paper's one-pot
#: reactions are equimolar; 20 fmol in a 10 µL reaction is the usual working
#: point and leaves room for a dozen components.
DEFAULT_FMOL = 20.0
#: Total reaction volume, µL.
DEFAULT_VOLUME = 10.0

#: Fixed additions, in µL per `DEFAULT_VOLUME`, with what each is for.
REAGENTS: tuple[tuple[str, float, str], ...] = (
    ("T4 DNA ligase buffer, 10×", 1.0, "supplies the ATP the ligase needs"),
    ("T4 DNA ligase (400 U/µL)", 0.5, "seals each junction as it forms"),
)

#: Cutting temperature by enzyme. Everything else about the cycle is shared.
CUT_TEMPERATURE = {"BsaI": 37, "BsmBI": 42, "Esp3I": 42, "BbsI": 37}
DEFAULT_CUT_TEMPERATURE = 37

#: Cycles of cut-and-ligate. The YTK paper uses 25 for a full eight-part
#: assembly; fewer is enough for two or three pieces but costs nothing.
CYCLES = 25


@dataclass(frozen=True)
class Component:
    """One line of the reaction table."""

    name: str
    kind: str
    """``part``, ``destination``, ``reagent`` or ``water``."""
    path: str = ""
    """Where the plasmid is filed, so a caller can record a concentration
    against it without parsing the display name back apart."""
    length: int | None = None
    conc_ng_ul: float | None = None
    fmol: float | None = None
    ng: float | None = None
    volume_ul: float | None = None
    note: str = ""

    @property
    def known(self) -> bool:
        return self.volume_ul is not None


@dataclass(frozen=True)
class Step:
    """One line of the thermocycler programme."""

    label: str
    detail: str


@dataclass
class Reaction:
    """A complete, pipettable setup - or an honest account of what is missing."""

    name: str
    enzyme: str
    total_ul: float
    fmol_each: float
    components: list[Component] = field(default_factory=list)
    steps: list[Step] = field(default_factory=list)
    issues: list[str] = field(default_factory=list)
    selection: str = ""

    @property
    def ok(self) -> bool:
        """True when every DNA component could be turned into a volume."""
        return not self.issues

    @property
    def dna_ul(self) -> float:
        return sum(c.volume_ul or 0.0 for c in self.components if c.kind in ("part", "destination"))

    @property
    def missing(self) -> list[str]:
        return [c.name for c in self.components
                if c.kind in ("part", "destination") and not c.known]


def molecular_weight(length_bp: int) -> float:
    """Molar mass of a double-stranded molecule of `length_bp`, g/mol."""
    return length_bp * BP_MASS + END_MASS


def ng_for(fmol: float, length_bp: int) -> float:
    """Nanograms of a `length_bp` molecule that come to `fmol` femtomoles."""
    return fmol * molecular_weight(length_bp) * 1e-6


def fmol_for(ng: float, length_bp: int) -> float:
    """The reverse: femtomoles in `ng` nanograms of a `length_bp` molecule."""
    return ng / (molecular_weight(length_bp) * 1e-6)


def volume_for(ng: float, conc_ng_ul: float) -> float:
    """Microlitres of a `conc_ng_ul` prep that deliver `ng` nanograms."""
    if conc_ng_ul <= 0:
        raise ValueError("a concentration must be greater than zero")
    return ng / conc_ng_ul


def _round(value: float) -> float:
    """To the nearest 0.01 µL - finer than anything you can pipette anyway."""
    return round(value + 0.0, 2)


def _dna_component(
    name: str,
    length: int,
    conc_ng_ul: float | None,
    fmol: float,
    kind: str,
    path: str = "",
) -> Component:
    ng = ng_for(fmol, length)
    if not conc_ng_ul or conc_ng_ul <= 0:
        return Component(
            name=name, kind=kind, path=path, length=length, conc_ng_ul=None,
            fmol=fmol, ng=_round(ng), volume_ul=None,
            note="no concentration recorded - measure the prep and enter it in the Library",
        )
    return Component(
        name=name, kind=kind, path=path, length=length, conc_ng_ul=conc_ng_ul,
        fmol=fmol, ng=_round(ng), volume_ul=_round(volume_for(ng, conc_ng_ul)),
        note=f"{_round(ng)} ng at {conc_ng_ul:g} ng/µL",
    )


def cycling(enzyme: Enzyme | str, cycles: int = CYCLES) -> list[Step]:
    """The thermocycler programme for a one-pot digest-ligate reaction.

    The cutting temperature is the enzyme's own; everything else is the same
    whichever Type IIS enzyme is driving the assembly.
    """
    name = enzyme if isinstance(enzyme, str) else enzyme.name
    cut = CUT_TEMPERATURE.get(name, DEFAULT_CUT_TEMPERATURE)
    return [
        Step(f"{cycles}× cycle",
             f"{cut} °C 5 min → 16 °C 5 min "
             f"({name} cuts, ligase seals)"),
        Step("Final digest",
             f"{cut} °C 10 min (linearises anything that re-formed uncut)"),
        Step("Inactivate", "80 °C 10 min"),
        Step("Hold", "4 °C"),
    ]


def reaction(
    name: str,
    enzyme: Enzyme | str,
    parts: list[tuple[str, int, float | None] | tuple[str, int, float | None, str]],
    destination: tuple[str, int, float | None] | tuple[str, int, float | None, str] | None = None,
    fmol_each: float = DEFAULT_FMOL,
    total_ul: float = DEFAULT_VOLUME,
    selection: str = "",
) -> Reaction:
    """Build the reaction table for one assembly.

    `parts` and `destination` are ``(name, length_bp, conc_ng_ul)``; a `None`
    concentration is carried through as a gap rather than guessed at.
    """
    enzyme_name = enzyme if isinstance(enzyme, str) else enzyme.name
    out = Reaction(
        name=name,
        enzyme=enzyme_name,
        total_ul=total_ul,
        fmol_each=fmol_each,
        selection=selection,
        steps=cycling(enzyme_name),
    )

    if destination:
        label, length, conc, *rest = destination
        out.components.append(
            _dna_component(label, length, conc, fmol_each, "destination", rest[0] if rest else "")
        )
    for label, length, conc, *rest in parts:
        out.components.append(
            _dna_component(label, length, conc, fmol_each, "part", rest[0] if rest else "")
        )

    scale = total_ul / DEFAULT_VOLUME
    for label, volume, why in REAGENTS:
        out.components.append(
            Component(name=label, kind="reagent", volume_ul=_round(volume * scale), note=why)
        )
    out.components.append(
        Component(
            name=f"{enzyme_name} (10 U/µL)",
            kind="reagent",
            volume_ul=_round(0.5 * scale),
            note="cuts every site, releasing the fragments that ligate",
        )
    )

    used = sum(c.volume_ul or 0.0 for c in out.components)
    water = total_ul - used
    out.components.append(
        Component(
            name="Nuclease-free water",
            kind="water",
            volume_ul=_round(max(water, 0.0)),
            note=f"to {total_ul:g} µL",
        )
    )

    missing = out.missing
    if missing:
        out.issues.append(
            f"{len(missing)} component(s) have no recorded concentration "
            f"({', '.join(missing[:4])}{'…' if len(missing) > 4 else ''}): "
            "enter ng/µL on the Library screen and the volumes fill in."
        )
    if water < 0 and not missing:
        out.issues.append(
            f"the DNA alone comes to {_round(used)} µL, more than the "
            f"{total_ul:g} µL reaction. Dilute the most concentrated preps, "
            f"or scale the reaction up."
        )
    return out


def as_text(rx: Reaction) -> str:
    """The reaction as a plain-text block, for a notebook or a protocol sheet."""
    width = max((len(c.name) for c in rx.components), default=10)
    lines = [
        f"{rx.name} - {rx.enzyme} Golden Gate assembly",
        f"{'=' * (len(rx.name) + len(rx.enzyme) + 26)}",
        "",
        f"{rx.fmol_each:g} fmol of each piece in {rx.total_ul:g} uL.",
        "",
        f"{'Component'.ljust(width)}  {'uL':>6}  {'ng':>8}  Notes",
        f"{'-' * width}  {'-' * 6}  {'-' * 8}  {'-' * 40}",
    ]
    for c in rx.components:
        volume = f"{c.volume_ul:.2f}" if c.volume_ul is not None else "  ?  "
        mass = f"{c.ng:.1f}" if c.ng is not None else ""
        lines.append(f"{c.name.ljust(width)}  {volume:>6}  {mass:>8}  {c.note}")

    lines += ["", "Thermocycler", "-" * 12]
    for step in rx.steps:
        lines.append(f"  {step.label}: {step.detail}")

    if rx.selection:
        lines += ["", f"Selection: {rx.selection}"]
    if rx.issues:
        lines += ["", "Before you set this up", "-" * 22]
        lines += [f"  ! {issue}" for issue in rx.issues]
    return "\n".join(lines) + "\n"
