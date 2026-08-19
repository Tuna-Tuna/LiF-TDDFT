"""Render Octopus inputs and enforce prescribed-path physical gates."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from .campaign import RunSpec


def render_input(
    template_path: str | Path,
    spec: RunSpec,
    config: dict[str, Any],
    *,
    relative_velocity_au: float | None = None,
    from_scratch: str = "no",
    move_ions: bool = True,
    constant_velocity: bool = False,
    calculation_mode: str = "td",
    propagation_time_expression: str | None = None,
) -> str:
    numerics = config["numerics"]
    velocity = float(spec.relative_velocity_au if relative_velocity_au is None else relative_velocity_au)
    pf = float(config.get("projectile_fraction", -0.5))
    cf = float(config.get("cluster_fraction", 0.5))
    outputs = list(config.get("outputs", []))
    spatial_outputs = [item for item in outputs if item in {"density", "current"}]
    td_outputs = [item for item in outputs if item not in spatial_outputs]
    projectile_velocity = -velocity if spec.variant == "isolated_projectile" else pf * velocity
    active_site_x_bohr = -0.8 * 3.7958
    initial_separation = float(config["initial_projectile_x_bohr"])
    final_separation = float(config["final_separation_bohr"])
    context = dict(
        run_id=spec.run_id,
        cluster_model=spec.cluster_model,
        variant=spec.variant,
        from_scratch=from_scratch,
        calculation_mode=calculation_mode,
        mpi_ranks=int(numerics["mpi_ranks"]),
        box_radius_angstrom=float(numerics["box_radius_angstrom"]),
        grid_spacing_angstrom=float(numerics["grid_spacing_angstrom"]),
        move_ions="yes" if move_ions else "no",
        ions_constant_velocity="yes" if constant_velocity else "no",
        height_bohr=spec.height_bohr,
        relative_velocity_au=velocity,
        projectile_velocity_au=projectile_velocity,
        cluster_velocity_au=cf * velocity,
        projectile_x_bohr=active_site_x_bohr + initial_separation,
        propagation_distance_bohr=initial_separation + final_separation,
        cap_strength=float(numerics["cap_strength"]),
        cap_start_angstrom=float(numerics["cap_start_angstrom"]),
        cap_end_angstrom=float(numerics["cap_end_angstrom"]),
        output_interval=int(numerics["output_interval"]),
        time_step_expression=str(numerics["time_step_expression"]),
        propagation_time_expression=(
            propagation_time_expression
            if propagation_time_expression is not None
            else ("0" if calculation_mode == "gs" else f"({initial_separation + final_separation}/Fv)")
        ),
        xc_functional=str(config["xc_functional"]),
        spatial_outputs=spatial_outputs,
        td_outputs=td_outputs,
    )
    source = Path(template_path).read_text(encoding="utf-8")
    try:
        from jinja2 import Environment, StrictUndefined

        env = Environment(undefined=StrictUndefined, keep_trailing_newline=True)
        text = env.from_string(source).render(**context)
    except ImportError:
        text = _minimal_jinja_render(source, context)
    expected_velocity_lines = 1 if spec.variant == "isolated_projectile" else 19
    assert_input_contract(
        text,
        spec.height_bohr,
        velocity,
        move_ions,
        constant_velocity,
        expected_velocity_lines=expected_velocity_lines,
    )
    return text


def _minimal_jinja_render(source: str, context: dict[str, Any]) -> str:
    """Render the small variable/for-loop subset used by our Octopus template."""
    loop_pattern = re.compile(
        r"{%\s*for\s+(\w+)\s+in\s+(\w+)\s*%}(.*?){%\s*endfor\s*%}", re.DOTALL
    )

    def expand_loop(match: re.Match[str]) -> str:
        item_name, sequence_name, body = match.groups()
        if sequence_name not in context:
            raise ValueError(f"undefined template sequence: {sequence_name}")
        return "".join(
            re.sub(r"{{\s*" + re.escape(item_name) + r"\s*}}", str(item), body)
            for item in context[sequence_name]
        )

    rendered = loop_pattern.sub(expand_loop, source)

    def substitute(match: re.Match[str]) -> str:
        name = match.group(1)
        if name not in context:
            raise ValueError(f"undefined template variable: {name}")
        return str(context[name])

    return re.sub(r"{{\s*(\w+)\s*}}", substitute, rendered)


def _assignment(text: str, name: str) -> str:
    match = re.search(rf"^\s*{re.escape(name)}\s*=\s*([^#\n]+)", text, re.MULTILINE)
    if not match:
        raise ValueError(f"missing Octopus assignment: {name}")
    return match.group(1).strip()


def assert_input_contract(
    text: str,
    height_bohr: float,
    relative_velocity_au: float,
    move_ions: bool = True,
    constant_velocity: bool = False,
    expected_velocity_lines: int = 19,
) -> None:
    expected_move = "yes" if move_ions else "no"
    if _assignment(text, "MoveIons").lower() != expected_move:
        raise ValueError("MoveIons does not match stage definition")
    expected = "yes" if constant_velocity else "no"
    if _assignment(text, "IonsConstantVelocity").lower() != expected:
        raise ValueError("IonsConstantVelocity does not match stage definition")
    if abs(float(_assignment(text, "height")) - height_bohr) > 1.0e-10:
        raise ValueError("rendered height does not match campaign")
    if abs(float(_assignment(text, "Fv")) - relative_velocity_au) > 1.0e-10:
        raise ValueError("rendered relative velocity does not match campaign")
    velocity_lines = re.findall(
        r'^\s*"[^\"]+"\s*\|\s*([^|\n]+)\|\s*([^|\n]+)\|\s*([^|\n]+)$',
        text,
        re.MULTILINE,
    )
    if len(velocity_lines) < expected_velocity_lines:
        raise ValueError("velocity block is incomplete")
    for _, vy, vz in velocity_lines[-expected_velocity_lines:]:
        if abs(float(vy)) > 1.0e-12 or abs(float(vz)) > 1.0e-12:
            raise ValueError("prescribed flyby requires zero y/z velocity")
