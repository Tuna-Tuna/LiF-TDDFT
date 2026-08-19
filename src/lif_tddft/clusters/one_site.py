"""Deterministic one-active-site C9/C25 topology models."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Atom:
    species: str
    coordinates_bohr: tuple[float, float, float]


@dataclass(frozen=True)
class ClusterModel:
    name: str
    extent: str
    active_site_id: str
    active_site_coordinates_bohr: tuple[float, float, float]
    path_direction: str
    explicit_collision_sites: int
    soft_potential_gamma: float
    replace_path_anions_with_soft_potential: bool
    atoms: tuple[Atom, ...]

    def octopus_coordinate_block(self) -> str:
        return "\n".join(
            f' "{atom.species}" | {atom.coordinates_bohr[0]:.8f} | '
            f'{atom.coordinates_bohr[1]:.8f} | {atom.coordinates_bohr[2]:.8f}'
            for atom in self.atoms
        )


def build_one_active_site_cluster(
    extent: str,
    active_site_id: str = "central_surface_F",
    path_direction: str = "[001]",
    replace_path_anions_with_soft_potential: bool = True,
) -> ClusterModel:
    if extent not in {"small", "large"}:
        raise ValueError("extent must be 'small' or 'large'")
    width = 3 if extent == "small" else 5
    center = width // 2
    spacing = 3.7958
    atoms: list[Atom] = []
    for ix in range(width):
        for iy in range(width):
            x = (ix - center) * spacing
            y = (iy - center) * spacing
            surface_is_f = (ix + iy) % 2 == 0
            is_active = ix == center and iy == center
            if surface_is_f:
                species = "F" if is_active else (
                    "udfF" if replace_path_anions_with_soft_potential else "F"
                )
            else:
                species = "Li"
            atoms.append(Atom(species, (x, y, 0.0)))
            atoms.append(Atom("Li" if surface_is_f else "F", (x, y, -spacing)))
    return ClusterModel(
        name="C9_one_site" if extent == "small" else "C25_one_site",
        extent=extent,
        active_site_id=active_site_id,
        active_site_coordinates_bohr=(0.0, 0.0, 0.0),
        path_direction=path_direction,
        explicit_collision_sites=1,
        soft_potential_gamma=0.566,
        replace_path_anions_with_soft_potential=replace_path_anions_with_soft_potential,
        atoms=tuple(atoms),
    )


def topology_check(small: ClusterModel, large: ClusterModel) -> dict[str, object]:
    result = {
        "explicit_collision_sites": small.explicit_collision_sites,
        "active_site_coordinates_match": small.active_site_coordinates_bohr == large.active_site_coordinates_bohr,
        "trajectory_matches": small.path_direction == large.path_direction,
        "soft_potential_definition_matches": (
            small.soft_potential_gamma == large.soft_potential_gamma
            and small.replace_path_anions_with_soft_potential == large.replace_path_anions_with_soft_potential
        ),
        "expected_atom_counts": len(small.atoms) == 18 and len(large.atoms) == 50,
    }
    result["passed"] = (
        result["explicit_collision_sites"] == 1
        and all(bool(result[key]) for key in result if key not in {"explicit_collision_sites", "passed"})
    )
    return result
