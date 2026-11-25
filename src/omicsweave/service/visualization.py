"""在既有结果包上共享通路图与摘要图渲染。"""
from __future__ import annotations

from pathlib import Path

from .schemas import ErrorCode, ServiceError


def _bundle(project, experiment: str, source: str | None) -> Path:
    from omicsweave.contract.export import find_run

    root = project.path("results")
    if source:
        return find_run(root, experiment, source)
    base = root / experiment
    results = sorted(path.parents[3] for path in base.glob(
        "run_*/integration/reaction_dossier/tables/reaction_assessment.parquet",
    ))
    if not results:
        raise FileNotFoundError(f"{base} 下没有含反应档案评估表的整合结果包")
    return results[-1]


def render_map(root: Path, experiment: str, source: str | None, output: Path | None, settings: dict) -> dict:
    from omicsweave.core.config import load_experiment, load_project
    from omicsweave.integration.pathway_map import redraw

    project = load_project(root)
    exp = load_experiment(project, experiment)
    bundle = _bundle(project, experiment, source)
    configured = dict((exp.config.integration.parameters.get("reaction_dossier") or {}).get("pathway_map") or {})
    for name, value in settings.items():
        if isinstance(value, dict) and value.get("enabled") is not False and isinstance(configured.get(name), dict):
            configured[name] = {**configured[name], **value}
        else:
            configured[name] = value
    configured["enabled"] = True
    directory, notes = redraw(bundle, output, configured, exp.root)
    return {"output_dir": str(directory), "notes": notes,
            "files": [str(path) for path in sorted(directory.rglob("*")) if path.is_file()]}


def render_abstract(root: Path, *, experiment: str | None = None, source: str | None = None,
                    spec_path: Path | None = None, output: Path | None = None, style: str = "metro_atlas",
                    color_mode: str = "semantic", journal: str = "cell_press", layout: str | None = None,
                    formats: list[str] | None = None, tissues: list[str] | None = None,
                    contrasts: list[str] | None = None, layout_path: str = "") -> dict:
    formats = formats or ["pdf", "png"]
    if spec_path is not None:
        from omicsweave.visualization.graphical_abstract import load_abstract_spec, render_graphical_abstract

        spec = load_abstract_spec(spec_path)
        stem = spec_path.with_name(f"{spec_path.stem}_{style}")
        if output is not None:
            stem = (output / stem.name) if output.is_dir() or not output.suffix else output.with_suffix("")
        written = render_graphical_abstract(spec, stem, style=style, color_mode=color_mode,
                                           journal=journal, layout=layout, formats=formats)
        return {"files": [str(path) for path in written], "style": style}
    if not experiment:
        raise ServiceError(ErrorCode.CLI_USAGE_ERROR, "请提供 --spec <YAML> 或 --experiment <实验ID>")

    from omicsweave.core.config import load_experiment, load_project
    from omicsweave.integration.pathway_map import render_abstracts

    project = load_project(root)
    exp = load_experiment(project, experiment)
    directory, written = render_abstracts(
        _bundle(project, experiment, source), output, tissues=tissues, contrasts=contrasts,
        style=style, color_mode=color_mode, journal=journal, ga_layout=layout, formats=formats,
        layout_path=layout_path, root=exp.root,
    )
    return {"output_dir": str(directory), "files": [str(path) for path in written], "style": style}
