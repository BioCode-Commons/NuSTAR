#!/usr/bin/env python3
"""Compile public headers independently and enforce module dependency boundaries."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
import subprocess
import sys

# Only runtime assembles the entire application. Algorithm modules consume
# narrow value/interface headers, never runtime or another layout's pipeline.
DEPENDENCIES = {
    "common": set(),
    "reference": {"common"},
    "annotation": {"common", "reference"},
    "index": {"common", "reference", "third_party"},
    "reads": {"common"},
    "seed_search": {"common", "index", "reads", "reference"},
    "window": {"common", "reads", "reference", "seed_search"},
    "stitching": {"common", "reads", "reference", "window"},
    "selection": {"common", "reads", "stitching"},
    "sam": {"common", "reads", "reference", "selection", "stitching"},
    "splice_junction": {"common", "reads", "reference", "selection", "stitching"},
    "statistics": {"common", "reads", "sam", "selection", "stitching"},
    "output": {"common", "reference", "sam", "statistics"},
    "pipeline": {"common", "reads", "reference", "sam", "seed_search", "selection",
                 "splice_junction", "statistics", "stitching", "window"},
    "bam": {"common", "sam", "reads", "reference", "output", "statistics", "third_party"},
    "solo": {"common", "reads", "reference", "selection", "stitching", "statistics", "pipeline", "sam", "third_party"},
    "profiling": set(),
    "cli": {"runtime", "reference", "index", "solo", "sam"},
}
DEPENDENCIES["runtime"] = set(DEPENDENCIES) - {"cli"}


def main() -> int:
    compiler, root_arg = sys.argv[1:3]
    project = Path(root_arg).resolve()
    headers = sorted(p for module in DEPENDENCIES
                     for p in (project / module).rglob("*.hpp"))
    files = sorted(p for module in (*DEPENDENCIES, "tests", "tools")
                   for p in (project / module).rglob("*")
                   if p.suffix in {".hpp", ".cpp"})
    errors = []
    graph = {}
    for path in files:
        text = path.read_text()
        rel = path.relative_to(project)
        if re.search(r'#include [<"](?:nustar/|STAR|[^"\n]*core_contracts)', text):
            errors.append(f"{rel}: forbidden external/legacy include")
        targets = set()
        for include in re.findall(r'^\s*#include\s+"([^"]+)"', text, re.MULTILINE):
            target = (path.parent / include).resolve()
            if Path(include).is_absolute() or not target.is_relative_to(project) or not target.is_file():
                errors.append(f"{rel}: missing or external relative include {include}")
                continue
            targets.add(target)
            source_module, target_module = rel.parts[0], target.relative_to(project).parts[0]
            if (source_module in DEPENDENCIES and target_module != source_module
                    and target_module not in DEPENDENCIES[source_module]):
                errors.append(f"{rel}: forbidden dependency on {target.relative_to(project)}")
        graph[path] = targets

    def closure(path):
        seen, pending = set(), [project / path]
        while pending:
            current = pending.pop()
            if current not in seen:
                seen.add(current)
                pending.extend(graph.get(current, ()))
        return {p.relative_to(project).as_posix() for p in seen}

    for entry in closure("window/interface.hpp"):
        if entry.startswith(("index/", "runtime/")) or entry == "seed_search/plan.hpp":
            errors.append(f"window/interface.hpp depends on forbidden {entry}")
    for header in ("pipeline/single_end/pipeline.hpp", "runtime/single_runner.hpp"):
        for entry in closure(header):
            if entry.startswith("pipeline/paired_end/") or entry in {
                    "runtime/paired_plan.hpp", "runtime/paired_runner.hpp"}:
                errors.append(f"{header} depends on forbidden {entry}")
    for entry in closure("pipeline/paired_end/overlap.hpp"):
        if entry in {"pipeline/single_end/pipeline.hpp", "runtime/single_runner.hpp"}:
            errors.append(f"paired overlap depends on forbidden {entry}")

    def compile_header(header):
        relative = header.relative_to(project).as_posix()
        result = subprocess.run(
            [compiler, "-std=c++20", "-fsyntax-only", "-x", "c++", "-I", str(project), "-"],
            input=f'#include "{relative}"\nint main() {{}}\n',
            text=True, capture_output=True, check=False)
        return f"{relative}:\n{result.stderr}" if result.returncode else None

    with ThreadPoolExecutor(max_workers=4) as pool:
        errors.extend(error for error in pool.map(compile_header, headers) if error)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"Checked {len(headers)} headers, {len(files)} source files and module/SE/PE/window boundaries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
