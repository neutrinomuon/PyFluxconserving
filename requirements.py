# save as generate_requirements.py at the root
import os
import ast
from pathlib import Path

project_dirs = ["src/python", "src/fortran"]  # directories to scan

all_imports = set()

for d in project_dirs:
    for py_file in Path(d).rglob("*.py"):
        try:
            with open(py_file, "r", encoding="utf-8") as f:
                tree = ast.parse(f.read(), filename=str(py_file))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for n in node.names:
                        all_imports.add(n.name.split('.')[0])
                elif isinstance(node, ast.ImportFrom):
                    if node.module:
                        all_imports.add(node.module.split('.')[0])
        except Exception as e:
            print(f"Skipping {py_file}: {e}")

# Exclude standard library modules
import sys
import importlib.util
stdlib_modules = set(sys.builtin_module_names)

# Minimal filter: remove standard lib modules
requirements = sorted(m for m in all_imports if m not in stdlib_modules)

# Save to requirements.txt
with open("requirements.txt", "w") as f:
    for r in requirements:
        f.write(r + "\n")

print("Generated requirements.txt:")
print("\n".join(requirements))
