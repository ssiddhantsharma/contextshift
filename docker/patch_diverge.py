"""Make DIVERGE 4.1.0 buildable and Type-I usable.

Three upstream problems, each with a PR open:
  setup.py compiles src/, which is MSVC-only            -> zjupgx/diverge4#10
  `#define version` collides with a pybind11 member     -> zjupgx/diverge4#10
  Gu99/Rvs/TypeOneAnalysis call an undefined get_colnames -> zjupgx/diverge4#8
"""

import sys
from pathlib import Path

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")

setup = root / "setup.py"
s = setup.read_text()
if "src_linux/{module}" not in s:
    s = s.replace("glob(f'src/{module}/*.c*')", "glob(f'src_linux/{module}/*.c*')")
    setup.write_text(s)
    print("setup.py -> src_linux")

macro = b'#define version "V1.0"\n'
for p in sorted((root / "src_linux").rglob("common.h")):
    b = p.read_bytes()
    if macro in b:
        p.write_bytes(b.replace(macro, b""))
        print(f"dropped version macro: {p.relative_to(root)}")

binding = root / "diverge" / "binding.py"
b = binding.read_text()
if "def get_colnames" not in b:
    binding.write_text(
        b.replace("def load_tree_file(",
                  "def get_colnames(names):\n    return list(names)\n\n\ndef load_tree_file(", 1))
    print("defined get_colnames")
