"""GG Assembler - Golden Gate assembly design for the Yeast Toolkit (YTK) standard.

Layering (enforced by tests/test_layering.py):

    core/    pure logic; knows nothing about levels, the API or the web front end
    levels/  one module per assembly level; may use core, never a sibling level
    api/     FastAPI routes; routes_levelN touches only levels/levelN_* and core
    web/     static front end; one self-contained module directory per level
"""

__version__ = "0.1.0"
