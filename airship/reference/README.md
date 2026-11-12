# reference/

`Final.py` is the conference-paper code, copied unchanged from Mana's upload.
It must never be edited, so the paper's results stay reproducible.

- SHA-256: `15075bb53b2f431a01b9e42c132b9c7f585fa6592bebd2bc33de9c278dc9f650`
- Settings: the paper's cost function (line 141 `EXCLUDED_OBJ_INDICES = (0, 1, 6, 8)`,
  so direction symmetry is on), 50 m path, `SIM_DT = 0.05` s. The header comment
  still says J6 is excluded; that comment is out of date, line 141 is what runs.
- The brief mentions `Final2.py`; Phase 0 was built from this `Final.py`
  because it reproduces the paper's layouts (N = 16 / 13 / 16).
