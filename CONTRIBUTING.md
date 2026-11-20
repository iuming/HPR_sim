# Contributing to hpr-sim

Thanks for considering a contribution! This is a research-grade, open
simulation tool for the SRF community; every improvement matters.

## Ways to contribute

- **Bug reports** — open an issue with the config, the version and the full
  CLI output.
- **New cavity cases** — a YAML + (small) mesh or a STEP path with the
  relevant mechanical parameters is a perfect contribution.
- **Physics improvements** — better drag correlations, droplet-size
  distributions, secondary-breakup models, or a wall-particle adhesion
  sub-model.
- **CI / packaging / docs** — always welcome.

## Development workflow

```bash
pip install -e . pytest
python -m pytest tests            # 24 tests, must stay green
```

Code conventions: PEP 8, type hints, no third-party dependencies for the
core `hpr_sim` package beyond numpy/matplotlib/PyYAML (tools may use more).

## Commit messages

Keep them imperative and short; reference the issue when one exists.

## Release checklist (maintainers)

1. `pytest` green on 3.9–3.12 (CI).
2. `examples/tesla_9cell.yaml` reproduces a ≤ ~5 % coverage change from the
   README value.
3. `CHANGELOG.md` updated; version bumped in `pyproject.toml` and
   `hpr_sim/__init__.py`.
4. Tag the commit (`git tag v<ver>`).

## Code of conduct

All participants must follow [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).
