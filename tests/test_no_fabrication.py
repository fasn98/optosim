"""Nothing that produces a scientific output may reach a random number generator.

The guard chemdisco uses, carried over for the same reason. A simulator is
entitled to be stochastic -- synaptic noise, channel noise, stochastic gating are
all legitimate -- but not silently. An RNG reached from a scoring or integration
path, with nothing declaring it, is how a trace acquires variability that its
provenance does not mention and that nobody can reproduce.

So: randomness is allowlisted per module, with a reason. Adding an import of
``random`` or ``numpy.random`` to a module not on the list fails the build. The
list is currently empty, which is the honest state of a deterministic simulator.
"""

from __future__ import annotations

import ast
import pathlib
import unittest

PACKAGE_ROOT = pathlib.Path(__file__).resolve().parent.parent / "optosim"

#: Modules permitted to reach a random number generator, each with the reason.
#: Empty by design: every model here is deterministic, so a seeded trace is
#: reproducible to the last digit and an unseeded one cannot exist.
RNG_ALLOWLIST: dict[str, str] = {}

RNG_MODULES = {"random", "numpy.random", "secrets", "scipy.stats"}
RNG_CALLS = {"rand", "randn", "random", "randint", "uniform", "normal", "choice", "shuffle", "seed"}


def python_files() -> list[pathlib.Path]:
    return sorted(PACKAGE_ROOT.rglob("*.py"))


def module_name(path: pathlib.Path) -> str:
    relative = path.relative_to(PACKAGE_ROOT.parent)
    return str(relative.with_suffix("")).replace("/", ".")


class TestNoHiddenRandomness(unittest.TestCase):
    def test_no_module_imports_an_rng_without_being_allowlisted(self) -> None:
        offenders: list[str] = []
        for path in python_files():
            name = module_name(path)
            tree = ast.parse(path.read_text())
            for node in ast.walk(tree):
                imported: set[str] = set()
                if isinstance(node, ast.Import):
                    imported = {alias.name for alias in node.names}
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imported = {node.module}
                    imported |= {f"{node.module}.{alias.name}" for alias in node.names}
                if imported & RNG_MODULES and name not in RNG_ALLOWLIST:
                    offenders.append(f"{name} imports {sorted(imported & RNG_MODULES)}")
        self.assertEqual(
            offenders,
            [],
            "a module reached an RNG without being allowlisted in RNG_ALLOWLIST "
            "with a reason. Stochasticity is allowed; undeclared stochasticity "
            "is not: " + "; ".join(offenders),
        )

    def test_no_module_calls_a_numpy_random_method(self) -> None:
        # Catches `np.random.normal(...)` where only `numpy` was imported, which
        # the import check above cannot see.
        offenders: list[str] = []
        for path in python_files():
            name = module_name(path)
            if name in RNG_ALLOWLIST:
                continue
            tree = ast.parse(path.read_text())
            for node in ast.walk(tree):
                if not isinstance(node, ast.Attribute):
                    continue
                if (
                    node.attr in RNG_CALLS
                    and isinstance(node.value, ast.Attribute)
                    and node.value.attr == "random"
                ):
                    offenders.append(f"{name} calls .random.{node.attr}")
        self.assertEqual(offenders, [], "; ".join(offenders))

    def test_the_allowlist_entries_all_carry_a_reason(self) -> None:
        for name, reason in RNG_ALLOWLIST.items():
            self.assertTrue(
                reason.strip(),
                f"{name} is allowlisted with no reason, which defeats the point "
                "of an allowlist",
            )

    def test_the_allowlist_only_names_modules_that_exist(self) -> None:
        # A stale entry would silently permit randomness in a module that was
        # renamed into something not on the list.
        existing = {module_name(path) for path in python_files()}
        for name in RNG_ALLOWLIST:
            self.assertIn(name, existing, f"allowlist names a module that no longer exists: {name}")


class TestDeterminism(unittest.TestCase):
    """The practical consequence: two identical runs agree to the last digit."""

    def test_two_identical_runs_are_bit_identical(self) -> None:
        from optosim.light import IrradianceProtocol
        from optosim.simulate import simulate_hh

        protocol = IrradianceProtocol.single(10.0, 20.0, 1.0)
        first = simulate_hh(protocol, duration_ms=50.0, dt_ms=0.02)
        second = simulate_hh(protocol, duration_ms=50.0, dt_ms=0.02)
        self.assertEqual(first.voltage.values, second.voltage.values)
        self.assertEqual(first.spike_times_ms, second.spike_times_ms)


class TestNoOutputClaimsToBeMeasured(unittest.TestCase):
    """A simulator that could label an output MEASURED would defeat the package."""

    def test_a_trace_refuses_a_non_predicted_origin(self) -> None:
        from optosim.provenance import Origin
        from optosim.simulate import Trace

        for origin in (Origin.MEASURED, Origin.DERIVED, Origin.HEURISTIC):
            with self.assertRaises(ValueError, msg=f"{origin} was accepted"):
                Trace(
                    name="V",
                    unit="mV",
                    times_ms=(0.0,),
                    values=(-65.0,),
                    model="test",
                    parameter_sets=("test",),
                    origin=origin,
                )

    def test_the_audit_passes_on_a_real_run(self) -> None:
        from optosim.light import IrradianceProtocol
        from optosim.report import provenance_audit
        from optosim.simulate import simulate_hh

        ok, message = provenance_audit(
            simulate_hh(IrradianceProtocol.single(5.0, 10.0, 1.0), duration_ms=30.0, dt_ms=0.05)
        )
        self.assertTrue(ok, message)

    def test_every_summary_quantity_is_predicted(self) -> None:
        from optosim.light import IrradianceProtocol
        from optosim.provenance import Origin
        from optosim.simulate import simulate_hh

        result = simulate_hh(
            IrradianceProtocol.single(5.0, 10.0, 1.0), duration_ms=30.0, dt_ms=0.05
        )
        for quantity in (
            result.voltage.peak(),
            result.voltage.trough(),
            result.open_fraction.peak(),
            result.photocurrent_pa.trough(),
            result.n_spikes,
            result.mean_rate_hz(),
        ):
            self.assertIs(quantity.origin, Origin.PREDICTED)
