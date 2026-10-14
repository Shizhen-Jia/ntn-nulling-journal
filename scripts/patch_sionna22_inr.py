#!/usr/bin/env python3
"""Apply the project's aggregate INR metric to official Sionna RT 2.2.0.

Run with the target environment's Python. Only two vendor files are changed;
their original bytes are backed up, and unexpected versions/edits are rejected.
Use --check to verify an installation without writing anything.
"""

import argparse
import hashlib
import importlib.metadata
from pathlib import Path
import tempfile


STOCK_SHA256 = {
    "radio_map.py": "5ef72086356857a325620711457c10044f4f2039f3cbb4d76a645de7b3f61072",
    "planar_radio_map.py": "6635acfd0521353d5e58858ffe2f7e201e1c8c83cbc1c616b3bd98e3fc51d313",
}
MARKER = "# NTN-Nulling aggregate INR patch for Sionna RT 2.2.0 (v1)\n"
INR = '("inr_ntn", "inr")'

PROPERTY = '''    @property
    def inr_ntn(self):
        r"""Aggregate terrestrial interference-to-noise ratio, in linear units.

        All transmitters contribute interference: sum(RSS, axis=0) / noise.
        The shape is ``[1, ...]`` for both planar and mesh radio maps,
        preserving the original NTN-Nulling project's singleton leading axis.
        This metric has no serving-transmitter association.
        """
        total_power = dr.sum(self.rss, axis=0)
        total_power = dr.reshape(mi.TensorXf, total_power.array,
                                 [1] + list(total_power.shape))
        return total_power / self._thermal_noise_power

    @property
    def inr(self):
        """Alias for :attr:`inr_ntn`."""
        return self.inr_ntn

'''


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise ValueError(f"Unexpected source layout for {old[:90]!r}")
    return source.replace(old, new, 1)


def patched_source(name, source):
    """Transform only the validated upstream 2.2.0 source."""
    source = MARKER + source
    for old, new in [
        ('["path_gain", "rss", "sinr"]',
         '["path_gain", "rss", "sinr", "inr_ntn", "inr"]'),
        ('("path_gain", "rss", "sinr")',
         '("path_gain", "rss", "sinr", "inr_ntn", "inr")'),
        ('["path_gain", "sinr"]', '["path_gain", "sinr", "inr_ntn", "inr"]'),
        ('One of ``"path_gain"``, ``"rss"``, or ``"sinr"``.',
         'One of ``"path_gain"``, ``"rss"``, ``"sinr"``, ``"inr_ntn"``, or ``"inr"``.'),
        ('if (tx is None) & (self.num_tx > 1):',
         f'if metric not in {INR} and tx is None and self.num_tx > 1:'),
        ('elif tx is not None:', f'elif metric not in {INR} and tx is not None:'),
    ]:
        source = source.replace(old, new)

    if name == "radio_map.py":
        source = replace_once(source, '    def tx_association(', PROPERTY + '    def tx_association(')
        source = replace_once(
            source,
            '        # No transmitter assignment for the cells with no coverage',
            f'''        if metric in {INR}:
            raise ValueError("Aggregate INR has no transmitter association; "
                             "use metric='rss' for strongest-source association.")

        # No transmitter assignment for the cells with no coverage''',
        )
        source = replace_once(
            source,
            '        cm = getattr(self, metric)\n        cm = dr.reshape',
            f'''        cm = getattr(self, metric)
        if metric in {INR}:
            if tx_association:
                raise ValueError("Aggregate INR sampling requires "
                                 "tx_association=False.")
            # Reuse the aggregate map for each TX's distance constraints.
            cm = dr.tile(cm.array, num_tx)
        cm = dr.reshape''',
        )
        association = '''        # Get cells for which the tx association is valid
        # [num_tx, cells_count]
        tx_ids = dr.arange(mi.UInt, num_tx)
        tx_ids = dr.reshape(mi.TensorXu, tx_ids, [num_tx, 1])
        tx_a = self.tx_association(metric)
        tx_a = dr.reshape(mi.TensorXu, tx_a, [1, cells_count])
        association_mask = tx_ids == tx_a

        # Compute combined mask
        # [num_tx, cells_count]
        active_cells = distance_mask & cm_mask
        if tx_association:
            active_cells = active_cells & association_mask'''
        source = replace_once(source, association, '''        # Compute combined mask; aggregate INR has no TX association.
        # [num_tx, cells_count]
        active_cells = distance_mask & cm_mask
        if tx_association:
            tx_ids = dr.arange(mi.UInt, num_tx)
            tx_ids = dr.reshape(mi.TensorXu, tx_ids, [num_tx, 1])
            tx_a = self.tx_association(metric)
            tx_a = dr.reshape(mi.TensorXu, tx_a, [1, cells_count])
            active_cells = active_cells & (tx_ids == tx_a)''')
        source = replace_once(
            source,
            '        else:\n            xlabel = "Signal-to-interference-plus-noise ratio (SINR) [dB]"',
            f'''        elif metric in {INR}:
            xlabel = "Interference-to-noise ratio (INR) [dB]"
            title = "Aggregate INR across all TXs"
        else:
            xlabel = "Signal-to-interference-plus-noise ratio (SINR) [dB]"''',
        )
        source = replace_once(source, '            tensor = tensor[tx]',
                              f'            tensor = tensor[0 if metric in {INR} else tx]')
        source = source.replace(
            '        accross the transmitters.\n',
            '        accross the transmitters.\n\n'
            '        For INR, a valid ``tx`` is accepted for API compatibility, but\n'
            '        the result always sums all transmitters.\n',
        )
        source = source.replace(
            '        :param num_cells: Number of returned random cells for each transmitter',
            '        Aggregate INR requires ``tx_association=False``. Its shared\n'
            '        map is filtered by each transmitter\'s distance constraints.\n\n'
            '        :param num_cells: Number of returned random cells for each transmitter',
        )
    else:
        source = replace_once(
            source,
            '        else:\n            colorbar_label = "Signal-to-interference-plus-noise ratio (SINR)"',
            f'''        elif metric in {INR}:
            colorbar_label = "Interference-to-noise ratio (INR) [dB]"
            title = "Aggregate INR across all TXs"
        else:
            colorbar_label = "Signal-to-interference-plus-noise ratio (SINR)"''',
        )
    compile(source, name, "exec")
    return source


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    dist = importlib.metadata.distribution("sionna-rt")
    if dist.version != "2.2.0":
        parser.error(f"Requires sionna-rt==2.2.0; found {dist.version}")
    root = Path(dist.locate_file("sionna/rt/radio_map_solvers"))
    changes = []
    # Validate all files before writing any file.
    for name, expected_hash in STOCK_SHA256.items():
        target = root / name
        backup = target.with_suffix(".py.pre-ntn-inr")
        current = target.read_bytes()
        original = backup.read_bytes() if backup.exists() else current
        if sha256(original) != expected_hash:
            raise RuntimeError(f"Unexpected upstream source: {target}; refusing to overwrite")
        patched = patched_source(name, original.decode()).encode()
        if current == patched:
            print(f"OK (already patched): {target}")
            continue
        if current != original:
            raise RuntimeError(f"Unrecognized edits in {target}; refusing to overwrite")
        changes.append((target, backup, original, patched))
    if args.check:
        if changes:
            parser.exit(1, "INR patch is not installed. Run without --check.\n")
        print("INR patch verified.")
        return
    for target, backup, original, patched in changes:
        if not backup.exists():
            with backup.open("xb") as f:
                f.write(original)
        # Replace the file rather than editing a possibly shared hardlink.
        with tempfile.NamedTemporaryFile(dir=root, delete=False) as f:
            f.write(patched)
            temporary = Path(f.name)
        temporary.chmod(target.stat().st_mode & 0o777)
        temporary.replace(target)
        print(f"Patched: {target}")
    print("INR patch installed. Restart any running notebook kernel.")


if __name__ == "__main__":
    main()
