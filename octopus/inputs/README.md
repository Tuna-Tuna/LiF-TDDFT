# Octopus inputs

Generate inputs after filling the required height nodes, CAP magnitudes and
actual pseudopotential filenames in `config/campaign_revision.yaml`.
The manuscript specifies the height range, not every calculated node or the
CAP magnitude for each velocity, so no generated input is published with
guessed values. Run `python octopus/generate_reference_inputs.py` for the
v=0.30 examples after completing those entries, or use `generate_inp_files.py`
for the full paired campaign. The local production stage prescribes the
projectile velocity and keeps all surface ions fixed.
