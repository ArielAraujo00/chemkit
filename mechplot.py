from dataclasses import dataclass, field

@dataclass(slots=True)
class Step:
    rc: int
    label: str
    codes: list[str]
    step_type: str = 'INT'
    color: str | None = None
    ddG: float | None = None


@dataclass(slots=True)
class Block:
    name: str
    steps: list[Step]
    color: str = 'black'
    linestyle: str = '-'
    linewidth: float = 2.0
    visible: bool = True
    def compute_energies(self, df):
        for step in self.steps:
            G = sum(
                df.loc[df['code'] == c, 'G'].iloc[0]
                for c in step.codes
            )
            step.ddG = G


@dataclass(slots=True)
class ReactionProfile:

    blocks: list[Block]

    connections: list[tuple[str, str]] = field(default_factory=list)

    def compute_relative_energies(self, reference_rc=0):

        all_steps = [
            step
            for block in self.blocks
            for step in block.steps
        ]

        G_ref = min(
            step.ddG
            for step in all_steps
            if step.rc == reference_rc
        )

        for step in all_steps:
            step.ddG -= G_ref