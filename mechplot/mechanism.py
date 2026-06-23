from dataclasses import dataclass, field
import pandas as pd

from chemkit.mechplot.styles import *

@dataclass(slots=True)
class Step:
    id: str
    rc: float
    label: str
    codes: list[str]
    style: StepStyle | None = None
    block: str | None = None
    x: float | None = None
    visible: bool = True
    def get_x(self) -> float:
        return self.rc if self.x is None else self.x

@dataclass(slots=True)
class Block:
    name: str
    step_ids: list[str]
    style: BlockStyle = field(default_factory=BlockStyle)
    @classmethod
    def from_steps(cls, name: str, steps: list[Step], **kwargs):
        for step in steps:
            step.block = name
        return cls(name=name, step_ids=[s.id for s in steps], **kwargs)

    def linear_connections(self, style: ConnectionStyle | None = None) -> list['Connection']:
        return [Connection(source=a, target=b, style=style)
                for a, b in zip(self.step_ids[:-1], self.step_ids[1:])]

@dataclass(slots=True)
class Connection:
    source: str
    target: str
    style: ConnectionStyle | None = None
    curved: bool = False
    
    @classmethod
    def from_blocks(cls, block1: Block, block2: Block, steps: list[Step],
                    style: ConnectionStyle | None = None) -> list['Connection']:
        step_map = {s.id: s for s in steps}
        b1 = [step_map[sid] for sid in block1.step_ids]
        b2 = [step_map[sid] for sid in block2.step_ids]
        b2_rc_map = {s.rc: s for s in b2}
        b1_rcs = {s.rc for s in b1}
        connections = []
        for s1 in b1:
            prev_rc = s1.rc - 1
            if prev_rc in b2_rc_map and prev_rc not in b1_rcs:
                connections.append(cls(source=s1.id,
                                       target=b2_rc_map[prev_rc].id,
                                       style=style))
            next_rc = s1.rc + 1
            if next_rc in b2_rc_map and next_rc not in b1_rcs:
                connections.append(cls(source=s1.id,
                                       target=b2_rc_map[next_rc].id,
                                       style=style))
        return connections

@dataclass
class Mechanism:
    steps: list[Step]
    connections: list[Connection] = field(default_factory=list)
    blocks: list[Block] = field(default_factory=list)
    auto_connect: bool = False
    def __post_init__(self):
        self.step_map = {step.id: step for step in self.steps}
        self.block_map = {block.name: block for block in self.blocks}

    def get_step(self, step_id: str) -> Step:
        return self.step_map[step_id]

    def get_block(self, block_name: str) -> Block:
        return self.block_map[block_name]

    def add_step(self, step: Step | list[Step]):
        if isinstance(step, list):
            self.step_map.update({s.id: s for s in step})
        else:
            self.step_map.update({step.id: step})
    
    def add_block(self, block: Block | list[Block]):
        if isinstance(block, list):
            self.block_map.update({b.name: b for b in block})
        else:
            self.block_map.update({block.name: block})

    def add_connection(self, connection: Connection | list[Connection]):
        if isinstance(connection, list):
            self.connections.extend(connection)
        else:
            self.connections.append(connection)
    def pop_connection(self, pairs):
        if not isinstance(pairs, list):
            pairs = [pairs]
        pairs = {tuple(p) for pair in pairs for p in (pair, pair[::-1])}
        self.connections = [c for c in self.connections
                            if (c.source, c.target) not in pairs]
            
        
    def compute_energies(self, df: pd.DataFrame,
                         energy_col: str,
                         reference: str | None = None) -> pd.DataFrame:

        # Calculate Absolute
        records = []
        for step in self.steps:
            energy = sum(df.loc[df['code'] == code, energy_col].iloc[0]
                         for code in step.codes)
            records.append({'step_id': step.id,
                            'label': step.label,
                            'rc': step.rc,
                            'x': step.get_x(),
                            'energy': energy,
                            'block': step.block})
        df_plot = pd.DataFrame(records)

        # Calculate Relative
        if reference is None:
            reference_energy = df_plot['energy'].min()
        else:
            reference_energy = df_plot.loc[df_plot['step_id'] == reference, 'energy'].iloc[0]
        df_plot['relative_energy'] = (df_plot['energy'] - reference_energy)

        return df_plot
    
    @classmethod
    def from_blocks(cls, steps: list[Step], blocks: list[Block],
                    connections: list[Connection] | None = None):
        internal_connections = []
        for block in blocks:
            internal_connections.extend(block.linear_connections())
        if connections is not None:
            internal_connections.extend(connections)
        return cls(steps=steps, blocks=blocks,
                   connections=internal_connections)