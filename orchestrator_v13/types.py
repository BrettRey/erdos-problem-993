from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, Optional


class AuthorityMode(str, Enum):
    AUTHORITATIVE = "AUTHORITATIVE"
    NON_AUTHORITATIVE = "NON_AUTHORITATIVE"


class PartitionId(str, Enum):
    PI0 = "PI0"
    PI1 = "PI1"
    PI2 = "PI2"
    PI3 = "PI3"
    PI4 = "PI4"


PARTITION_ID_ORDER: dict[PartitionId, int] = {
    PartitionId.PI0: 0,
    PartitionId.PI1: 1,
    PartitionId.PI2: 2,
    PartitionId.PI3: 3,
    PartitionId.PI4: 4,
}


class BucketId(str, Enum):
    S1 = "S1"
    S2 = "S2"
    S3 = "S3"


class ModuleTag(str, Enum):
    ENV = "ENV"
    A2_LOCAL = "A2_LOCAL"
    HARDSTEP = "HARDSTEP"


class ClosureStatus(str, Enum):
    CLOSED = "CLOSED"
    OPEN = "OPEN"
    UNKNOWN = "UNKNOWN"
    NO_DATA = "NO_DATA"


class RegimeCore(str, Enum):
    NO_BREAK = "NO_BREAK"
    CROSS_LINE_MISMATCH = "CROSS_LINE_MISMATCH"
    INTERNAL_LINE_DEFICIT = "INTERNAL_LINE_DEFICIT"


class ActionType(str, Enum):
    P = "P"
    M_A2 = "M_A2"
    B = "B"
    A = "A"
    C = "C"


ACTION_TYPE_ORDER: dict[ActionType, int] = {
    ActionType.P: 0,
    ActionType.M_A2: 1,
    ActionType.B: 2,
    ActionType.A: 3,
    ActionType.C: 4,
}


@dataclass(frozen=True)
class WitnessRaw:
    witness_id: str
    n: int
    a: int
    b: int
    k: int
    m: int
    D: float
    sum_err: float
    C10: float
    C01: float
    C11: float
    sum_all: float
    Lambda: float
    Gamma_L: float


@dataclass(frozen=True)
class Witness:
    raw: WitnessRaw
    R_shift: float
    Defect: float
    Gamma_G: float
    Alpha_req: float
    Gamma_M_star: float
    Gamma_D: float


@dataclass(frozen=True)
class Config:
    N_target: int
    N_authoritative: int
    epsilon_buffer: float
    eps: float
    alpha_19: float
    alpha_star: float
    Delta_alpha: float
    rho_max: float
    lambda_max: float
    K_G: int
    K_A: int
    K_L: int
    kappa_P: float
    kappa_A: float
    kappa_B: float
    kappa_C: float
    kappa_M: float
    K0: float = 1000.0


@dataclass(frozen=True, order=True)
class ClassKey:
    a: int
    b: int


@dataclass(frozen=True)
class ClassAgg:
    key: ClassKey
    alpha_c: float
    lambda_c: float
    local_c: float
    margin_c_star: float
    drift_c: float
    w_alpha_id: str
    w_lambda_id: str
    w_local_id: str


@dataclass(frozen=True)
class LineAgg:
    a: int
    alpha_line: float
    lambda_line: float
    local_line: float
    Delta_line: float
    c_alpha_line: ClassKey
    c_lambda_line: ClassKey
    c_local_line: ClassKey


@dataclass(frozen=True)
class GlobalAgg:
    alpha_front: float
    lambda_front: float
    local_front: float
    gap: float
    Delta_global: float
    w_alpha_id: str
    w_lambda_id: str
    w_local_id: str
    c_alpha: ClassKey
    c_lambda: ClassKey
    c_local: ClassKey


@dataclass(frozen=True)
class Diagnostics:
    regime_core: RegimeCore
    a_bad: Optional[int]
    Delta_line_max: float
    Delta_withinclass_max: float


@dataclass(frozen=True)
class Aggregates:
    global_: GlobalAgg
    classes: dict[ClassKey, ClassAgg]
    lines: dict[int, LineAgg]
    diagnostics: Diagnostics


@dataclass(frozen=True)
class PartitionSpec:
    id: PartitionId
    pred_S1: Callable[[Witness], bool]
    pred_S2: Callable[[Witness], bool]
    module_S1: ModuleTag
    module_S2: ModuleTag
    module_S3: ModuleTag
    pred_S1_str: str
    pred_S2_str: str
    pred_S3_str: str


@dataclass(frozen=True)
class Bucket:
    id: BucketId
    module: ModuleTag
    records: tuple[Witness, ...]


@dataclass(frozen=True)
class A2Calibration:
    bucket_id: BucketId
    rho_min: float
    rho_max: float
    lambda_max: float
    rho_hat: float
    lambda_hat: float
    active_witnesses: tuple[str, ...]
    feasible_domain: bool


@dataclass(frozen=True)
class HardstepEntry:
    witness_id: str
    a: int
    b: int
    n: int
    k: int
    Gamma_G: float
    Alpha_req: float
    Gamma_L: float
    flags: tuple[str, ...]


@dataclass(frozen=True)
class HardstepList:
    bucket_id: BucketId
    must_include: tuple[str, ...]
    entries: tuple[HardstepEntry, ...]


@dataclass(frozen=True)
class BucketStats:
    bucket_id: BucketId
    module: ModuleTag
    size: int
    alpha_s: float
    lambda_s: float
    local_s: float
    deficit_s: float
    w_alpha_s: Optional[str]
    w_lambda_s: Optional[str]
    w_local_s: Optional[str]


@dataclass(frozen=True)
class PartitionEval:
    partition: PartitionSpec
    buckets: tuple[Bucket, Bucket, Bucket]
    stats: tuple[BucketStats, BucketStats, BucketStats]
    a2: Optional[A2Calibration]
    hardstep: Optional[HardstepList]
    Phi: float
    Complexity: int
    bottleneck_bucket_ids: tuple[BucketId, ...]


@dataclass(frozen=True)
class RunMeta:
    version: str
    run_id: str
    N_target: int
    N_authoritative: int
    N_auth_used: int
    authority_mode: AuthorityMode


@dataclass(frozen=True)
class FrontierState:
    run_meta: RunMeta
    global_frontier: GlobalAgg
    line_summary: tuple[LineAgg, ...]
    diagnostics: Diagnostics
    provisional: Optional[dict[str, object]]


@dataclass(frozen=True)
class SelectionKey:
    Phi: float
    Complexity: int
    PartitionIdOrder: int


@dataclass(frozen=True)
class PartitionPlan:
    run_meta: RunMeta
    selected_partition_id: PartitionId
    selection_key: SelectionKey
    buckets: tuple[BucketStats, BucketStats, BucketStats]
    a2_calibration: Optional[A2Calibration]
    hardstep_summary: Optional[dict[str, object]]


@dataclass(frozen=True)
class RepairTarget:
    partition_id: Optional[PartitionId] = None
    bucket_id: Optional[BucketId] = None
    class_key: Optional[ClassKey] = None


@dataclass(frozen=True)
class RepairQueueItem:
    rank: int
    action_type: ActionType
    target: RepairTarget
    bucket_id: Optional[BucketId]
    estimated_phi_reduction: float
    estimated_cost: float
    priority: float
    witness_justification: tuple[str, ...]
    notes: Optional[str]


@dataclass(frozen=True)
class BucketObligation:
    bucket_id: BucketId
    module: ModuleTag
    status: ClosureStatus
    deficit: float
    constants: Optional[dict[str, float]]
    frontier_witnesses: Optional[dict[str, str]]
    a2_params: Optional[dict[str, float]]
    active_witnesses: Optional[tuple[str, ...]]
    hardstep_list: Optional[tuple[HardstepEntry, ...]]
    must_include: Optional[tuple[str, ...]]


@dataclass(frozen=True)
class GlobalClosure:
    status: ClosureStatus
    Phi: float
    epsilon_buffer: float
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class V13Obligations:
    run_meta: RunMeta
    partition_id: PartitionId
    buckets: tuple[BucketObligation, BucketObligation, BucketObligation]
    global_closure: GlobalClosure


@dataclass(frozen=True)
class RunOutputs:
    frontier_state: FrontierState
    partition_plan: PartitionPlan
    repair_queue: tuple[RepairQueueItem, ...]
    v13_obligations: V13Obligations


@dataclass(frozen=True)
class VerifierError:
    code: str
    path: str
    expected: object
    got: object
